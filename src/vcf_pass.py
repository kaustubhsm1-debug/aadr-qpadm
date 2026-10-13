"""Single streaming pass over a VCF, including merged small-variant + structural-variant
output of the kind long-read pipelines produce (e.g. Clair3 + Sniffles2).

Per sample writes into work/<name>/:
  snv_auto.npz  key(int64=chrom*2^32+pos), gt(int8 1=het,2=homalt,3=other), dp,gq(int16), ra(uint8 ref*4+alt)
  snv_X.npz / snv_Y.npz  same, for chrX / chrY
  clinvar_hits.tsv  exact ClinVar P/LP allele matches (all chromosomes)
  cv_posonly.tsv    same position as a ClinVar P/LP variant but different allele
  chrM.tsv, sv.tsv
  stats.json
"""
import sys, os, gzip, json, pickle, collections
import numpy as np

BASES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
CHROM_IDX = {('chr%d' % i): i for i in range(1, 23)}
CHROM_IDX.update({'chrX': 23, 'chrY': 24, 'chrM': 25})
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}


def info_dict(s):
    d = {}
    for kv in s.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def run(name, path, clinvar_pkl, outdir):
    with open(clinvar_pkl, 'rb') as fh:
        cv = pickle.load(fh)
    cv_keep, cv_pos = cv['keep'], cv['pos']
    os.makedirs(outdir, exist_ok=True)

    arrs = {k: {'key': [], 'gt': [], 'dp': [], 'gq': [], 'ra': []} for k in ('auto', 'X', 'Y')}
    st = {'n_records': 0, 'by_type': collections.Counter(), 'by_chrom': collections.Counter(),
          'gt_counts': collections.Counter(), 'titv': collections.Counter(), 'phased': 0,
          'dp_hist': collections.Counter(), 'sv_types': collections.Counter(),
          'het_auto': 0, 'hom_auto': 0, 'X_het': 0, 'X_hom': 0, 'Y_n': 0, 'Y_dp': [],
          'auto_dp_sum': 0, 'auto_dp_n': 0, 'X_dp_sum': 0, 'X_dp_n': 0, 'multiallelic': 0}
    cv_hits, cv_posonly, chrM_rows, sv_rows = [], [], [], []

    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            chrom, pos_s, ref, alt, info, fmt, smp = p[0], p[1], p[3], p[4], p[7], p[8], p[9]
            ci = CHROM_IDX.get(chrom)
            if ci is None:
                continue
            pos = int(pos_s)
            st['n_records'] += 1
            st['by_chrom'][chrom] += 1
            g = dict(zip(fmt.split(':'), smp.split(':')))
            gtf = g.get('GT', './.')
            phased = '|' in gtf
            if phased:
                st['phased'] += 1
            a1a2 = gtf.replace('|', '/').split('/')
            try:
                h = sorted(int(x) for x in a1a2)
            except ValueError:
                h = None
            dp = int(g['DP']) if g.get('DP', '.').isdigit() else -1
            if dp < 0 and g.get('DR', '.').isdigit() and g.get('DV', '.').isdigit():
                dp = int(g['DR']) + int(g['DV'])
            gq = int(g['GQ']) if g.get('GQ', '.').lstrip('-').isdigit() else -1

            I = info_dict(info) if ('SVTYPE=' in info) else None
            is_sv = I is not None and 'SVTYPE' in I
            multi = ',' in alt
            if multi:
                st['multiallelic'] += 1

            # ---- classify
            if is_sv:
                st['by_type']['SV'] += 1
                st['sv_types'][I['SVTYPE']] += 1
                sv_rows.append('\t'.join([chrom, pos_s, I.get('SVTYPE', ''), str(I.get('SVLEN', '')),
                                          str(I.get('END', '')), gtf, str(gq), str(g.get('DR', '')),
                                          str(g.get('DV', '')), alt[:40] if alt.startswith('<') else ('seq%d' % len(alt))]))
            elif multi:
                st['by_type']['multiallelic'] += 1
            elif len(ref) == 1 and len(alt) == 1:
                st['by_type']['SNV'] += 1
                st['titv']['ti' if (ref, alt) in TRANSITIONS else 'tv'] += 1
            else:
                st['by_type']['indel'] += 1

            # ---- genotype bookkeeping
            if h is None:
                code = 3
            elif h == [0, 1]:
                code = 1
            elif h == [1, 1]:
                code = 2
            elif h == [0, 0]:
                code = 0
            else:
                code = 3

            if dp >= 0:
                st['dp_hist'][min(dp, 120)] += 1

            # ---- arrays (biallelic SNV only)
            if (not is_sv) and (not multi) and len(ref) == 1 and len(alt) == 1 and ref in BASES and alt in BASES:
                bucket = 'auto' if ci <= 22 else ('X' if ci == 23 else ('Y' if ci == 24 else None))
                if bucket:
                    a = arrs[bucket]
                    a['key'].append((ci << 32) | pos)
                    a['gt'].append(code)
                    a['dp'].append(min(dp, 32000))
                    a['gq'].append(min(max(gq, -1), 32000))
                    a['ra'].append(BASES[ref] * 4 + BASES[alt])
                    if bucket == 'auto':
                        st['het_auto'] += (code == 1)
                        st['hom_auto'] += (code == 2)
                        if dp > 0:
                            st['auto_dp_sum'] += dp
                            st['auto_dp_n'] += 1
                    elif bucket == 'X':
                        st['X_het'] += (code == 1)
                        st['X_hom'] += (code == 2)
                        if dp > 0:
                            st['X_dp_sum'] += dp
                            st['X_dp_n'] += 1
                    else:
                        st['Y_n'] += 1
                        if dp > 0:
                            st['Y_dp'].append(dp)

            if chrom == 'chrM':
                chrM_rows.append('\t'.join([pos_s, ref, alt, gtf, str(gq), str(dp),
                                            g.get('AF', ''), g.get('AD', '')]))

            # ---- ClinVar
            cvc = 'MT' if chrom == 'chrM' else chrom[3:]
            ps = cv_pos.get(cvc)
            if ps is not None and pos in ps and not is_sv:
                alts = alt.split(',')
                matched = False
                for ai, al in enumerate(alts, start=1):
                    rec = cv_keep.get((cvc, pos, ref, al))
                    if rec:
                        matched = True
                        carries = (h is not None and ai in h)
                        cv_hits.append('\t'.join([chrom, pos_s, ref, al, gtf, str(gq), str(dp),
                                                  g.get('AF', ''), g.get('AD', ''), str(int(carries))] + [str(x) for x in rec]))
                if not matched:
                    cv_posonly.append('\t'.join([chrom, pos_s, ref, alt, gtf, str(gq), str(dp), g.get('AD', '')]))

    for b in arrs:
        a = arrs[b]
        np.savez_compressed(os.path.join(outdir, 'snv_%s.npz' % b),
                            key=np.array(a['key'], dtype=np.int64), gt=np.array(a['gt'], dtype=np.int8),
                            dp=np.array(a['dp'], dtype=np.int16), gq=np.array(a['gq'], dtype=np.int16),
                            ra=np.array(a['ra'], dtype=np.uint8))
    hdr = {'clinvar_hits.tsv': 'chrom\tpos\tref\talt\tgt\tgq\tdp\taf\tad\tcarries\tgene\tclnsig\tclnrevstat\tclndn\talleleid\tmc\tclnvc\torigin\trs',
           'cv_posonly.tsv': 'chrom\tpos\tref\talt\tgt\tgq\tdp\tad',
           'chrM.tsv': 'pos\tref\talt\tgt\tgq\tdp\taf\tad',
           'sv.tsv': 'chrom\tpos\tsvtype\tsvlen\tend\tgt\tgq\tdr\tdv\talt'}
    for fn, rows in (('clinvar_hits.tsv', cv_hits), ('cv_posonly.tsv', cv_posonly),
                     ('chrM.tsv', chrM_rows), ('sv.tsv', sv_rows)):
        with open(os.path.join(outdir, fn), 'w') as fh:
            fh.write(hdr[fn] + '\n')
            fh.write('\n'.join(rows) + ('\n' if rows else ''))

    st['Y_dp_mean'] = float(np.mean(st['Y_dp'])) if st['Y_dp'] else 0.0
    st['Y_dp_median'] = float(np.median(st['Y_dp'])) if st['Y_dp'] else 0.0
    del st['Y_dp']
    st['sample'] = name
    st['n_clinvar_hits'] = len(cv_hits)
    out = {k: (dict(v) if isinstance(v, collections.Counter) else v) for k, v in st.items()}
    with open(os.path.join(outdir, 'stats.json'), 'w') as fh:
        json.dump(out, fh)
    return name, st['n_records'], len(cv_hits)


if __name__ == '__main__':
    print(run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
