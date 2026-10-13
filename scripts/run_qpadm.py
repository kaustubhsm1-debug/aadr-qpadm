"""Run a qpAdm model, or a set of robustness variants of it, from a JSON config.

    python scripts/run_qpadm.py config/example_swat.json

Replaces the one-off runner scripts that an analysis tends to accumulate. Every
threshold that changes the answer is named in the config and written into the
output, because the single most expensive failure in this kind of work is a
number you cannot reproduce six weeks later: a per-population minimum-count
threshold that lived only in a notebook session is enough to move a SNP count by
60,000 and an estimate by 0.01.

Config keys
  panel          directory holding group_freqs.npz, Fv.npy, Nc.npy (aadr_tfetch.py)
  sources        list of population labels, or {"name": [labels...]} to pool several
  rights         list of outgroup labels
  targets        list of population labels to fit
  min_n          non-missing individuals required per population per SNP  (default 2)
  min_n_pool     ... for pooled populations                               (default 4)
  block_mb       jackknife block size in Mb                               (default 5)
  contrasts      optional [[t1, t2], ...] pairs for the joint difference estimator
  variants       optional named overrides, each a partial config, run in addition
                 to the base model -- this is how a robustness table is produced
                 under ONE estimator instead of several ad hoc scripts
"""
import json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import qpadm                                                    # noqa: E402

DEFAULTS = dict(min_n=2, min_n_pool=4, block_mb=5, contrasts=[], variants={})


def load_panel(d):
    D = np.load(os.path.join(d, 'group_freqs.npz'), allow_pickle=True)
    Fv = np.load(os.path.join(d, 'Fv.npy'), mmap_mode='r')
    Nc = np.load(os.path.join(d, 'Nc.npy'), mmap_mode='r')
    return D, Fv, Nc, {g: i for i, g in enumerate(list(D['gnames']))}


def frequencies(cfg, D, Fv, Nc, gidx):
    """Allele frequency per population, NaN where too few individuals are observed.

    The packed AADR value counts the a1 (reference) allele, so the alternate
    frequency is 1 - value/2. Pooled sources are summed over individuals, not
    averaged over group frequencies: averaging would weight a group of 2 equally
    with a group of 40.
    """
    def one(label, min_n):
        f = 1 - np.asarray(Fv[:, gidx[label]], float) / 2.0
        f[np.asarray(Nc[:, gidx[label]]) < min_n] = np.nan
        return f

    pops = {}
    for spec in list(cfg['sources']) + list(cfg['rights']) + list(cfg['targets']):
        if isinstance(spec, dict):
            (name, members), = spec.items()
            idx = [gidx[m] for m in members]
            s = np.nansum(np.asarray(Fv[:, idx], float) * np.asarray(Nc[:, idx], float), axis=1)
            n = np.asarray(Nc[:, idx], float).sum(1)
            pops[name] = np.where(n >= cfg['min_n_pool'], 1 - s / np.maximum(2 * n, 1), np.nan)
        else:
            pops[spec] = one(spec, cfg['min_n'])
    return pops


def names(specs):
    return [list(s)[0] if isinstance(s, dict) else s for s in specs]


def blocks(key, mb):
    """Contiguous physical blocks for the jackknife. key packs (chrom<<32)|pos.

    Block boundaries are arbitrary -- any contiguous partition much larger than
    the LD scale works -- but they must be derived from ONE coordinate system
    and recorded, because they change the estimate slightly. The example config
    includes a 10 Mb variant so that sensitivity is visible rather than assumed.
    """
    return (key >> 32) * 1000000 + (key & 0xFFFFFFFF) // int(mb * 1e6)


def run_one(cfg, D, Fv, Nc, gidx, tag):
    pops = frequencies(cfg, D, Fv, Nc, gidx)
    Q = qpadm.Qpadm(pops, blocks(D['key38'], cfg['block_mb']))
    src, rgt = names(cfg['sources']), names(cfg['rights'])
    rows = []
    for t in names(cfg['targets']):
        try:
            r = Q.fit(t, src, rgt)
        except Exception as exc:                 # a target absent from the panel
            print('  skipped %s (%s)' % (t, type(exc).__name__)); continue
        row = dict(variant=tag, target=t, p=r['p'], chisq=r['chi2'], df=r['df'],
                   n_snps=r['n_snps'], n_blocks=r['n_blocks'],
                   cov_rank=r['cov_rank'], cov_cond=r['cov_cond'])
        for i, s in enumerate(src):
            row['w_' + s], row['se_' + s] = r['weights'][i], r['se'][i]
        rows.append(row)
    diffs = []
    for t1, t2 in cfg['contrasts']:
        d = Q.fit_diff(t1, t2, src, rgt)
        row = dict(variant=tag, contrast='%s - %s' % (t1, t2), n_snps=d['n_snps'])
        for i, s in enumerate(src):
            row['d_' + s], row['se_' + s], row['z_' + s] = d['diff'][i], d['se'][i], d['z'][i]
            row['wj1_' + s], row['wj2_' + s] = d['wj1'][i], d['wj2'][i]
        diffs.append(row)
    return rows, diffs


def main(path):
    cfg = dict(DEFAULTS, **json.load(open(path)))
    D, Fv, Nc, gidx = load_panel(cfg['panel'])
    wanted = set()
    for spec in cfg['sources'] + cfg['rights'] + cfg['targets']:
        wanted.update(list(spec.values())[0] if isinstance(spec, dict) else [spec])
    for over in cfg['variants'].values():
        for spec in over.get('sources', []) + over.get('rights', []) + over.get('targets', []):
            wanted.update(list(spec.values())[0] if isinstance(spec, dict) else [spec])
    missing = sorted(w for w in wanted if w not in gidx)
    if missing:
        # A label that matches nothing is silently dropped by every downstream
        # step and produces a plausible wrong answer. Refuse to start.
        raise SystemExit('labels absent from the panel: %s' % ', '.join(sorted(missing)))

    allrows, alldiffs = [], []
    for tag, over in [('base', {})] + list(cfg['variants'].items()):
        rows, diffs = run_one(dict(cfg, **over), D, Fv, Nc, gidx, tag)
        allrows += rows; alldiffs += diffs
        print('%-28s %d targets, %d contrasts' % (tag, len(rows), len(diffs)), flush=True)

    out = os.path.splitext(os.path.basename(path))[0]
    pd.DataFrame(allrows).to_csv(out + '_fits.csv', index=False)
    if alldiffs:
        pd.DataFrame(alldiffs).to_csv(out + '_contrasts.csv', index=False)
    print('wrote %s_fits.csv%s' % (out, ' and %s_contrasts.csv' % out if alldiffs else ''))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'config/example_swat.json')
