"""Extract selected individuals from an AADR packed *transposed* .geno (TGENO).

File layout, verified against the server's own byte count for both AADR datasets:
  48-byte ASCII header: 'TGENO <nind> <nsnp> <hash1> <hash2>'
  then one fixed-width record per INDIVIDUAL, of max(48, ceil(nsnp/4)) bytes,
  holding 2 bits per SNP, MSB-first within each byte, value 3 = missing.
    1240K: 48 + 23089 * 308254 = 7117276654  (matches Content-Range total)
    HO:    48 + 27594 * 146033 = 4029634650  (matches Content-Range total)

Individual-major layout means each wanted individual is ONE contiguous range
request, so only the kept individuals are transferred (288 MB of 7.1 GB for
934 individuals) and concurrent workers touch strictly disjoint output columns.
An md5 of every raw record is stored so a sequential re-fetch can be compared
byte-for-byte afterwards.

NB: the non-transposed 'GENO' format is SNP-major with record max(48,ceil(nind/4));
reading a TGENO file under that assumption yields structurally valid but entirely
wrong genotypes (no heterozygotes, uniform ~48% apparent missingness).
"""
import numpy as np, requests, sys, os, time, hashlib, json
from concurrent.futures import ThreadPoolExecutor

DATASETS = {
    '1240K': dict(geno='13994829', nind=23089, nsnp=1233013),
    'HO':    dict(geno='13994808', nind=27594, nsnp=584131),
}
BASE = 'https://dataverse.harvard.edu/api/access/datafile/'

# Which AADR populations to pull. Supplied as a file (one Group ID per line,
# '#' comments ignored) so the module carries no study-specific selection:
#     AADR_GROUPS=config/groups_example.txt  python src/aadr_tfetch.py ...
# See docs/PANEL_SELECTION.md for why these must be exact Group IDs taken from
# the .anno structured fields and never keyword-matched.
def load_groups(path=None):
    path = path or os.environ.get('AADR_GROUPS', 'config/groups_example.txt')
    with open(path) as fh:
        return [ln.split('#')[0].strip() for ln in fh
                if ln.split('#')[0].strip()]


KEEP_GROUPS = load_groups()


def decode_record(buf, nsnp):
    """2 bits per SNP, MSB-first; 3 -> -1. Independent of any column arithmetic."""
    bits = np.unpackbits(np.frombuffer(buf, np.uint8))
    v = (bits[0::2].astype(np.int8) << 1) | bits[1::2]
    v = v[:nsnp].copy()
    v[v == 3] = -1
    return v


def cache_fetch(which='1240K', workers=3, pause=0.4):
    """Fetch each kept individual's raw record to its own file under ref/aadr/rec_<ds>/.

    Resumable: a row whose cached file already has the exact record size is skipped.
    Workers touch strictly disjoint files, so the only shared resource is the server,
    which rate-limits (403 + 118-byte body) above roughly 5 requests/s -- hence the
    pause and the long backoff on 403 specifically.
    """
    workers = int(workers); pause = float(pause)
    cfg = DATASETS[which]
    nind, nsnp = cfg['nind'], cfg['nsnp']
    rec = max(48, (nsnp + 3) // 4)
    ind = [l.split() for l in open('ref/aadr/%s.ind' % which)]
    assert len(ind) == nind, 'ind file %d != header %d' % (len(ind), nind)
    groups = np.array([x[2] for x in ind]); ids = np.array([x[0] for x in ind])
    keep = np.where(np.isin(groups, KEEP_GROUPS))[0]
    cdir = 'ref/aadr/rec_%s' % which
    os.makedirs(cdir, exist_ok=True)
    todo = [int(i) for i in keep
            if not (os.path.exists('%s/%d.bin' % (cdir, i))
                    and os.path.getsize('%s/%d.bin' % (cdir, i)) == rec)]
    url = BASE + cfg['geno']
    print('%s: %d kept, %d cached, %d to fetch (%.0f MB), %d workers'
          % (which, len(keep), len(keep) - len(todo), len(todo),
             len(todo) * rec / 1e6, workers), flush=True)
    t0 = time.time(); done = [0]; fail = []

    def grab(i):
        off = 48 + i * rec
        sess = requests.Session()
        for attempt in range(6):
            try:
                r = sess.get(url, headers={'Range': 'bytes=%d-%d' % (off, off + rec - 1)},
                             timeout=600)
                if r.status_code == 206 and len(r.content) == rec:
                    tmp = '%s/%d.tmp' % (cdir, i)
                    open(tmp, 'wb').write(r.content)
                    os.replace(tmp, '%s/%d.bin' % (cdir, i))
                    done[0] += 1
                    if done[0] % 50 == 0:
                        print('  %d/%d  %.0f s' % (done[0], len(todo), time.time() - t0),
                              flush=True)
                    time.sleep(pause)
                    return
                time.sleep(30 * (attempt + 1) if r.status_code == 403 else 2 + attempt)
            except Exception:
                time.sleep(2 + attempt)
        fail.append(i)
        print('  FAILED %s (row %d)' % (ids[i], i), flush=True)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(grab, todo))
    print('fetched %d, failed %d, %.0f s' % (done[0], len(fail), time.time() - t0), flush=True)
    return fail


def assemble(which='1240K'):
    """Single-threaded, deterministic assembly of the cached records into one array."""
    cfg = DATASETS[which]
    nind, nsnp = cfg['nind'], cfg['nsnp']
    rec = max(48, (nsnp + 3) // 4)
    ind = [l.split() for l in open('ref/aadr/%s.ind' % which)]
    snp = [l.split() for l in open('ref/aadr/%s.snp' % which)]
    assert len(snp) == nsnp, 'snp file %d != header %d' % (len(snp), nsnp)
    groups = np.array([x[2] for x in ind]); ids = np.array([x[0] for x in ind])
    keep = np.where(np.isin(groups, KEEP_GROUPS))[0]
    cdir = 'ref/aadr/rec_%s' % which
    have = [j for j, i in enumerate(keep)
            if os.path.exists('%s/%d.bin' % (cdir, int(i)))
            and os.path.getsize('%s/%d.bin' % (cdir, int(i))) == rec]
    keep = keep[have]
    G = np.empty((nsnp, len(keep)), np.int8); md5 = {}
    for j, i in enumerate(keep):
        buf = open('%s/%d.bin' % (cdir, int(i)), 'rb').read()
        G[:, j] = decode_record(buf, nsnp)
        md5[str(int(i))] = hashlib.md5(buf).hexdigest()
    np.savez('ref/aadr/%s_T.npz' % which, G=G, ids=ids[keep], groups=groups[keep],
             rows=keep.astype(np.int64),
             snp_id=np.array([x[0] for x in snp]),
             snp_chrom=np.array([x[1] for x in snp]),
             snp_pos=np.array([int(float(x[3])) for x in snp], dtype=np.int64),
             snp_a1=np.array([x[4] for x in snp]), snp_a2=np.array([x[5] for x in snp]))
    json.dump(md5, open('ref/aadr/%s_T_md5.json' % which, 'w'))
    print('assembled %s from %d cached records' % (str(G.shape), len(keep)), flush=True)
    return G.shape


def main(which='1240K', workers=8):
    workers = int(workers)
    cfg = DATASETS[which]
    nind, nsnp = cfg['nind'], cfg['nsnp']
    rec = max(48, (nsnp + 3) // 4)
    ind = [l.split() for l in open('ref/aadr/%s.ind' % which)]
    assert len(ind) == nind, 'ind file %d != header %d' % (len(ind), nind)
    snp = [l.split() for l in open('ref/aadr/%s.snp' % which)]
    assert len(snp) == nsnp, 'snp file %d != header %d' % (len(snp), nsnp)
    groups = np.array([x[2] for x in ind]); ids = np.array([x[0] for x in ind])
    keep = np.where(np.isin(groups, KEEP_GROUPS))[0]
    url = BASE + cfg['geno']
    print('%s: %d ind, %d SNPs, record %d B; keeping %d ind / %d groups, %.0f MB'
          % (which, nind, nsnp, rec, len(keep), len(set(groups[keep])),
             len(keep) * rec / 1e6), flush=True)

    G = np.full((nsnp, len(keep)), -1, np.int8)
    md5 = [None] * len(keep)
    t0 = time.time(); done = [0]

    def grab(j):
        i = int(keep[j]); off = 48 + i * rec
        sess = requests.Session()
        for attempt in range(5):
            try:
                r = sess.get(url, headers={'Range': 'bytes=%d-%d' % (off, off + rec - 1)},
                             timeout=600)
                if r.status_code == 206 and len(r.content) == rec:
                    G[:, j] = decode_record(r.content, nsnp)
                    md5[j] = hashlib.md5(r.content).hexdigest()
                    done[0] += 1
                    if done[0] % 100 == 0:
                        print('  %d/%d  %.0f s' % (done[0], len(keep), time.time() - t0),
                              flush=True)
                    return
                print('  retry %s: status %s len %d' % (ids[i], r.status_code,
                      len(r.content)), flush=True)
            except Exception as e:
                print('  retry %s: %s' % (ids[i], type(e).__name__), flush=True)
            time.sleep(1 + attempt)
        print('  FAILED %s (row %d)' % (ids[i], i), flush=True)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(grab, range(len(keep))))
    nfail = sum(m is None for m in md5)
    print('failed records: %d' % nfail, flush=True)
    np.savez('ref/aadr/%s_T.npz' % which, G=G, ids=ids[keep], groups=groups[keep],
             rows=keep.astype(np.int64),
             snp_id=np.array([x[0] for x in snp]),
             snp_chrom=np.array([x[1] for x in snp]),
             snp_pos=np.array([int(float(x[3])) for x in snp], dtype=np.int64),
             snp_a1=np.array([x[4] for x in snp]), snp_a2=np.array([x[5] for x in snp]))
    json.dump({str(int(keep[j])): md5[j] for j in range(len(keep))},
              open('ref/aadr/%s_T_md5.json' % which, 'w'))
    print('saved %s in %.0f s' % (str(G.shape), time.time() - t0), flush=True)
    return nfail


if __name__ == '__main__':
    main(*(sys.argv[1:]))
