"""f-statistics with block jackknife standard errors.

Conventions (Patterson et al. 2012, Genetics 192:1065):
  f4(A,B;C,D) = mean_i (pA_i - pB_i)(pC_i - pD_i)
  f3(C; A,B)  = mean_i (pC_i - pA_i)(pC_i - pB_i)   [outgroup-f3 shared drift]
All frequencies are of the SAME allele at each SNP (orientation is the
caller's responsibility). SNPs with missing data in any population used by a
statistic are dropped for that statistic only.
"""
import numpy as np


class FreqPanel:
    """pops: dict name -> (n_snps,) float array of allele frequencies, NaN = missing.
       blk: (n_snps,) integer block id for the jackknife."""

    def __init__(self, pops, blk):
        self.pops = pops
        self.blk = np.asarray(blk)
        self.blocks = np.unique(self.blk)

    def _terms(self, a, b, c, d):
        pa, pb, pc, pd = (self.pops[x] for x in (a, b, c, d))
        m = ~(np.isnan(pa) | np.isnan(pb) | np.isnan(pc) | np.isnan(pd))
        return (pa[m] - pb[m]) * (pc[m] - pd[m]), self.blk[m]

    def _jack(self, t, blk):
        """delete-one-block jackknife on a mean of per-SNP terms."""
        n = len(t)
        tot = t.sum()
        est = tot / n
        vals, counts = [], []
        for b in self.blocks:
            sel = blk == b
            nb = int(sel.sum())
            if nb == 0 or nb == n:
                continue
            vals.append((tot - t[sel].sum()) / (n - nb))
            counts.append(nb)
        vals = np.array(vals); counts = np.array(counts)
        g = len(vals)
        h = n / counts                       # Busing et al. weights for unequal blocks
        pseudo = h * est - (h - 1) * vals
        var = np.sum((pseudo - pseudo.mean()) ** 2 / (h - 1)) / g
        return est, np.sqrt(var), n

    def f4(self, a, b, c, d):
        t, blk = self._terms(a, b, c, d)
        e, se, n = self._jack(t, blk)
        return dict(stat=f'f4({a},{b};{c},{d})', est=e, se=se, z=e / se if se > 0 else np.nan, n_snps=n)

    def f3out(self, c, a, b):
        pa, pb, pc = self.pops[a], self.pops[b], self.pops[c]
        m = ~(np.isnan(pa) | np.isnan(pb) | np.isnan(pc))
        t = (pc[m] - pa[m]) * (pc[m] - pb[m])
        e, se, n = self._jack(t, self.blk[m])
        return dict(stat=f'f3({c};{a},{b})', est=e, se=se, z=e / se if se > 0 else np.nan, n_snps=n)

    def f4_ratio(self, num, den):
        """ratio of two f4 statistics, jackknifed jointly (ratio of block-deleted means)."""
        tn, bn = self._terms(*num)
        td, bd = self._terms(*den)
        # NOTE: numerator and denominator are each computed on their OWN complete
        # SNP set -- _terms() masks independently and no intersection is taken.
        # The sets differ by the populations entering one statistic and not the
        # other. The jackknife below is still paired block-by-block, and blocks
        # missing from either statistic are dropped, but this is not a ratio
        # over a single common SNP set and the earlier comment claimed it was.
        est = tn.sum() / len(tn) / (td.sum() / len(td))
        vals, counts = [], []
        for b in self.blocks:
            sn, sd = bn == b, bd == b
            cn, cd = int(sn.sum()), int(sd.sum())
            if cn == 0 or cd == 0 or cn == len(tn) or cd == len(td):
                continue   # a block absent from either statistic makes the ratio undefined
            num_m = (tn.sum() - tn[sn].sum()) / (len(tn) - sn.sum())
            den_m = (td.sum() - td[sd].sum()) / (len(td) - sd.sum())
            vals.append(num_m / den_m); counts.append(sn.sum())
        vals = np.array(vals); counts = np.array(counts, float)
        n = len(tn); g = len(vals)
        h = n / np.clip(counts, 1, None)
        pseudo = h * est - (h - 1) * vals
        var = np.sum((pseudo - pseudo.mean()) ** 2 / (h - 1)) / g
        return dict(est=est, se=np.sqrt(var), n_snps=n)
