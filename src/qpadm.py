"""qpAdm-style admixture-proportion estimation with block-jackknife inference.

Model: the target's allele frequency is a mixture of K source frequencies,
    p_T = sum_j w_j p_Sj,   sum_j w_j = 1.
Because the weights sum to one, for any pair of right (outgroup) populations
(R_1, R_i) and any fixed reference population X,
    f4(R_1, R_i; T, X) = sum_j w_j f4(R_1, R_i; S_j, X).
Taking X = S_K (the last source) gives the estimating equations
    y_i = f4(R_1, R_{i+1}; T, S_K),        i = 1 .. n-1
    A_ij = f4(R_1, R_{i+1}; S_j, S_K),     j = 1 .. K-1
    y = A w',   w_K = 1 - sum(w')
so there are n-1 equations for K-1 free weights, and the model is testable
whenever n-1 > K-1, i.e. more right populations than sources.

Inference is by delete-one-block jackknife over genomic blocks: the covariance
of y is estimated from the block-deleted replicates, the weights are obtained by
generalised least squares against that covariance, and the weight standard
errors come from refitting the whole GLS problem on each block-deleted
replicate. The model-fit statistic is the GLS residual quadratic form, compared
to chi-square with (n-1)-(K-1) degrees of freedom.

Two deliberate choices, both stated rather than hidden:
  * A single SNP set is used for every statistic in a model -- the intersection
    of non-missing sites across target, sources and rights (ADMIXTOOLS
    'useallsnps: NO'). Per-statistic SNP sets would use more data but make the
    statistics non-comparable and can bias the fit when ancients differ in
    coverage, which here they do by a factor of ten.
  * Weights are estimated UNCONSTRAINED. A weight outside [0,1] is diagnostic
    of a wrong source set and is reported as such, not clipped out of sight.
  * A is treated as fixed in the GLS. ADMIXTOOLS' likelihood accounts for
    sampling error in A as well; ignoring it would be expected to make these
    standard errors optimistic. Measured on simulated data (sec 3.3 of the
    report) they come out ~35% WIDER than the empirical spread of the estimate,
    i.e. conservative overall -- the block jackknife more than absorbs the
    omission. The fit test, separately, over-rejects about twofold.
"""
import numpy as np


def _block_sums(t, blk, nb):
    """Per-block sums of per-SNP terms, for fast leave-one-block-out means."""
    return np.bincount(blk, weights=t, minlength=nb)


class Qpadm:
    def __init__(self, freqs, blk):
        """freqs: name -> (n_snp,) alt-allele frequency, NaN = missing.
           blk:   (n_snp,) integer genomic block id."""
        self.freqs = freqs
        self.blk = np.asarray(blk)

    def _setup(self, target, sources, rights):
        names = [target] + list(sources) + list(rights)
        P = [self.freqs[n] for n in names]
        m = np.ones(len(P[0]), bool)
        for p in P:
            m &= np.isfinite(p)
        p = {n: self.freqs[n][m] for n in names}
        blk = self.blk[m]
        ub, blk_i = np.unique(blk, return_inverse=True)
        nb = len(ub)
        R1 = rights[0]; SK = sources[-1]
        K = len(sources); n = len(rights)
        yT = np.empty((n - 1, nb)); AT = np.empty((n - 1, K - 1, nb))
        for i in range(n - 1):
            d = p[R1] - p[rights[i + 1]]
            yT[i] = _block_sums(d * (p[target] - p[SK]), blk_i, nb)
            for j in range(K - 1):
                AT[i, j] = _block_sums(d * (p[sources[j]] - p[SK]), blk_i, nb)
        cnt = np.bincount(blk_i, minlength=nb).astype(float)
        return yT, AT, cnt, int(m.sum()), nb

    @staticmethod
    def _gls(y, A, Sinv):
        AtS = A.T @ Sinv
        w = np.linalg.solve(AtS @ A, AtS @ y)
        r = y - A @ w
        return w, float(r @ Sinv @ r)

    def fit_diff(self, t1, t2, sources, rights, S_fixed=None):
        """Weight difference between two targets, jackknifed JOINTLY.

        Both targets are fitted on the same SNP set (the intersection including
        both) and the same block-deleted replicates, so the errors -- which are
        shared through the sources, the right populations and the SNPs -- cancel
        in the difference instead of adding. Fitting separately and subtracting
        would use two independent error bars and overstate the uncertainty.
        """
        names = [t1, t2] + list(sources) + list(rights)
        m = np.ones(len(self.freqs[names[0]]), bool)
        for n in names:
            m &= np.isfinite(self.freqs[n])
        sub = Qpadm({n: self.freqs[n][m] for n in names}, self.blk[m])
        r1, r2 = sub.fit(t1, sources, rights), sub.fit(t2, sources, rights)
        y1, A1_, c1, nsnp, nb = sub._setup(t1, sources, rights)
        y2, A2_, _, _, _ = sub._setup(t2, sources, rights)
        N = c1.sum()
        # ONE weighting matrix, used for BOTH the point estimate and every
        # replicate, and built symmetrically from the two targets' jackknife
        # covariances. Previously the point estimate came from each target's own
        # iterated covariance while the replicates used a single matrix taken
        # from t1 alone: the SE then described a different estimator from the
        # one reported, and swapping the argument order changed it.
        # ... and ITERATED to self-consistency exactly as fit() does, so that the
        # difference is a difference of the same estimator as the weights it
        # differences. Without the iteration the two disagreed by 0.005 and the
        # weight table in the report would not subtract to the reported value.
        A_full = A1_.sum(2) / N
        yj1 = (y1.sum(1)[:, None] - y1) / (N - c1)[None, :]
        yj2 = (y2.sum(1)[:, None] - y2) / (N - c1)[None, :]
        Aj = (A1_.sum(2)[:, :, None] - A1_) / (N - c1)[None, None, :]

        def _shared(e1, e2):
            k1 = np.atleast_2d(np.cov(e1, bias=True) * (nb - 1))
            k2 = np.atleast_2d(np.cov(e2, bias=True) * (nb - 1))
            return np.linalg.pinv(0.5 * (k1 + k2))

        # S_fixed: reuse a weighting matrix estimated elsewhere (e.g. genome-wide)
        # instead of re-estimating it from this SNP subset. Used to separate a real
        # per-chromosome difference from the noise of estimating a 7x7 covariance
        # off ~25 blocks -- see report sec 3.2.
        S = _shared(yj1, yj2) if S_fixed is None else S_fixed
        wf1 = self._gls(y1.sum(1) / N, A_full, S)[0]
        wf2 = self._gls(y2.sum(1) / N, A_full, S)[0]
        for _ in range(0 if S_fixed is not None else 4):
            S = _shared(yj1 - np.einsum('ijb,j->ib', Aj, wf1),
                        yj2 - np.einsum('ijb,j->ib', Aj, wf2))
            wf1 = self._gls(y1.sum(1) / N, A_full, S)[0]
            wf2 = self._gls(y2.sum(1) / N, A_full, S)[0]
        d = np.empty((nb, len(sources)))
        for b in range(nb):
            yb1 = (y1.sum(1) - y1[:, b]) / (N - c1[b])
            yb2 = (y2.sum(1) - y2[:, b]) / (N - c1[b])
            Ab = (A1_.sum(2) - A1_[:, :, b]) / (N - c1[b])
            w1b = self._gls(yb1, Ab, S)[0]; w2b = self._gls(yb2, Ab, S)[0]
            d[b] = np.append(w1b - w2b, (1 - w1b.sum()) - (1 - w2b.sum()))
        diff = np.append(wf1 - wf2, (1 - wf1.sum()) - (1 - wf2.sum()))
        se = np.sqrt((nb - 1) / nb * ((d - d.mean(0)) ** 2).sum(0))
        return dict(t1=t1, t2=t2, sources=list(sources), diff=diff, se=se,
                    z=diff / np.where(se > 0, se, np.nan), n_snps=nsnp, n_blocks=nb,
                    w1=r1['weights'], w2=r2['weights'], p1=r1['p'], p2=r2['p'],
                    # the weights the difference is actually a difference OF.
                    # w1/w2 above are the separate single-target fits reported
                    # per population; these two subtract exactly to `diff`.
                    wj1=np.append(wf1, 1 - wf1.sum()),
                    wj2=np.append(wf2, 1 - wf2.sum()), S=S)

    def fit(self, target, sources, rights):
        yT, AT, cnt, nsnp, nb = self._setup(target, sources, rights)
        N = cnt.sum()
        y = yT.sum(1) / N
        A = AT.sum(2) / N
        # leave-one-block-out replicates
        yj = (yT.sum(1)[:, None] - yT) / (N - cnt)[None, :]
        Aj = (AT.sum(2)[:, :, None] - AT) / (N - cnt)[None, None, :]
        g = nb
        cov = np.cov(yj, bias=True) * (g - 1)        # jackknife covariance of y
        cov = np.atleast_2d(cov)
        Sinv = np.linalg.pinv(cov)
        w1, chi2 = self._gls(y, A, Sinv)
        # Iterate to a self-consistent weighting: the residual covariance must be
        # evaluated at the weights that minimise the quadratic form it defines,
        # otherwise the statistic is not at its own minimum and over-rejects.
        for _ in range(4):
            e_b = yj - np.einsum('ijb,j->ib', Aj, w1)
            cov_i = np.atleast_2d(np.cov(e_b, bias=True) * (g - 1))
            Sinv = np.linalg.pinv(cov_i)
            w1, chi2 = self._gls(y, A, Sinv)
        # The chi2 p-value below assumes df = (n_rights-1)-(n_sources-1), which is
        # only right if the residual covariance has full rank: pinv would quietly
        # drop directions and leave the df too large. Report it rather than assume.
        cov_rank = int(np.linalg.matrix_rank(cov_i))
        cov_cond = float(np.linalg.cond(cov_i))
        wj = np.empty((g, len(w1)))
        for b in range(g):
            wj[b] = self._gls(yj[:, b], Aj[:, :, b], Sinv)[0]
        se1 = np.sqrt((g - 1) / g * ((wj - wj.mean(0)) ** 2).sum(0))
        w = np.append(w1, 1 - w1.sum())
        se = np.append(se1, np.sqrt((g - 1) / g * (((1 - wj.sum(1)) -
                       (1 - wj.sum(1)).mean()) ** 2).sum()))
        df = (len(rights) - 1) - (len(sources) - 1)
        from scipy.stats import chi2 as chi2d
        pval = float(chi2d.sf(chi2, df)) if df > 0 else np.nan
        # Fit test on the RESIDUAL vector instead of on y alone. e_b = y_b - A_b w
        # uses each replicate's own A, so its jackknife covariance carries the
        # sampling error of the sources as well; the y-only covariance above
        # treats A as fixed and therefore over-rejects.
        eb = yj - np.einsum('ijb,j->ib', Aj, w1)
        e = y - A @ w1
        cov_e = np.atleast_2d(np.cov(eb, bias=True) * (g - 1))
        chi2_r = float(e @ np.linalg.pinv(cov_e) @ e)
        p_r = float(chi2d.sf(chi2_r, df)) if df > 0 else np.nan
        # Full block-deleted replicate weight matrix (n_blocks x K), so that any
        # DERIVED quantity -- a ratio of two components, a sum, a difference --
        # can be jackknifed directly instead of propagated by the delta method,
        # which would need the weight covariance and usually ignores it.
        w_reps = np.column_stack([wj, 1 - wj.sum(1)])
        return dict(target=target, sources=list(sources), weights=w, se=se,
                    chi2=chi2, df=df, p=pval, chi2_resid=chi2_r, p_resid=p_r,
                    n_snps=nsnp, n_blocks=nb, w_reps=w_reps,
                    cov_rank=cov_rank, cov_cond=cov_cond)

    @staticmethod
    def jack_derived(res, fn):
        """Jackknife any scalar function of the weight vector, e.g. a ratio.
        fn maps a length-K weight vector to a scalar."""
        reps = np.array([fn(w) for w in res['w_reps']])
        g = len(reps)
        return fn(res['weights']), np.sqrt((g - 1) / g * ((reps - reps.mean()) ** 2).sum())
