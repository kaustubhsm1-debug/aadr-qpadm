"""Validate the qpAdm implementation against data where the truth is known.

Everything here is simulated from an explicit admixture graph -- no reference
panel, no downloads, no private data. Clone the repository and run it:

    python simulation/simulate_qpadm.py

Populations drift from a common root under the Balding-Nichols model: an allele
at frequency p in the parent is drawn in the child from Beta(p(1-F)/F,
(1-p)(1-F)/F), which has mean p and variance F*p(1-p). F is therefore the
expected drift on that edge. Individuals are then sampled binomially, so the
observed frequencies carry sampling error on top of drift, as real ones do.

The graph gives the outgroups the differential relatedness qpAdm needs: O3, O4
and O5 each share an internal node with one source, so the f4 statistics can
tell the three sources apart. Without that the model is unidentifiable and the
weights would be arbitrary -- a property worth seeing in a simulation before
trusting it on real data.

                      root
         +-------+------+------+------+
         |       |      |      |      |
        O1      O2     nA     nB     nC
                      /  \\   /  \\   /  \\
                    O3   S1 O4  S2 O5   S3

    target = w1*S1 + w2*S2 + w3*S3, followed by its own post-admixture drift.

What it measures
  1. bias            mean(estimate - truth) over replicates
  2. SE calibration  mean jackknife SE / empirical SD of the estimate
  3. CI coverage     fraction of 95% intervals containing the truth
  4. fit test size   fraction of p < 0.05 when the model is CORRECT (want ~5%)
  5. negative control  p-values when a source is replaced by an unrelated
                       population (want them tiny)
  6. the two-target difference estimator: bias, SE calibration, and type-I
     error when the two targets are in truth identical
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import qpadm                                                    # noqa: E402

NSNP = 60000
NBLOCK = 120
NIND = 20            # diploid individuals sampled per population
REPS = 40
TRUE_W = np.array([0.55, 0.25, 0.20])
SOURCES = ['S1', 'S2', 'S3']
RIGHTS = ['O1', 'O2', 'O3', 'O4', 'O5']


def drift(p, F, rng):
    """One edge of Balding-Nichols drift. F=0 returns p unchanged."""
    if F <= 0:
        return p.copy()
    a = p * (1 - F) / F
    b = (1 - p) * (1 - F) / F
    out = rng.beta(np.maximum(a, 1e-6), np.maximum(b, 1e-6))
    return np.clip(out, 0.0, 1.0)


def sample(p, nind, rng):
    """Observed frequency in nind diploids."""
    return rng.binomial(2 * nind, p) / (2.0 * nind)


def simulate(rng, w=TRUE_W, nsnp=NSNP, nind=NIND, extra_target=None):
    """Return a dict of observed allele frequencies keyed by population."""
    root = rng.uniform(0.05, 0.95, nsnp)
    nA, nB, nC = (drift(root, f, rng) for f in (0.05, 0.06, 0.07))
    true = {'O1': drift(root, 0.30, rng), 'O2': drift(root, 0.25, rng),
            'O3': drift(nA, 0.15, rng), 'O4': drift(nB, 0.14, rng),
            'O5': drift(nC, 0.16, rng),
            'S1': drift(nA, 0.10, rng), 'S2': drift(nB, 0.12, rng),
            'S3': drift(nC, 0.11, rng),
            # an unrelated population, for the negative control
            'X': drift(root, 0.28, rng)}
    mix = w[0] * true['S1'] + w[1] * true['S2'] + w[2] * true['S3']
    true['T'] = drift(mix, 0.02, rng)
    # T0 is the placebo: the SAME true weights and the same sources, but its own
    # independent post-admixture drift and its own sampled individuals. Any
    # difference T0 - T recovers is noise by construction. (Drawing it from a
    # freshly simulated graph instead would compare mixtures of different source
    # realisations and is not a null at all -- that mistake shows up as a ~70%
    # false-positive rate, which is how it was caught here.)
    true['T0'] = drift(mix, 0.02, rng)
    if extra_target is not None:
        mix2 = (extra_target[0] * true['S1'] + extra_target[1] * true['S2']
                + extra_target[2] * true['S3'])
        true['T2'] = drift(mix2, 0.02, rng)
    return {k: sample(v, nind, rng) for k, v in true.items()}


def blocks(nsnp=NSNP, nblock=NBLOCK):
    return np.repeat(np.arange(nblock), int(np.ceil(nsnp / nblock)))[:nsnp]


def main():
    rng = np.random.default_rng(20260413)
    blk = blocks()
    W, SE, P, PNEG = [], [], [], []
    D, DSE, DNULL, DNULLSE = [], [], [], []
    for r in range(REPS):
        pops = simulate(rng, extra_target=np.array([0.45, 0.35, 0.20]))
        Q = qpadm.Qpadm(pops, blk)
        f = Q.fit('T', SOURCES, RIGHTS)
        W.append(f['weights']); SE.append(f['se']); P.append(f['p'])
        # negative control: swap a true source for the unrelated population
        PNEG.append(Q.fit('T', ['S1', 'X', 'S3'], RIGHTS)['p'])
        # two-target difference, true difference in w2 = 0.35 - 0.25 = +0.10
        d = Q.fit_diff('T2', 'T', SOURCES, RIGHTS)
        D.append(d['diff'][1]); DSE.append(d['se'][1])
        # ... and the placebo: two targets with identical true weights
        d0 = Q.fit_diff('T0', 'T', SOURCES, RIGHTS)
        DNULL.append(d0['diff'][1]); DNULLSE.append(d0['se'][1])
        if (r + 1) % 10 == 0:
            print('  %d/%d replicates' % (r + 1, REPS), flush=True)

    W, SE, P = np.array(W), np.array(SE), np.array(P)
    D, DSE = np.array(D), np.array(DSE)
    DNULL, DNULLSE = np.array(DNULL), np.array(DNULLSE)
    cov = np.mean(np.abs(W - TRUE_W) < 1.96 * SE, axis=0)

    print('\nSingle-target fit, %d replicates, %d SNPs, %d blocks, %d individuals/population'
          % (REPS, NSNP, NBLOCK, NIND))
    print('%-10s %8s %9s %9s %9s %9s' % ('component', 'truth', 'mean', 'bias', 'SE/SD', 'cover95'))
    for i, nm in enumerate(SOURCES):
        print('%-10s %8.3f %9.4f %+9.4f %9.2f %9.2f'
              % (nm, TRUE_W[i], W[:, i].mean(), W[:, i].mean() - TRUE_W[i],
                 SE[:, i].mean() / W[:, i].std(ddof=1), cov[i]))
    print('\nFit test when the model is CORRECT:  %.1f%% of replicates reject at p<0.05'
          ' (nominal 5%%)' % (100 * np.mean(P < 0.05)))
    print('Negative control (unrelated source):  median p = %.3g, %d/%d rejected'
          % (np.median(PNEG), int(np.sum(np.array(PNEG) < 0.05)), REPS))
    print('\nTwo-target difference estimator (fit_diff)')
    print('  true difference +0.100:  mean %+.4f  bias %+.4f  SE/SD %.2f  cover95 %.2f'
          % (D.mean(), D.mean() - 0.10, DSE.mean() / D.std(ddof=1),
             np.mean(np.abs(D - 0.10) < 1.96 * DSE)))
    print('  true difference  0.000:  mean %+.4f  bias %+.4f  SE/SD %.2f  '
          'false positives %.1f%% (nominal 5%%)'
          % (DNULL.mean(), DNULL.mean(), DNULLSE.mean() / DNULL.std(ddof=1),
             100 * np.mean(np.abs(DNULL / DNULLSE) > 1.96)))
    print('\nRead SE/SD as: >1 conservative (error bars wider than the estimate actually '
          'varies), <1 optimistic.')


if __name__ == '__main__':
    main()
