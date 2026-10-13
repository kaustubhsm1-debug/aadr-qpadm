# Validating a qpAdm implementation

An ancestry-proportion estimator will return three numbers that sum to one for
almost any input. Nothing about the output tells you whether it is right. These
are the checks that do, ordered by how much they are worth.

## 1. Simulation with a known answer

`simulation/simulate_qpadm.py` builds populations from an explicit admixture
graph with Balding-Nichols drift and binomial sampling, so the true mixture
weights are known. It needs no data and runs in about fifteen seconds:

```
python simulation/simulate_qpadm.py
```

Current output (40 replicates, 60,000 SNPs, 120 blocks, 20 diploids per
population, true weights 0.55 / 0.25 / 0.20):

```
component     truth      mean      bias     SE/SD   cover95
S1            0.550    0.5504   +0.0004      0.85      0.90
S2            0.250    0.2499   -0.0001      0.89      0.97
S3            0.200    0.1997   -0.0003      0.94      0.88

Fit test when the model is CORRECT:  7.5% of replicates reject at p<0.05 (nominal 5%)
Negative control (unrelated source):  median p = 1.96e-63, 40/40 rejected

Two-target difference estimator (fit_diff)
  true difference +0.100:  mean +0.0997  bias -0.0003  SE/SD 1.00  cover95 0.97
  true difference  0.000:  mean +0.0004  bias +0.0004  SE/SD 1.01  false positives 2.5%
```

Read it as follows. **Bias** is negligible at the fourth decimal, so the GLS is
solving the right problem. **SE/SD** compares the jackknife standard error with
how much the estimate actually varies across replicates: 0.85–0.94 on single
targets means the error bars are mildly *optimistic* there, which is why the
confidence intervals cover at 0.88–0.97 rather than 0.95. The **fit test
over-rejects** — 7.5 % where 5 % is nominal — so a borderline p-value is weaker
evidence against a model than it looks. The **negative control** is the check
that the fit test works at all: substitute an unrelated population for a true
source and every replicate is rejected at p ~ 1e-63. An estimator that cannot
reject a wrong model is not measuring anything.

The last two lines matter most if you are comparing two targets. The paired
difference estimator is unbiased, its standard error matches the empirical
spread exactly (SE/SD 1.00–1.01), and under a true null — two targets with
identical true weights, differing only by their own drift and sampling — it
produces a false positive 2.5 % of the time against a nominal 5 %.

**A simulation can only test what it also simulates.** This one cannot detect a
bias shared by the model and the generating process: if your real source
population is a proxy for an unsampled one, no amount of simulation will tell
you, because the simulation assumes the sources are the sources. That failure
mode needs the next two checks.

### The placebo has to be built correctly

The null above compares two targets drawn from the *same* sources with
independent post-admixture drift. An earlier version of this file instead drew
the second target from a freshly simulated graph — a mixture of different source
realisations entirely. That is not a null, and it reported a 70 % false-positive
rate. The symptom looked like a broken estimator; the cause was a broken
placebo. If your null test fails catastrophically, suspect the null first.

## 2. A target whose answer is already published

Simulation shares your assumptions. A published result does not. Fit a
population someone else has already decomposed with the same sources and see
whether you recover their proportions.

`config/example_swat.json` does this with the Swat Valley Iron Age cohort
(`Pakistan_IA`, n = 99) from Narasimhan et al. (2019). Expect Steppe ~0.25,
matching the published range, **and a rejected fit** (p ~ 0.005–0.01 depending
on block definition). That combination is the informative outcome: the
proportions reproduce, so the arithmetic is right, while the fit test says the
three-source model is incomplete for this target — which it is, because the
published Indus Periphery source is a subgroup that does not exist as a separate
population in the archive. Both halves are worth having. A control that only
ever passes teaches you nothing.

## 3. An independent implementation

Run the same models through ADMIXTOOLS 2 (Maier et al. 2023), the reference
implementation. `crosscheck/admixtools_crosscheck.R` does this from the PLINK
export; `crosscheck/sensitivity_maxmiss0.R` repeats it on the common-SNP
convention.

Two warnings from doing it here. First, compare on a **like-for-like SNP set**:
the two tools' default missing-data conventions differ, and that alone moved
individual weights by 2–2.7 standard errors while leaving differences between
targets almost unchanged. Second, `qpadm()` takes one label per source and will
not pool a vector of labels — pooled sources must be relabelled at the
individual level when `extract_f2` runs.

Expect agreement on weights to a fraction of a standard error. If a
disagreement survives the SNP-set alignment, the reference implementation is
right.

## 4. Pre-register the conditions, then honour them

Write down before you look: what effect size, what Z, what robustness
conditions, and **what you will say if they are not met**. Then apply that.

This is the check most easily failed, and it was failed here: of four
pre-registered conditions, one was never run, one was evaluated against a mask
that was never built, and one failed outright. The consequence had been written
in advance — label the result "suggestive, not established" — and it took an
external review to notice it had not been applied.

A related trap, worth naming because it is tempting and invisible: when a
pre-registered condition fails and you then find a better-powered diagnostic
that passes, you have not satisfied the condition. Report both, and say which
one you chose in advance.

## 5. Cheap checks that keep catching things

- **Assert every population label exists** before fetching or fitting. A label
  that matches nothing is silently dropped and the model runs with fewer
  outgroups than you think, against the wrong degrees of freedom.
- **Check the rank of the residual covariance.** The chi-square p-value assumes
  `df = (n_rights - 1) - (n_sources - 1)`; a pseudo-inverse of a rank-deficient
  covariance quietly invalidates that. `run_qpadm.py` reports rank and condition
  number beside every fit.
- **Vary the jackknife block size.** Block boundaries are arbitrary; if the
  answer moves materially between 2, 5 and 10 Mb, the error bar is not real.
- **Use a placebo pair** — two targets that should not differ — and check it
  returns null on real data, not just in simulation.
- **Derive figures and reports from the result tables**, never from numbers
  pasted into plotting code. Two external reviews of this project each caught a
  stale hand-typed number that had been superseded by a re-run.
