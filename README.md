# aadr-qpadm

A small, dependency-light qpAdm implementation in Python, the machinery to read
the Allen Ancient DNA Resource directly, and — the part most likely to save
someone a week — written-down notes on the ways this analysis fails silently.

Nothing here needs a C++ toolchain, and the validation suite needs no data at
all. It is not a replacement for [ADMIXTOOLS 2](https://github.com/uqrmaie1/admixtools),
which is the reference implementation and should be your arbiter; this exists
because having a second, independent implementation is how you find out that
one of them is wrong.

```
src/        qpadm.py, f4_tools.py, aadr_tfetch.py, export_plink.py, vcf_pass.py
scripts/    run_qpadm.py        config-driven model runner
simulation/ simulate_qpadm.py   self-contained validation, no data needed
crosscheck/ ADMIXTOOLS 2 scripts for checking results against the reference tool
config/     example model and population list
docs/       the things that cost time to learn
```

## Start here

```bash
pip install numpy scipy pandas requests
python simulation/simulate_qpadm.py
```

Fifteen seconds, no downloads. It builds populations from a known admixture
graph and checks that the estimator recovers the weights, that the error bars
mean what they claim, that the fit test rejects a wrong source, and that the
paired-difference estimator does not manufacture differences between targets
that are in truth identical. Output and interpretation: [docs/VALIDATION.md](docs/VALIDATION.md).

To fit real data you need the AADR, which `src/aadr_tfetch.py` will pull by
range request without downloading the whole 7 GB file. Then:

```bash
python scripts/run_qpadm.py config/example_swat.json
```

The example reproduces a published positive control — the Swat Valley Iron Age
cohort from Narasimhan et al. (2019) — and is set up to demonstrate the right
habit: a base model plus named robustness variants in one config, all run
through one estimator, every threshold recorded in the output.

## What is actually useful here

**`docs/AADR_TGENO_FORMAT.md`.** Recent AADR `.geno` files are individual-major
`TGENO`, not the SNP-major `GENO` that most documentation describes. The two are
byte-for-byte plausible as each other: decoding one as the other yields a
complete, well-formed, entirely wrong genotype matrix, with no error raised
anywhere. It cost two full extractions here before biology caught it — uniform
47.7 % missingness across every population, and zero heterozygotes in 20,000 ×
934 genotypes. The document gives the layout, the arithmetic that distinguishes
the formats before you decode anything, and five verification gates that each
caught something real.

**`docs/PANEL_SELECTION.md`.** Population selection decides results and is
rarely written down. Select on the annotation's structured fields, never by
keyword-matching Group IDs: a keyword pass here missed a 99-individual Iron Age
cohort and two 1000 Genomes cohorts, while four outgroup labels matched nothing
at all and were silently dropped from a model that still produced confident
numbers. Includes the specific traps — "Malta" matches the Mediterranean island
before the Siberian site, "Iran" also matches Tirana, and Political Entity
records where a sample was collected rather than where its ancestry is from, so
the Gujarati 1000 Genomes cohort is filed under USA.

**`docs/VALIDATION.md`.** Five checks ordered by what they are worth, with the
honest output of each, including the ones that did not pass.

## `qpadm.py`

Generalised least squares on the f4 matrix, following Haak et al. (2015) and
Harney et al. (2021). Differences from a naive implementation, each of which
changed an answer:

- The residual covariance is **iterated to self-consistency**. Evaluating it at
  weights that do not minimise the quadratic form it defines leaves the
  statistic off its own minimum, and the fit test over-rejects.
- Weights are **unconstrained**. A weight outside [0, 1] is diagnostic of a
  wrong source set and is reported as such, not clipped out of sight.
- The fit test runs on the **residual vector** with each replicate's own design
  matrix, so the covariance carries the sources' sampling error rather than
  treating them as fixed.
- `fit_diff(t1, t2, ...)` compares two targets **jointly**: same SNP
  intersection, same block-deleted replicates, and one weighting matrix built
  symmetrically from both targets' residual covariances and iterated the same
  way. Errors shared through the sources, outgroups and sites then cancel in
  the difference instead of adding. Fitting separately and subtracting
  overstates the uncertainty; using an asymmetric weighting makes the standard
  error depend on argument order, which is a bug this went through before it
  was caught.
- `jack_derived()` jackknifes any scalar function of the weights — a ratio, a
  sum, a difference — directly from the block-deleted replicates, instead of
  propagating by a delta method that would need the full weight covariance.
- Every fit reports the **rank and condition number** of its residual
  covariance, because the chi-square degrees of freedom assume full rank and a
  pseudo-inverse will quietly not tell you.

## Provenance and scope

The code was written for an analysis of four whole genomes against the AADR.
That analysis, its data and its results are not part of this repository and
will not be; what is here is the methodology, the format and selection notes,
and the validation suite, all of which are general.

Results here were cross-checked against ADMIXTOOLS 2 and agreed to a fraction
of a standard error on a like-for-like SNP set. Where the two disagree, believe
ADMIXTOOLS 2.

## References

- Haak W, et al. (2015). Massive migration from the steppe was a source for
  Indo-European languages in Europe. *Nature* 522:207–211. (qpAdm.)
- Harney É, Patterson N, Reich D, Wakeley J (2021). Assessing the performance of
  qpAdm. *Genetics* 217:iyaa045.
- Maier R, Flegontov P, Flegontova O, Işıldak U, Changmai P, Reich D (2023). On
  the limits of fitting complex models of population history to f-statistics.
  *eLife* 12:e85492. (ADMIXTOOLS 2.)
- Mallick S, et al. (2024). The Allen Ancient DNA Resource (AADR): a curated
  compendium of ancient human genomes. *Scientific Data* 11:182.
- Narasimhan VM, et al. (2019). The formation of human populations in South and
  Central Asia. *Science* 365:eaat7487.
- Busing FMTA, Meijer E, van der Leeden R (1999). Delete-m jackknife for unequal
  m. *Statistics and Computing* 9:3–8. (Weighted block jackknife.)
- Patterson N, et al. (2012). Ancient admixture in human history. *Genetics*
  192:1065–1093. (f-statistics.)

## Licence

MIT — see [LICENSE](LICENSE).
