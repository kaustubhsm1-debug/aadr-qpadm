# Choosing populations from the AADR annotation file

Population selection is the step most likely to quietly determine your result,
and the least likely to be written down. These are the rules this repository
enforces in code, and why each exists.

## Select on structured fields, never by keyword

The `.anno` file has the information you need in specific columns:

| column | field |
|---|---|
| 14 | Group ID — the label everything downstream joins on |
| 15 | Locality |
| 16 | Political Entity |
| 10 | date (BP) |
| 26 | SNPs hit on the target set |

Searching Group IDs for substrings looks equivalent and is not. A keyword pass
over Group IDs missed, in one real attempt here: a 99-individual Iron Age
cohort, two 1000 Genomes cohorts of 99 and 103, every Gujarati and Punjabi
cohort, and most of the Steppe Middle-to-Late Bronze Age groups — while
**four outgroup labels matched nothing at all** and were silently dropped,
leaving a model that still ran and still produced confident numbers.

Specific traps, all encountered:

- **"Malta"** matches the Mediterranean island before it matches the Upper
  Palaeolithic Siberian site (`Russia_MA1_HG`).
- **"Iran"** as a substring also matches `Tirana` and `Kairan-1`.
- **Political Entity is where a sample was collected, not where its ancestry
  is from.** The Gujarati 1000 Genomes cohort (`GIH`) has Political Entity
  "USA" (Houston). Filtering South Asian cohorts by country loses it.
- Label spellings drift between AADR releases. `Russia_Ust_Ishim_HG_published`,
  `Russia_Kostenki14`, `Russia_MA1_HG_published` and
  `Italy_North_Villabruna_HG` were all valid in earlier releases and match
  nothing in v66.p1.

## Assert every label before fetching

The single most valuable line in this codebase:

```python
missing = sorted(w for w in wanted if w not in gidx)
if missing:
    raise SystemExit('labels absent from the panel: %s' % ', '.join(missing))
```

A label that matches nothing does not raise. It drops out of the model, which
then runs with fewer outgroups than you believe it has, and reports a p-value
computed against the wrong degrees of freedom. `scripts/run_qpadm.py` refuses to
start if any label in the config — including labels used only in a robustness
variant — is absent.

## Check coverage, not just presence

A population being in the archive does not make it usable. Median SNPs hit, per
group, decides this:

- `India_Harappan` (Rakhigarhi) is a single individual with 31,760 SNPs on the
  1240K set. It is the most *relevant* ancient genome for South Asian work and
  it cannot carry a model.
- `Turkmenistan_BA1-1` (Gonur) has 34 individuals at a median of 14,866 SNPs.
  Large n does not rescue per-individual coverage.
- By contrast `Uzbekistan_BA` (730,788) and `Iran_ShahTepe_BA` (941,334) are
  usable at the same date and region.

Pull the per-group inventory (n, date range, median SNPs hit) before choosing
anything, and record it alongside the results.

## Pooling

Several small related groups are often pooled into one source. Pool over
**individuals**, not over group frequencies:

```python
s = nansum(F[:, idx] * N[:, idx], axis=1)      # allele counts
n = N[:, idx].sum(1)                           # observed individuals
pooled = where(n >= min_n_pool, 1 - s / (2*n), nan)
```

Averaging the group frequencies instead would weight a group of 2 equally with
a group of 40.

Pooling is a modelling assumption, not a convenience: pooling populations the
data can distinguish is a model violation. Test it with an f4 statistic before
accepting it. In this project, pooling Western and Central Steppe looked
harmless and was not — f4(Mbuti, Han; WSteppe, CSteppe) = +0.000404, Z = +4.98,
and splitting them improved the positive control's fit substantially.

## Choosing a SNP set

The 1240K set is usually the right choice over Human Origins — roughly twice the
sites, and the ancient samples are captured on it. Larger "2M"-style sets can be
a trap: their additional sites may exist only on a newer capture reagent, so the
older ancient samples have no data there at all and the apparent gain is empty
columns.
