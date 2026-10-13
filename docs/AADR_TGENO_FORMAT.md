# Reading AADR `.geno` files: the transposed (`TGENO`) layout

The Allen Ancient DNA Resource ships its genotypes in EIGENSTRAT *packed*
format. Most documentation and most example code describes the **SNP-major**
`GENO` layout. Recent AADR releases are **individual-major** `TGENO`, and the
two are byte-for-byte plausible as each other: a `TGENO` file decoded as `GENO`
produces a complete, well-formed, entirely wrong genotype matrix with no error
anywhere in the pipeline.

This cost two full extractions here. The notes below are what it takes to not
repeat that.

## Layout

```
offset 0      48-byte ASCII header, beginning "TGENO" (or "GENO")
offset 48     one fixed-size record per INDIVIDUAL  (TGENO)
                                   ... per SNP      (GENO)
```

Record size is the same rule in both cases, applied to a different axis:

```
TGENO:  reclen = max(48, ceil(nsnp / 4))     records = nind
GENO:   reclen = max(48, ceil(nind / 4))     records = nsnp
```

Within a record, genotypes are packed two bits each, **most significant bits
first**: element 0 occupies bits 7–6 of byte 0, element 1 bits 5–4, and so on.
Note this is the opposite order from PLINK `.bed`, which fills from the low bits
up — getting these two the same way round is a second, independent way to
produce a silently wrong matrix.

Values are counts of the **first** allele in column 5 of the `.snp` file, with
`3` meaning missing:

| packed value | meaning |
|---|---|
| 0 | homozygous for the `.snp` column-6 allele (a2) |
| 1 | heterozygous |
| 2 | homozygous for the `.snp` column-5 allele (a1) |
| 3 | missing |

## Confirm the layout arithmetically before decoding anything

Both AADR sets satisfy the `TGENO` arithmetic exactly, and neither satisfies the
`GENO` arithmetic:

```
1240K:  48 + 23089 x 308254 = 7,117,276,654 bytes   (nsnp = 1,233,013)
HO:     48 + 27594 x 146033 = 4,029,634,650 bytes   (nsnp =   584,131)
```

Check the file size against both formulas before writing a decoder. If one of
them matches to the byte, that is your answer; if neither does, stop.

## How a wrong decode announces itself — and how it does not

The failed extractions were structurally perfect: correct size, no exceptions,
plausible allele frequencies, populations separating sensibly in a PCA. What
gave them away was **biology, not software**:

- **Uniform missingness.** Every population, ancient and modern, came out at
  47.7–47.8 % missing. Real ancient-DNA capture missingness varies enormously
  between samples and sites; a constant is the signature of reading a bit
  pattern that has nothing to do with the data.
- **Zero heterozygotes** in 20,000 x 934 genotypes. Not "few" — none.

Neither would be caught by a unit test on the decoder. Check the distributional
properties of what comes out, per population, against what the biology requires.

## Verification gates worth keeping

These are cheap and each one has caught a real error:

1. **Publisher checksums** on every downloaded metadata file.
2. **Byte-identical re-fetch** of a random subset of records, decoded by an
   independent code path (here `numpy.unpackbits` against the hand-rolled
   shift-and-mask), reproducing every column.
3. **Cross-dataset concordance.** Overlap the panel with an independent source
   covering some of the same individuals — SGDP works well. Correct orientation
   gave 99.988 % agreement over 3.36M genotype comparisons; the reversed
   orientation gave 24.7 %. A test that cannot fail is not a test, and this one
   fails loudly.
4. **Allele polarity.** Check `.snp` column 5 against a reference genome at
   sites you have independent calls for. In AADR v66.p1 it is the GRCh37
   reference allele at 99.71 % of coinciding sites; exclude the rest rather than
   assuming.
5. **Heterozygosity sanity.** Diploid whole genomes in the archive sit near a
   25.5 % median. A panel of diploids that disagrees with that is wrong.

## Liftover

AADR 1240K coordinates are hg19/GRCh37. Lifting to GRCh38 with the UCSC chain
file maps 99.91 % of sites uniquely. Keep the unmapped and multi-mapping sites
out rather than approximating them, and keep the mapping itself cached: it is
deterministic and re-deriving it invites inconsistency between runs.

One consequence worth stating because it bit this project: if jackknife blocks
are derived from the *pre*-liftover coordinates and the panel is keyed on
*post*-liftover coordinates, the blocks are still contiguous and still valid,
but they are not the blocks a reader recomputing from the shipped keys will get.
Derive block identifiers from one coordinate system, and record which.
