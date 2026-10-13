"""Export a merged AADR + study-sample panel as PLINK 1 binary, for ADMIXTOOLS 2.

Allele convention, established empirically rather than assumed (docs/AADR_TGENO_FORMAT.md):

  * AADR .snp column 5 (a1) is the GRCh37 reference allele at 99.71 % of sites
    that coincide with a diploid variant call in our test data.
  * The AADR packed genotype value counts the a1 (reference) allele. Fixed by
    99.988 % concordance with SGDP over 3.36M genotype comparisons; the opposite
    orientation gives 24.7 %, so this is a real test and not a formality.
  * Declaring PLINK A1 := AADR a1 therefore makes the AADR value equal the PLINK
    A1-count directly, with no transformation.
  * For samples read from a VCF, PLINK A1-count = 2 - (ALT dosage), restricted to
    sites that lifted uniquely to GRCh38 and where a1 is the reference allele.

PLINK 1 .bed layout: magic 6c 1b 01, then SNP-major blocks of ceil(nind/4) bytes;
two bits per individual, low bits first within each byte (individual 0 in the low
2 bits), codes 00 = hom A1, 10 = het, 11 = hom A2, 01 = missing.

CAVEAT for anyone exporting VCF-derived samples: a VCF holding only non-reference
calls encodes every site with no record as homozygous reference. Without a gVCF or
BAM there is no way to distinguish "confidently reference" from "not covered", and
the sample is biased toward the reference allele wherever coverage failed. Nothing
in this file mitigates that -- it has to be handled upstream, with a callability
mask built from the alignments, or bounded downstream by a dropout simulation.
"""
import numpy as np

CODE = np.array([3, 2, 0], np.uint8)       # A1-count 0,1,2 -> PLINK 11, 10, 00
MISS = np.uint8(1)                          # PLINK 01


def write_bed(path, Gmat, chunk=20000):
    """Gmat: (n_snp, n_ind) int8 of A1 counts 0/1/2, negative = missing."""
    nsnp, nind = Gmat.shape
    rowbytes = (nind + 3) // 4
    pad = rowbytes * 4 - nind
    with open(path, 'wb') as fh:
        fh.write(bytes([0x6c, 0x1b, 0x01]))
        for i in range(0, nsnp, chunk):
            blk = Gmat[i:i + chunk]
            c = np.where((blk >= 0) & (blk <= 2), CODE[np.clip(blk, 0, 2)], MISS).astype(np.uint8)
            if pad:
                c = np.hstack([c, np.full((c.shape[0], pad), MISS, np.uint8)])
            c = c.reshape(c.shape[0], rowbytes, 4)
            packed = (c[:, :, 0] | (c[:, :, 1] << 2) | (c[:, :, 2] << 4) | (c[:, :, 3] << 6))
            fh.write(packed.astype(np.uint8).tobytes())
    return nsnp, nind


def write_bim(path, snp_id, chrom, pos, a1, a2):
    with open(path, 'w') as fh:
        for s, c, p, x, y in zip(snp_id, chrom, pos, a1, a2):
            fh.write('%s\t%s\t0\t%d\t%s\t%s\n' % (c, s, p, x, y))


def write_fam(path, iids, groups):
    with open(path, 'w') as fh:
        for i, g in zip(iids, groups):
            fh.write('%s\t%s\t0\t0\t0\t1\n' % (g, i))
