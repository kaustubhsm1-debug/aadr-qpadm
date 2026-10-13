# Cross-check of the three-way ancestry proportions with ADMIXTOOLS 2
# (Maier, Flegontov, Flegontova, Isildak, Changmai & Reich 2023, eLife 12:e85492),
# the reference implementation of qpAdm (Haak et al. 2015, Nature 522:207, supp.).
#
# Input: merged.{bed,bim,fam}, written by src/export_plink.py --
#   the AADR 1240K set (hg19 coordinates) plus any study samples, one population
#   label per individual in the .fam population column. PLINK A1 is the AADR .snp
#   reference allele.
#
# Purpose: run the SAME models through the reference implementation so the
# in-house estimator (src/qpadm.py) can be checked against it. Targets are taken
# from the QPADM_TARGETS environment variable (comma-separated) so this script
# carries no study-specific sample names; the default is the published positive
# control plus modern reference cohorts. Run on any machine with a C++ toolchain:
#   install.packages("remotes"); remotes::install_github("uqrmaie1/admixtools")
#   Rscript admixtools_crosscheck.R

# Use the personal library if one is set (Ubuntu's /usr/local/lib/R/site-library
# is root-owned, so packages normally land in ~/R/library instead).
if (nzchar(Sys.getenv("R_LIBS_USER"))) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("admixtools", quietly = TRUE)) {
  cat("admixtools not found. .libPaths() currently:\n"); print(.libPaths())
  stop("Install it first: Rscript -e 'remotes::install_github(\"uqrmaie1/admixtools\", upgrade=\"never\")'")
}

library(admixtools)
library(dplyr)

PREFIX <- "merged"   # files sit next to this script in admixtools_check/
F2DIR  <- "ref/f2"

STEPPE_C <- c("Kazakhstan_MLBA", "Kazakhstan_MLBA_Alakul", "Kazakhstan_MLBA_Fedorovo",
              "Russia_KrasnoyarskKrai_MLBA", "Russia_KrasnoyarskKrai_MLBA_Andronovo",
              "Russia_Andronovo")
IPC    <- "Iran_ShahriSokhta_BA2-2"     # the AASI-bearing Shahr-i-Sokhta cluster
RIGHT  <- c("Mbuti", "Papuan", "Han", "Iran_GanjDareh_N",
            "Georgia_KotiasKlde_Mesolithic", "Turkey_N", "Russia_Kostenki_UP", "Karitiana")
TARGETS <- strsplit(Sys.getenv("QPADM_TARGETS",
            "Pakistan_IA,VLR,STU,ITU,UBR,Brahmin,ILA,Kalash,Pathan"), ",")[[1]]

# f2 blocks: 5 Mb to match the in-house jackknife; pooling the Central Steppe
# groups into one source population label.
pops <- c(TARGETS, IPC, STEPPE_C, RIGHT, "ONG")
# Central Steppe as a single source. qpadm() takes one label per source and does
# not pool a vector of labels, so the six groups are relabelled to one population
# ("Steppe_C") at the individual level when f2 is extracted -- which is what the
# in-house `pool()` does.
fam <- read.table(paste0(PREFIX, ".fam"), stringsAsFactors = FALSE)
fam <- fam[fam$V1 %in% pops, ]
ind_pops <- ifelse(fam$V1 %in% STEPPE_C, "Steppe_C", fam$V1)
if (!dir.exists(F2DIR)) {
  extract_f2(PREFIX, F2DIR, inds = fam$V2, pops = ind_pops, blgsize = 5e6,
             maxmiss = 0.5, adjust_pseudohaploid = TRUE, overwrite = TRUE)
}
f2 <- f2_from_precomp(F2DIR)

res <- lapply(TARGETS, function(t) {
  out <- qpadm(f2, left = c(IPC, "Steppe_C", "ONG"), right = RIGHT, target = t)
  tibble(target = t,
         iranian = out$weights$weight[1], iranian_se = out$weights$se[1],
         steppe  = out$weights$weight[2], steppe_se  = out$weights$se[2],
         aasi    = out$weights$weight[3], aasi_se    = out$weights$se[3],
         p = out$popdrop$p[1], chisq = out$popdrop$chisq[1], dof = out$popdrop$dof[1],
         feasible = out$popdrop$feasible[1])
}) %>% bind_rows()

print(as.data.frame(res), digits = 3)
write.csv(res, "threeway_admixtools2.csv", row.names = FALSE)

# Pairwise contrasts. ADMIXTOOLS 2 does not expose a paired-difference jackknife,
# so this reports each target's own fit; the joint-jackknife difference (the
# correct error bar, since two fits share sources, rights and SNPs) stays with
# qpadm.py::fit_diff, and only the point estimates are comparable here. Set
# QPADM_CONTRAST="A,B" to print one.
ct <- strsplit(Sys.getenv("QPADM_CONTRAST", ""), ",")[[1]]
if (length(ct) == 2) cat("\nSteppe,", ct[1], "-", ct[2], "(point estimates):",
    round(res$steppe[res$target == ct[1]] - res$steppe[res$target == ct[2]], 4), "\n")
