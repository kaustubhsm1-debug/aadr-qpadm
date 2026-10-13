options(warn=1, width=200)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({library(admixtools); library(dplyr)})
STEPPE_C <- c("Kazakhstan_MLBA", "Kazakhstan_MLBA_Alakul", "Kazakhstan_MLBA_Fedorovo",
              "Russia_KrasnoyarskKrai_MLBA", "Russia_KrasnoyarskKrai_MLBA_Andronovo", "Russia_Andronovo")
IPC <- "Iran_ShahriSokhta_BA2-2"
RIGHT <- c("Mbuti", "Papuan", "Han", "Iran_GanjDareh_N", "Georgia_KotiasKlde_Mesolithic", "Turkey_N", "Russia_Kostenki_UP", "Karitiana")
TARGETS <- strsplit(Sys.getenv("QPADM_TARGETS",
            "Pakistan_IA,VLR,STU,ITU,UBR,Brahmin,ILA,Kalash,Pathan"), ",")[[1]]
pops <- c(TARGETS, IPC, STEPPE_C, RIGHT, "ONG")
# primary cache diagnostics
f2 <- f2_from_precomp("ref/f2", verbose=FALSE)
cat("PRIMARY: dim f2 =", dim(f2), "\n")
m <- apply(f2, 1:2, function(x) weighted.mean(x, as.numeric(sub(".*l","",dimnames(f2)[[3]]))))
neg <- which(apply(f2, 1:2, min) < 0 & upper.tri(m), arr.ind=TRUE)
cat("pairs with any negative block f2:", nrow(neg), "\n")
negm <- which(m < 0 & upper.tri(m), arr.ind=TRUE)
cat("pairs with negative genome-wide mean f2:", nrow(negm), "\n")
if (nrow(negm)) print(data.frame(p1=rownames(m)[negm[,1]], p2=colnames(m)[negm[,2]], f2=m[negm]))
d <- data.frame(p1=rownames(m)[neg[,1]], p2=colnames(m)[neg[,2]], mean_f2=m[neg], nneg=apply(neg,1,function(i) sum(f2[i[1],i[2],]<0)))
print(d[order(-d$nneg),][1:min(25,nrow(d)),], row.names=FALSE)
# sensitivity: SNPs present in every population (maxmiss = 0)
fam <- read.table("merged.fam", stringsAsFactors=FALSE); fam <- fam[fam$V1 %in% pops,]
ip <- ifelse(fam$V1 %in% STEPPE_C, "Steppe_C", fam$V1)
extract_f2("merged", "ref/f2_maxmiss0", inds=fam$V2, pops=ip, blgsize=5e6, maxmiss=0, adjust_pseudohaploid=TRUE, overwrite=TRUE)
g <- f2_from_precomp("ref/f2_maxmiss0")
cat("MAXMISS0: dim f2 =", dim(g), "\n")
res <- lapply(TARGETS, function(t) {
  out <- qpadm(g, left=c(IPC,"Steppe_C","ONG"), right=RIGHT, target=t, verbose=FALSE)
  tibble(target=t, iranian=out$weights$weight[1], iranian_se=out$weights$se[1], steppe=out$weights$weight[2], steppe_se=out$weights$se[2],
         aasi=out$weights$weight[3], aasi_se=out$weights$se[3], p=out$popdrop$p[1], chisq=out$popdrop$chisq[1], dof=out$popdrop$dof[1], feasible=out$popdrop$feasible[1])
}) %>% bind_rows()
print(as.data.frame(res), digits=3)
write.csv(res, "threeway_admixtools2_maxmiss0.csv", row.names=FALSE)
cat("Steppe S-R:", round(res$steppe[2]-res$steppe[1],4), "\n")
