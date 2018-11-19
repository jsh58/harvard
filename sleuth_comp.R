# JMG 11/2018
# comparing sleuth results, full vs. PE vs. SE

library(sleuth)

# get CL args
args <- commandArgs(trailingOnly=T)
if (length(args) < 7) {
  cat('Need 7 CL args: samples (2), replicates, folder, lengths (3)\n')
  q()
}
sample.1 <- args[1]
sample.2 <- args[2]
reps <- args[3]
folder <- args[4]
length.full <- args[5]
length.se <- args[6]
length.pe <- args[7]

sample <- c(paste0(rep(paste0(sample.1, '_rep'), reps), 1:reps),
  paste0(rep(paste0(sample.2, '_rep'), reps), 1:reps))
condition <- c(rep(sample.1, reps), rep(sample.2, reps))
pref <- paste0(folder, '/', sample, '_')
suf <- '_trimmed_quant'

# full analysis
cat('Full:', length.full, '\n')
tab <- data.frame(sample=sample, condition=condition,
  path=paste0(pref, length.full, suf),
  stringsAsFactors=F)
so <- sleuth_prep(tab)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth.res <- sleuth_results(so, 'reduced:full', 'lrt')
#sleuth.sig <- subset(sleuth.res, qval <= 0.05)

# SE analysis
cat('SE:', length.se, '\n')
tab.se <- data.frame(sample=sample, condition=condition,
  path=paste0(pref, length.se, suf),
  stringsAsFactors=F)
so.se <- sleuth_prep(tab.se)
so.se <- sleuth_fit(so.se, ~condition, 'full')
so.se <- sleuth_fit(so.se, ~1, 'reduced')
so.se <- sleuth_lrt(so.se, 'reduced', 'full')
sleuth.res.se <- sleuth_results(so.se, 'reduced:full', 'lrt')
#sleuth.sig.se <- subset(sleuth.res.se, qval <= 0.05)

# PE analysis
cat('PE:', length.pe, '\n')
tab.pe <- data.frame(sample=sample, condition=condition,
  path=paste0(pref, length.pe, suf),
  stringsAsFactors=F)
so.pe <- sleuth_prep(tab.pe)
so.pe <- sleuth_fit(so.pe, ~condition, 'full')
so.pe <- sleuth_fit(so.pe, ~1, 'reduced')
so.pe <- sleuth_lrt(so.pe, 'reduced', 'full')
sleuth.res.pe <- sleuth_results(so.pe, 'reduced:full', 'lrt')
#sleuth.sig.pe <- subset(sleuth.res.pe, qval <= 0.05)

# compare full to SE
merge.se <- merge(sleuth.res[, 1:2], sleuth.res.se[, 1:2], by='target_id')
table(merge.se$pval.x < 0.05, merge.se$pval.y < 0.05)

# compare full to PE
merge.pe <- merge(sleuth.res[, 1:2], sleuth.res.pe[, 1:2], by='target_id')
table(merge.pe$pval.x < 0.05, merge.pe$pval.y < 0.05)
