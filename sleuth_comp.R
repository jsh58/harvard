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

# compare full to SE
pval <- 0.05
merge.se <- merge(sleuth.res[, 1:2], sleuth.res.se[, 1:2], by='target_id')
merge.se <- subset(merge.se, ! is.na(merge.se$pval.x) & ! is.na(merge.se$pval.y))
# piecemeal, but cleaner than using table() in case of 0s
tp.se <- sum(merge.se$pval.x < pval & merge.se$pval.y < pval)
fp.se <- sum(merge.se$pval.x >= pval & merge.se$pval.y < pval)
tn.se <- sum(merge.se$pval.x >= pval & merge.se$pval.y >= pval)
fn.se <- sum(merge.se$pval.x < pval & merge.se$pval.y >= pval)
cat('SE', length.se, 'vs. full\n')
cat('  TPR:', ifelse(tp.se + fn.se > 0, tp.se/(tp.se + fn.se), 'n/a'), '\n')
cat('  FPR:', ifelse(fp.se + tn.se > 0, fp.se/(fp.se + tn.se), 'n/a'), '\n')
cat('  recall:', ifelse(tp.se + fn.se > 0, tp.se/(tp.se + fn.se), 'n/a'), '\n')
cat('  precision:', ifelse(tp.se + fp.se > 0, tp.se/(tp.se + fp.se), 'n/a'), '\n')

# compare full to PE
merge.pe <- merge(sleuth.res[, 1:2], sleuth.res.pe[, 1:2], by='target_id')
merge.pe <- subset(merge.pe, ! is.na(merge.pe$pval.x) & ! is.na(merge.pe$pval.y))
tp.pe <- sum(merge.pe$pval.x < pval & merge.pe$pval.y < pval)
fp.pe <- sum(merge.pe$pval.x >= pval & merge.pe$pval.y < pval)
tn.pe <- sum(merge.pe$pval.x >= pval & merge.pe$pval.y >= pval)
fn.pe <- sum(merge.pe$pval.x < pval & merge.pe$pval.y >= pval)
cat('PE', length.pe, 'vs. full\n')
cat('  TPR:', ifelse(tp.pe + fn.pe > 0, tp.pe/(tp.pe + fn.pe), 'n/a'), '\n')
cat('  FPR:', ifelse(fp.pe + tn.pe > 0, fp.pe/(fp.pe + tn.pe), 'n/a'), '\n')
cat('  recall:', ifelse(tp.pe + fn.pe > 0, tp.pe/(tp.pe + fn.pe), 'n/a'), '\n')
cat('  precision:', ifelse(tp.pe + fp.pe > 0, tp.pe/(tp.pe + fp.pe), 'n/a'), '\n')
