# JMG 11/2018
# comparing sleuth results, full vs. PE vs. SE
# plotting precision-recall

library(sleuth)

# get CL args
args <- commandArgs(trailingOnly=T)
if (length(args) < 8) {
  cat('Need 8 CL args: samples (2), replicates, folder, lengths (3), output\n')
  q()
}
sample.1 <- args[1]
sample.2 <- args[2]
reps <- args[3]
folder <- args[4]
length.full <- args[5]
length.se <- args[6]
length.pe <- args[7]
output <- args[8]

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

# merge full results with SE and PE
merge.se <- merge(sleuth.res[, 1:2], sleuth.res.se[, 1:2], by='target_id')
merge.se <- subset(merge.se, ! is.na(merge.se$pval.x) & ! is.na(merge.se$pval.y))
merge.pe <- merge(sleuth.res[, 1:2], sleuth.res.pe[, 1:2], by='target_id')
merge.pe <- subset(merge.pe, ! is.na(merge.pe$pval.x) & ! is.na(merge.pe$pval.y))

# create range of p-values
min.p <- min(merge.se[,2:3], merge.pe[,2:3], na.rm=T)
div <- 100
pvals <- numeric(div)
for (j in 1:div) {
  pvals[j] <- min.p^((div-j)/div)
}

# collect results over range of p-values
results.se <- data.frame()
results.pe <- data.frame()
for (pval in pvals) {
  tp.se <- sum(merge.se$pval.x <= pval & merge.se$pval.y <= pval)
  fp.se <- sum(merge.se$pval.x > pval & merge.se$pval.y <= pval)
  tn.se <- sum(merge.se$pval.x > pval & merge.se$pval.y > pval)
  fn.se <- sum(merge.se$pval.x <= pval & merge.se$pval.y > pval)
  recall.se <- ifelse(tp.se + fn.se > 0, tp.se/(tp.se + fn.se), NA)
  precision.se <- ifelse(tp.se + fp.se > 0, tp.se/(tp.se + fp.se), NA)
  results.se <- rbind(results.se, data.frame(recall=recall.se, precision=precision.se,
    pval=pval, length=paste0('SE_', length.se)))

  tp.pe <- sum(merge.pe$pval.x <= pval & merge.pe$pval.y <= pval)
  fp.pe <- sum(merge.pe$pval.x > pval & merge.pe$pval.y <= pval)
  tn.pe <- sum(merge.pe$pval.x > pval & merge.pe$pval.y > pval)
  fn.pe <- sum(merge.pe$pval.x <= pval & merge.pe$pval.y > pval)
  recall.pe <- ifelse(tp.pe + fn.pe > 0, tp.pe/(tp.pe + fn.pe), NA)
  precision.pe <- ifelse(tp.pe + fp.pe > 0, tp.pe/(tp.pe + fp.pe), NA)
  results.pe <- rbind(results.pe, data.frame(recall=recall.pe, precision=precision.pe,
    pval=pval, length=paste0('PE_', length.pe)))
}

# plot precision-recall
library(ggplot2)
results <- rbind(results.se, results.pe)
png(output, width=5.5, height=5, units="in", res=600, type='cairo')
ggplot(data=results, aes(x=recall, y=precision, color=length)) +
  geom_point(size=0.8) + geom_path() +
  xlab('Recall') + ylab('Precision') + xlim(0, 1) + ylim(0, 1)
dev.off()
