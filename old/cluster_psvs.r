library(ggplot2)
library(tidyr)
library(dplyr)

ref_fractions <- function(column, split='|', remove_last=T) {
  sapply(column, function(s) {
    if (s == '*') {
      NA
    } else {
      values = as.numeric(unlist(strsplit(s, split, fixed=T)))
      if (remove_last) {
        values[1] / sum(head(values, -1))
      } else {
        values[1] / sum(values)
      }
    }
  })
}

psv_counts <- read.csv('~/Data/hg38/jvc/runs/021/SMN1/depth/psv_matrix.csv',
         sep='\t', comment.char='#')
psv_counts <- subset(psv_counts, pos >= 70915001 & pos <= 70956598)

sample_genotypes <- read.csv(
  '~/Data/hg38/jvc/plots/em_opt/data/SMN1/021/region1/sample_weights.csv',
  sep='\t', comment.char='#')
samples <- subset(sample_genotypes, iteration == max(sample_genotypes$iteration))$sample
n_samples <- length(samples)

exp_ref_fractions <- ref_fractions(psv_counts$exp_gt_counts, ',', F)

psv_fractions <- psv_counts[match(samples, colnames(psv_counts))] %>%
  select(starts_with('HG')) %>%
  apply(1, ref_fractions)
psv_fractions <- psv_fractions / exp_ref_fractions
colnames(psv_fractions) <- psv_counts$pos

psv_fractions <- psv_fractions[, colSums(!is.na(psv_fractions)) >= 100]
dim(psv_fractions)

psv_cor <- cor(psv_fractions, use = "pairwise.complete.obs")
psv_clustering <- hclust(as.dist(1 - abs(psv_cor)))
plot(psv_clustering)
clusters <- cutree(psv_clustering, k=4)

fractions_long <- psv_fractions %>% as.data.frame %>%
  tibble::rownames_to_column('sample') %>%
  pivot_longer(!sample, names_to='psv', values_to='ref_fraq')
fractions_long$cluster <- clusters[fractions_long$psv]

ggplot(fractions_long) +
  geom_boxplot(aes(psv, 2 * ref_fraq, fill=factor(cluster)))

