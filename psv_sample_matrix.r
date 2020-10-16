#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pdf(NULL)

suppressMessages(library(dplyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(viridis))
suppressMessages(library(circlize))

if (length(args) == 0) {
  data_dir <- '~/Mapping/hg38/jvc/genes/012/SMN1/paralog_ploidy/region1/'
  plots_prefix <- '~/Downloads/'
} else {
  data_dir <- args[1]
  dir_split <- unlist(strsplit(data_dir, '/', fixed=T))
  gene <- dir_split[1]
  region_n <- dir_split[length(dir_split)]
  plots_prefix <- sprintf('%s/%s.%s.', args[2], gene, region_n)
}
cat(sprintf('Reading data from %s/*.csv\n', data_dir))
cat(sprintf('Saving plots to %smatrix.png\n', plots_prefix))

load <- function(name, ...) {
  read.csv(sprintf('%s/%s.csv', data_dir, name), sep='\t', comment.char='#', ...)
}

# Loading genotype matrix

gt_df <- load('sample_matrix', stringsAsFactors=F)
gt_df <- subset(gt_df, n_psvs > 0)
gt_df <- gt_df[, colSums(gt_df == '*') != nrow(gt_df)]
n_samples <- nrow(gt_df)
n_psvs <- ncol(gt_df) - 2
if (n_samples == 0) {
  cat(sprintf('No appropriate samples and PSVs [%s]\n', data_dir))
  stop(1)
}

# Changing format and sorting genotype matrix

gt_matrix <- gt_df[, c(-1, -2)] %>% as.matrix
rownames(gt_matrix) <- gt_df$sample
colnames(gt_matrix) <- gsub('\\.', ':', colnames(gt_df)[c(-1, -2)])
gt_matrix <- apply(gt_matrix, 2, function(x) sapply(strsplit(x, ' ', fixed=T), `[`, 1))
gt_matrix[gt_matrix == '*'] <- NA
gt_matrix <- t(gt_matrix)

# Loading PSV weights

psv_weights_df <- load('psv_weights', stringsAsFactors=F)
psv_weights_df <- subset(psv_weights_df, iteration == max(psv_weights_df$iteration))
psv_weights <- psv_weights_df[match(rownames(gt_matrix), psv_weights_df$psv), 'weight']

# Loading full sample genotypes

sample_gt <- load('sample_weights', stringsAsFactors=F)
sample_gt <- subset(sample_gt, iteration == max(sample_gt$iteration))
sample_gt$iteration <- NULL
genotypes <- gsub('\\.', ',', gsub('X', '', colnames(sample_gt)[-1]))
sample_gt <- sample_gt[match(colnames(gt_matrix), sample_gt$sample),]
sample_gt$best <- apply(sample_gt[,-1], 1, which.max)
sample_gt$best_gt <- genotypes[sample_gt$best]
sample_gt$best_prob <-
  sapply(1:nrow(sample_gt), function(i) sample_gt[i, 1 + sample_gt$best[i]])
sample_gt$best_prob <- 10 ** sample_gt$best_prob
all_values <- sort(unique(c(sample_gt$best_gt, do.call(c, apply(gt_matrix, 1, unique)))))

# Clustering samples.

sample_dist <- matrix(apply(expand.grid(1:n_samples, 1:n_samples), 1,
  function(ij) sum(gt_matrix[, ij[1]] != gt_matrix[, ij[2]], na.rm=T)), nrow=n_samples)
rownames(sample_dist) <- colnames(gt_matrix)
colnames(sample_dist) <- colnames(gt_matrix)
sample_dist <- as.dist(sample_dist)
sample_clust <- hclust(sample_dist)

# Creating annotations

gt_colors <- structure(viridis(length(all_values)), names = all_values)
psv_colors <- colorRamp2(c(0, 0.5, 1), c('red', 'white', 'blue'))
psv_annot <- HeatmapAnnotation(
  Weight = psv_weights, col = list('Weight' = psv_colors),
  which = 'row', show_annotation_name=F)

sample_prob_colors <- colorRamp2(c(0.4, 0.7, 1), c('red', 'white', 'blue'))
sample_annot <- HeatmapAnnotation(
  Genotype = sample_gt$best_gt,
  Probability = sample_gt$best_prob,
  col = list(Genotype = gt_colors, Probability = sample_prob_colors),
  which = 'column', show_legend=c(Genotype=F))

scale <- 2
png(sprintf('%smatrix.png', plots_prefix),
    width=2000 * scale, height=800 * scale, res=100 * scale)
h = Heatmap(gt_matrix,
        name = 'Genotype',
        col = gt_colors,
        use_raster = TRUE,
        column_names_gp = grid::gpar(fontsize = 5.5),
        left_annotation = psv_annot,
        bottom_annotation = sample_annot,
        cluster_columns = sample_clust)
draw(h, merge_legend=T)
ignore <- invisible(dev.off())
cat(sprintf('Success [%s]\n', data_dir))

# Loading genotype matrix

prob_df <- load('sample_matrix2', stringsAsFactors=F)
prob_df <- subset(prob_df, n_psvs > 0)
prob_df <- prob_df[, colSums(is.na(prob_df)) != nrow(prob_df)]
n_samples <- nrow(prob_df)
n_psvs <- ncol(prob_df) - 4
if (n_samples == 0) {
  cat(sprintf('No appropriate samples and PSVs [%s]\n', data_dir))
  stop(1)
}

# Changing format and sorting genotype matrix

prob_matrix <- prob_df[, -1:-4] %>% as.matrix
rownames(prob_matrix) <- prob_df$sample
colnames(prob_matrix) <- gsub('\\.', ':', colnames(prob_df)[-1:-4])
prob_matrix <- t(prob_matrix)

col_fn <- colorRamp2(c(0, 0.5, 1), c('red', 'white', 'blue'))
h = Heatmap(10 ** prob_matrix,
            name = 'Genotype',
            col = viridis(100),
            use_raster = TRUE,
            column_names_gp = grid::gpar(fontsize = 5.5),
            left_annotation = psv_annot,
            bottom_annotation = sample_annot,
            cluster_columns = sample_clust,
            cluster_rows = F)
draw(h, merge_legend=T)
