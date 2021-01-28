suppressMessages(library(dplyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(viridis))
suppressMessages(library(circlize))
library(RColorBrewer)
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
ht_opt$message <- F

colors_to_ramp <- function(colors, min_v=NULL, max_v=NULL, data=NULL) {
  if (is.null(min_v)) {
    min_v <- min(data)
  }
  if (is.null(max_v)) {
    max_v <- max(data)
  }
  n <- length(colors)
  colorRamp2(seq(min_v, max_v, length = n), colors)
}

draw_matrices <- function(region_name, sample_gts, sample_psv_support, likelihoods,
                          f_values, out_prefix, important_samples=c()) {
  start_time <- proc.time()['elapsed']
  region_group = likelihoods$region_group[1]
  cat(sprintf('[%s: %s] Start\n',
              region_name, region_group))
  stopifnot(all(likelihoods$region_group == region_group))
  
  # Select best cluster and load its likelihood.
  last_iter_lik <- likelihoods %>% group_by(cluster) %>% slice_tail(n = 1) %>% ungroup
  best_cluster <- last_iter_lik[which.max(last_iter_lik$likelihood),]$cluster[1]
  likelihoods <- filter(likelihoods, cluster == best_cluster)
  last_iteration <- tail(likelihoods, 1)$iteration
  last_likelihood <- tail(likelihoods, 1)$likelihood

  # Extract best sample genotypes and probabilities.
  stopifnot(all(sample_gts$region_group == region_group))
  sample_gts <- filter(sample_gts, cluster == best_cluster & iteration == last_iteration)
  sample_gts <- sample_gts[, !is.na(sample_gts[1,])]
  samples <- colnames(sample_gts)[6:ncol(sample_gts)]
  n_samples <- length(samples)
  if (n_samples == 0) {
    cat(sprintf('[%s: %s] No samples present.\n', region_name, region_group))
    return()
  }
  all_gts <- sample_gts$genotype

  best_gts <- all_gts[apply(sample_gts[samples], 2, which.max)]
  names(best_gts) <- samples
  best_gt_probs <- apply(sample_gts[samples], 2, max)
  best_gt_probs <- 10 ^ best_gt_probs
  names(best_gt_probs) <- samples

  # Extract best PSV genotypes and probabilities.
  stopifnot(all(sample_psv_support$region_group == region_group))
  sample_psv_support <- select(sample_psv_support, 2:3, all_of(samples))
  sample_psv_support <- sample_psv_support[
    rowSums(is.na(sample_psv_support[samples])) < n_samples,]
  psvs <- unique(sample_psv_support$psv)
  n_psvs <- length(psvs)
  if (n_psvs == 0) {
    cat(sprintf('[%s: %s] No PSVs present.\n', region_name, region_group))
    return()
  }
  psv_support_long <- pivot_longer(sample_psv_support, all_of(samples),
                                   names_to='sample', values_to='prob')
  
  # Best PSV genotypes.
  psv_best_gts <- aggregate(prob ~ psv + sample, psv_support_long, FUN=which.max)
  psv_best_gts$gt <- all_gts[psv_best_gts$prob]
  psv_best_gts$prob <- NULL
  psv_best_gts <- pivot_wider(psv_best_gts, names_from='sample', values_from='gt')
  psv_best_gts <- psv_best_gts[match(psvs, psv_best_gts$psv),]
  psv_best_gts <- column_to_rownames(psv_best_gts, 'psv') %>% as.matrix
  psv_best_gts <- psv_best_gts[psvs, , drop=F][, samples, drop=F]
  
  # Support for the best sample genotypes.
  psv_support_best <- psv_support_long[
    psv_support_long$genotype == best_gts[psv_support_long$sample],]
  psv_support_best$genotype <- NULL
  psv_support_best <- pivot_wider(psv_support_best, names_from='sample', values_from='prob')
  psv_support_best <- column_to_rownames(psv_support_best, 'psv') %>% as.matrix()
  psv_support_best <- 10 ^ psv_support_best
  psv_support_best <- psv_support_best[psvs, , drop=F][, samples, drop=F]
  
  # Clustering samples.
  sample_dist <- matrix(apply(expand.grid(1:n_samples, 1:n_samples), 1,
      function(ij) sum(psv_best_gts[, ij[1]] != psv_best_gts[, ij[2]], na.rm=T)),
      nrow=n_samples)
  rownames(sample_dist) <- samples
  colnames(sample_dist) <- samples
  sample_dist <- as.dist(sample_dist)
  sample_clust <- hclust(sample_dist)
  
  # Loading PSV f values.
  stopifnot(all(f_values$region_group == region_group))
  f_values <- filter(f_values, cluster == best_cluster & iteration == last_iteration)
  f_values <- f_values[match(psvs, f_values$psv),]
  f_values <- f_values[, apply(f_values, 2, function(x) sum(is.na(x))) < n_psvs]
  rownames(f_values) <- NULL
  f_values <- column_to_rownames(f_values, 'psv')
  f_matrix <- select(f_values, starts_with('copy')) %>% as.matrix
  colnames(f_matrix)
  colnames(f_matrix) <- paste('Copy', 1:ncol(f_matrix))
  
  # Sample annotation.
  gt_colors <- structure(viridis(length(all_gts)), names = all_gts)
  zero_one_colors <- colors_to_ramp(brewer.pal(9, 'YlOrRd'), 0.0, 1.0)
  
  sample_annot <- HeatmapAnnotation(
    Genotype = best_gts,
    Weight = best_gt_probs,
    col = list(Genotype = gt_colors, Weight = zero_one_colors),
    which = 'column', show_legend=c(Genotype=F))
  
  sample_annot_short <- HeatmapAnnotation(
    Weight = best_gt_probs,
    col = list(Weight = zero_one_colors),
    which = 'column')
  
  # PSV annotation.
  psv_colors <- colorRamp2(c(0, 0.5, 1), c('red', 'white', 'blue'))
  psv_annot <- HeatmapAnnotation(
    Weight = f_matrix, col = list('Weight' = zero_one_colors),
    which = 'row', show_annotation_name=T)
  
  # Additional information.
  reliable_psvs <- which(rowSums(f_matrix < 0.8) == 0)
  n_reliable <- length(reliable_psvs)
  info <- sprintf('Cluster %d,  likelihood %.3f,  %s iterations,  %d/%d reliable PSVs.',
                  best_cluster, last_likelihood, last_iteration, n_reliable, n_psvs)

  # Genotype heatmap.
  sample_colors <- ifelse(samples %in% important_samples, 'red', 'black')
  scale <- 1.8
  png(sprintf('%sgenotypes.png', out_prefix),
      width=2000 * scale, height=1000 * scale, res=100 * scale)
  h = Heatmap(psv_best_gts,
              name = 'Genotype',
              col = gt_colors,
              use_raster = TRUE,
              column_names_gp = grid::gpar(fontsize = 5.5, col = sample_colors),
              left_annotation = psv_annot,
              bottom_annotation = sample_annot,
              cluster_columns = sample_clust,
              column_title = sprintf('%s: %s. Genotype matrix.\n%s',
                                     region_name, region_group, info))
  draw(h, merge_legend=T)
  ignore <- invisible(dev.off())
  
  # Probabilities heatmap.
  prob_col_fn <- colorRamp2(c(0, 0.5, 1), c('red', 'white', 'blue'))
  png(sprintf('%sprob_all.png', out_prefix),
      width=2000 * scale, height=1000 * scale, res=100 * scale)
  h = Heatmap(psv_support_best,
              name = 'Probability',
              col = prob_col_fn,
              use_raster = TRUE,
              column_names_gp = grid::gpar(fontsize = 5.5, col = sample_colors),
              left_annotation = psv_annot,
              bottom_annotation = sample_annot_short,
              cluster_columns = sample_clust,
              cluster_rows = F,
              column_title = sprintf('%s: %s. PSV-sample concordance.\n%s',
                                     region_name, region_group, info))
  draw(h, merge_legend=T)
  ignore <- invisible(dev.off())
  
  # Probabilities heatmap (only reliable PSVs).
  if (n_reliable > 0) {
    good_psv_annot <- HeatmapAnnotation(
      Weight = f_matrix[reliable_psvs, , drop=F], col = list('Weight' = zero_one_colors),
      which = 'row', show_annotation_name=T)
    
    png(sprintf('%sprob_reliable.png', out_prefix),
        width=2000 * scale, height=1000 * scale, res=100 * scale)
    h = Heatmap(psv_support_best[reliable_psvs, , drop=F],
                name = 'Probability',
                col = prob_col_fn,
                use_raster = TRUE,
                column_names_gp = grid::gpar(fontsize = 5.5, col = sample_colors),
                left_annotation = good_psv_annot,
                bottom_annotation = sample_annot_short,
                cluster_columns = sample_clust,
                cluster_rows = F,
                column_title = sprintf('%s: %s. PSV-sample concordance (only reliable).\n%s',
                                       region_name, region_group, info))
    draw(h, merge_legend=T)
    ignore <- invisible(dev.off())
  }
  
  cat(sprintf('[%s: %s] Success (%.1f seconds)\n',
              region_name, region_group, proc.time()['elapsed'] - start_time))
}
