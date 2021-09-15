library(ggplot2)
library(viridis)
library(tidyverse)
library(ComplexHeatmap)

RELIABLE_THRESHOLD <- 0.95

transform_f_values <- function(f_values, likelihoods, min_samples=30) {
  f_values <- filter(f_values, !is.na(copy1) & n_samples >= min_samples)
  f_values$mean_fval <- select(f_values, starts_with('copy')) %>%
    apply(1, mean, na.rm=T)
  f_values$min_fval <- select(f_values, starts_with('copy')) %>%
    apply(1, min, na.rm=T)
  f_values$reliable <- f_values$min_fval >= RELIABLE_THRESHOLD
  f_values
}

load_f_values <- function(upper_dir, region_name, dataset_names) {
  all_f_values <- data.frame()
  for (i in seq_along(dataset_names)) {
    dir <- sprintf('%s/%s/%s/r009/extra', upper_dir, names(dataset_names)[i], region_name)
    if (!file.exists(sprintf('%s/em_likelihoods.csv', dir))) {
      cat(sprintf('WARN: Directory %s does not contain necessary files.\n', dir))
      next
    }
    f_values <- read.csv(sprintf('%s/em_f_values.csv', dir),
                            sep='\t', comment.char='#')
    f_values <- transform_f_values(f_values)
    if (nrow(f_values) == 0) {
      next
    }

    f_values$dataset <- dataset_names[[i]]
    all_f_values <- rbind(all_f_values, f_values)
  }
  all_f_values$dataset <- factor(all_f_values$dataset, levels=dataset_names)
  all_f_values
}

compare_f_values <- function(region_name, f_values, plots_dir) {
  region_group <- f_values$region_group[1]
  f_values <- f_values[, colSums(is.na(f_values)) < nrow(f_values)]
  datasets <- levels(f_values$dataset)

  stopifnot(all(f_values$region_group == region_group))
  f_values2 <- select(f_values, !'reliable') %>%
    pivot_longer('info_content' | starts_with('copy') | ends_with('fval'),
      names_to='stat', values_to='value')

  f_values2$stat <- as.character(f_values2$stat)
  f_values2 <- mutate(f_values2,
                         stat = case_when(stat == 'info_content' ~ 'Information content',
                                          stat == 'mean_fval' ~ 'Mean f-value',
                                          stat == 'min_fval' ~ 'Minimal f-value',
                                          stat == 'copy1' ~ 'Copy 1',
                                          stat == 'copy2' ~ 'Copy 2',
                                          stat == 'copy3' ~ 'Copy 3',
                                          stat == 'copy4' ~ 'Copy 4',
                                          stat == 'copy5' ~ 'Copy 5',
                                          T ~ stat))
  
  f_values3 <- select(f_values2, !'n_samples') %>%
    pivot_wider(names_from='dataset', values_from='value')

  dataset_dist <- data.frame()
  for (stat in unique(f_values3$stat)) {
    curr_f_values <- f_values3[f_values3$stat == stat,]
    curr_dist <- expand.grid(dataset1=datasets, dataset2=datasets, stat=stat)
    curr_dist$dist <- apply(curr_dist, 1,
      function(d2) mean(abs(curr_f_values[[d2[1]]] - curr_f_values[[d2[2]]]), na.rm=T))
    dataset_dist <- rbind(dataset_dist, curr_dist)
  }

  dist_matrix <- filter(dataset_dist, stat == 'Minimal f-value') %>%
    pivot_wider(names_from='dataset2', values_from='dist') %>%
    column_to_rownames('dataset1') %>%
    select(-1) %>% as.matrix
  if (nrow(dist_matrix) > 2 & all(!is.na(dist_matrix))) {
    clustering <- hclust(as.dist(dist_matrix))
    order <- rownames(dist_matrix)[clustering$order]
    dataset_dist$dataset1 <- factor(dataset_dist$dataset1, levels=order)
    dataset_dist$dataset2 <- factor(dataset_dist$dataset2, levels=order)
  }

  ggplot(dataset_dist) +
    geom_tile(aes(dataset1, dataset2, fill=dist)) +
    facet_wrap(~ stat) +
    scale_fill_viridis('Distance') +
    scale_x_discrete(limits = rev(levels(dataset_dist$dataset2))) +
    theme_bw() +
    theme(axis.title=element_blank()) +
    ggtitle(sprintf('%s: %s', region_name, region_group))
  ggsave(sprintf('%s/%s.%s.distance.png', plots_dir, region_name, region_group),
         width=10, height=6, scale=1.3)
  
  datasets_present <- length(unique(f_values$dataset))
  reliable_psvs <- lapply(datasets,
                          function(dat) subset(f_values, dataset == dat & reliable)$psv)
  reliable_all <- names(which(
    apply(list_to_matrix(reliable_psvs), 1, sum) == datasets_present))

  copy1_f_values3 <- filter(f_values3, stat == 'Copy 1')
  present_in_all <- copy1_f_values3[apply(is.na(copy1_f_values3), 1, sum) == 0, ]$psv
  unreliable_any <- names(which(
    apply(list_to_matrix(reliable_psvs), 1, sum) != datasets_present))
  unreliable_any <- intersect(unreliable_any, present_in_all)

  ggplot(filter(f_values2, startsWith(stat, 'Copy'))) +
    geom_vline(xintercept = reliable_all, size=2, color='blue', alpha=.1) +
    geom_vline(xintercept = unreliable_any, size=2, color='red', alpha=.1) +
    geom_point(aes(factor(psv), value, color=dataset),
               position=position_dodge(width=0.5)) +
    geom_hline(yintercept = RELIABLE_THRESHOLD) +
    scale_x_discrete('PSV position', labels=function(x) gsub('[^:]*:', '', x)) +
    scale_y_continuous('f-value') +
    scale_color_discrete('Dataset') +
    facet_wrap(~ stat, ncol=1) +
    theme_bw() +
    ggtitle(sprintf('%s: %s', region_name, region_group)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8))
  ggsave(sprintf('%s/%s.%s.fval.png', plots_dir, region_name, region_group),
         width=10, height=6)
  
  if (sum(sapply(reliable_psvs, length) > 0) < 2) {
    return()
  }
  
  scale <- 1
  png(sprintf('%s/%s.%s.upset.png', plots_dir, region_name, region_group),
      width=2000 * scale, height=1000 * scale, res=200 * scale)
  names(reliable_psvs) <- datasets
  reliable_psvs_comb <- make_comb_mat(reliable_psvs)
  ht <- draw(UpSet(reliable_psvs_comb,
                   column_title = sprintf('%s: %s. Reliable PSVs',
                                          region_name, region_group)))
  od <- column_order(ht)
  cs <- comb_size(reliable_psvs_comb)
  decorate_annotation("intersection_size", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
  })
  ignore <- invisible(dev.off())
}

dataset_names <- c('241.EAS' = 'East-Asian',
                   '242.EUR' = 'European',
                   '243.AFR' = 'African',
                   '244.AMR' = 'Admixed-American',
                   '245.SAS' = 'South-Asian')

upper_dir <- '~/Data/hg38/jvc/runs'
plots_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/f_val'


# for (path in Sys.glob(sprintf('%s/031.CHS/*/extra', upper_dir))) {
# 
#   path_split <- strsplit(path, '/')[[1]]
#   region_name <- path_split[length(path_split) - 1]

region_names <- 'SMN1'

for (region_name in region_names) {
  cat(sprintf('[%s]\n', region_name))
  
  all_f_values <- load_f_values(upper_dir, region_name, dataset_names)
  if (nrow(all_f_values) == 0) {
    next
  }
  
  region_groups <- unique(all_f_values$region_group)
  for (region_group in region_groups) {
    cat(sprintf('[%s: %s]\n', region_name, region_group))
    f_values <- all_f_values[all_f_values$region_group == region_group,]
    compare_f_values(region_name, f_values, plots_dir)
  }
}

region_name <- 'ADAMTS7'
region_group <- '02-01'

smn_caller_psvs <- data.frame(pos=factor(c(70950493, 70950966, 70951392, 70951463,
                70951897, 70951946, 70952094, 70952209)))

levels(factor(f_values2$psv_pos))
f_values2$psv_pos <- as.numeric(gsub('[^:]*:', '', f_values2$psv))

f_values_subs <- filter(f_values2, startsWith(stat, 'Copy') & dataset != 'Combined' &
                          psv_pos >= 70920000)
colors <- ifelse(levels(factor(f_values_subs$psv_pos)) %in% smn_caller_psvs$pos,
                 'red', 'black')
ggplot(f_values_subs) +
  geom_vline(xintercept = reliable_all, size=5, color='#ddddff') +
  geom_point(aes(factor(psv_pos), value, color=dataset),
             position=position_dodge(width=0.5)) +
  geom_hline(yintercept = RELIABLE_THRESHOLD) +
  scale_x_discrete('PSV position') +
  scale_y_continuous('f-value', breaks=c(seq(0, 1, 0.2), 0.95)) +
  scale_color_discrete('Dataset') +
  facet_wrap(~ stat, ncol=1) +
  theme_bw() +
  # ggtitle(sprintf('%s: %s', region_name, region_group)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8,
                                   color=colors))
ggsave('~/Data/hg38/jvc/plots/population_comparison/r002/smn_f_values.png',
       width=10, height=5, scale=1.2)
