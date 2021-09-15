#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pdf(NULL)

source('~/Code/homology-extra/psv_sample_matrix_inner.r')

if (length(args) == 0) {
  data_dir <- '~/Data/hg38/jvc/runs/201.han_g1k/./SMN1/extra'
  plots_dir <- '~/Tmp'
  important_samples <- c()
} else {
  data_dir <- trimws(args[1], which='right', whitespace='/')
  plots_dir <- trimws(args[2], which='right', whitespace='/')
  if (length(args) > 2) {
    important_samples <- args[3:length(args)]
  } else {
    important_samples <- c()
  }
}

dir_split <- unlist(strsplit(data_dir, '/', fixed=T))
stopifnot('.' %in% dir_split)
start <- match('.', dir_split)[1]
region_name <- dir_split[start + 1]
plots_prefix <- sprintf('%s/%s.', plots_dir, region_name)

cat(sprintf('[%s] %s/*.csv  ->  %s*.png\n', region_name, data_dir, plots_prefix))

load <- function(name, comment.char='#', ...) {
  filename <- sprintf('%s/%s.csv', data_dir, name)
  if (!file.exists(filename)) {
    cat(sprintf('Cannot load csv file "%s"\n', filename))
    stop(1)
  }
  read.csv(filename, sep='\t', comment.char=comment.char, ...)
}

all_sample_gts <- load('em_sample_gts')
all_likelihoods <- load('em_likelihoods')
all_f_values <- load('em_interm_f_values')
all_sample_psv_gts <- load('em_sample_psv_gts')
all_sample_psv_support <- load('em_sample_psv_support')

region_groups <- unique(all_likelihoods$region_group)
for (group in region_groups) {
  sample_gts <- filter(all_sample_gts, region_group == group)
  likelihoods <- filter(all_likelihoods, region_group == group)
  f_values <- filter(all_f_values, region_group == group)
  sample_psv_gts <- filter(all_sample_psv_gts, region_group == group)
  sample_psv_support <- filter(all_sample_psv_support, region_group == group)
  out_prefix <- sprintf('%s%s.', plots_prefix, group)
  
  tryCatch(draw_matrices(region_name, sample_gts, sample_psv_gts, sample_psv_support,
                         likelihoods, f_values, out_prefix, important_samples),
           error = function(e) {
             cat(sprintf('!!!\n[%s: %s] Could not finish: %s\n!!!\n', region_name, group, e))
           })
}
