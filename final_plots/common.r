library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')

load <- function(gene, filename, keep_paralog=T, keep_qual=F, keep_b=F) {
  df <- read.csv(sprintf('%s/%s/%s.csv', data_dir, gene, filename),
                 sep='\t', comment='#') %>%
    as_tibble
  if (!keep_b) {
    df <- select(df, !matches('b_|dist'))
  }
  if (!keep_paralog) {
    df <- select(df, !matches('paralog'))
  }
  if (!keep_qual) {
    df <- select(df, !matches('qual|filter'))
  }
  df
}

sep_column <- function(vec, ix=1, sep=',') {
  suppressWarnings(as.numeric(sapply(strsplit(vec, sep), `[`, ix)))
}

extend_paralog <- function(df, copy_num) {
  possible_colnames <- c('paralog', 'paralog_qual', 'b_paralog', 'psCN', 'psCN_qual')
  for (i in 1:copy_num) {
    for (colname in possible_colnames) {
      if (colname %in% colnames(df)) {
        df[[paste0(colname, i)]] <- sep_column(df[[colname]], i)
      }
    }
  }
  for (colname in possible_colnames) {
    df[colname] <- NULL
  }
  df
}

to_wide <- function(df) {
  pivot_wider(df, names_from = 'region',
              values_from = matches('copy_num|paralog|agCN|psCN'))
}

filter_df <- function(df, ...) {
  df_filt <- filter(df, ...)
  cat(sprintf('Removed %d samples', nrow(df) - nrow(df_filt)))
  df_filt
}

to_counts <- function(df, region_names, ..., names_prefix=NULL) {
  sample_contribution <- 100 / table(df$superpopulation)
  cn_counts <- select(df, c('sample', 'population', 'superpopulation') | all_of(...)) %>%
    pivot_longer(!c('sample', 'population', 'superpopulation'),
                 names_to='key', names_prefix) %>%
    count(superpopulation, key, value)
  cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])
  cn_counts$pop_name <- pop_names[cn_counts$superpopulation]
  cn_counts$region <- factor(region_names[cn_counts$key], levels=region_names)
  names(cn_counts$region) <- NULL
  cn_counts$value <- factor(cn_counts$value)
  cn_counts
}

ref_cns_df <- function(cn_counts, values) {
  regions <- levels(cn_counts$region)
  data.frame(region=factor(regions, levels=regions),
             cn=factor(values))
}

reorder_copies <- function(vec, ixs) {
  sapply(strsplit(vec, ','),
         function(row) do.call(paste, c(as.list(row[ixs]), sep=',')))
}

add_noise <- function(df, columns, sd=0.05) {
  n <- nrow(df)
  max_dev <- sd * 3
  for (col in columns) {
    df[[paste0(col, '_noise')]] <- df[[col]] +
      pmin(max_dev, pmax(-max_dev, rnorm(n, 0, sd)))
  }
  df
}

get_observations <- function(df, label, col1, col2, ...) {
  df_filt <- filter(df, ...)
  res <- select(df_filt, c('sample', 'population', 'superpopulation'))
  res$a_obs <- df_filt[[col1]]
  res$b_obs <- df_filt[[col2]]
  res$dist <- abs(res$a_obs - res$b_obs)
  res$label <- label
  res
}

combine_families <- function(df, pedigree) {
  df$status <- NA
  df$family <- NA
  ixs <- match(df$sample, pedigree$Individual.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'child'
  
  ixs <- match(df$sample, pedigree$Paternal.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'father'
  
  ixs <- match(df$sample, pedigree$Maternal.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'mother'
  
  df_wide <- df %>%
    select(!'sample') %>%
    pivot_wider(c('population', 'superpopulation', 'family', 'region'),
                names_from='status', values_from=matches('copy_num|paralog'))
  df_wide <- filter(df_wide, !is.na(copy_num_child) & !is.na(copy_num_father) &
                      !is.na(copy_num_mother))
  df_wide
}

perc_msg <- function(msg, x, y, fixed_size = T) {
  get_count <- function(z) {
    if (is.numeric(z)) {
      return(z)
    } else if (is.vector(z)) {
      return(length(z))
    } else {
      return(nrow(z))
    }
  }
  
  x <- get_count(x)
  y <- get_count(y)
  if (fixed_size) {
    cat(sprintf('%s %4d / %4d (%5.1f%%)\n', msg, x, y, ifelse(y > 0, x / y * 100, 0.0)))
  } else {
    cat(sprintf('%s %d / %d (%.1f%%)\n', msg, x, y, ifelse(y > 0, x / y * 100, 0.0)))
  }
}

double_filter <- function(df, ...) {
  n <- nrow(df)
  df2 <- filter(df, ...)
  m <- nrow(df2)
  perc_msg('', m, n)
}

cn_match <- function(df, ..., qual=20) {
  df <- filter(df, ...) %>%
    filter(!is.na(agCN) & !is.na(b_copy_num) & agCN_qual >= qual)
  double_filter(df, agCN == round(b_copy_num))
}

par_cn_match <- function(df, n_copies, ..., qual=20) {
  df <- filter(df, ...)
  n <- nrow(df)
  high_qual <- df$psCN_filter == 'PASS'
  par_match <- rep(T, n)

  for (i in 1:n_copies) {
    a_qual <- df[[sprintf('psCN_qual%d', i)]]
    a_par <- df[[sprintf('psCN%d', i)]]
    b_par <- df[[sprintf('b_paralog%d', i)]]
    curr_high_qual <- high_qual & a_qual >= qual & !is.na(a_par) & !is.na(b_par)
    curr_par_match <- curr_high_qual & a_par == round(b_par)
    cat(sprintf('Copy %d: ', i))
    perc_msg('', sum(curr_par_match), sum(curr_high_qual))
    
    high_qual <- high_qual & curr_high_qual
    par_match <- par_match & curr_par_match
  }
  cat('Total:  ')
  perc_msg('', sum(par_match), sum(high_qual))
}

add_axis <- function(g, x_name, y_name, title=NULL, ...) {
  if (!is.null(x_name)) {
    x_name <- textGrob(x_name, gp = gpar(...))
  }
  if (!is.null(y_name)) {
    y_name <- textGrob(y_name, rot=90, gp = gpar(...))
  }
  if (!is.null(title)) {
    title <- textGrob(title, gp = gpar(...))
  }
  grid.arrange(arrangeGrob(g, bottom = x_name, left = y_name, top = title))
}
