wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/genes/'

extract_a <- function(..., qual=20) {
  res <- data.frame()
  for (df in list(...)) {
    curr_res <- select(df, sample, population, superpopulation)
    curr_res$Parascopy_copy_num <- ifelse(df$agCN_filter == 'PASS' &
                                            df$agCN_qual >= qual,
                              df$agCN, NA)
    n_copies <- length(grep('^psCN[0-9]$', colnames(df)))
    if (n_copies > 0) {
      for (i in 1:n_copies) {
        values <- df[[paste0('psCN', i)]]
        qualities <- df[[paste0('psCN_qual', i)]]
        curr_res[paste0('Parascopy_paralog', i)] <- ifelse(df$psCN_filter == 'PASS' &
                                                      qualities >= qual, values, NA)
      }
    }
    res <- rbind(res, curr_res)
  }
  group_by(res, sample) %>% slice_head(n = 1) %>% ungroup()
}

extract_b <- function(df, method) {
  res <- select(df, sample, population, superpopulation)
  for (column in colnames(df)) {
    if (startsWith(column, 'b_')) {
      res[sub('^b', method , column)] = df[[column]]
    }
  }
  res
}

compare_to_exper <- function(df, fun=NULL, write_all=F) {
  n_copies <- length(grep('^Parascopy_paralog[0-9]$', colnames(df)))
  methods <- sub('_copy_num', '', colnames(df)[grep('_copy_num$', colnames(df))])
  primary <- methods[1]
  for (method in methods[-1]) {
    cat(sprintf('%s & %s:\n', primary, method))
    col_a <- df[paste0(primary, '_copy_num')]
    col_b <- df[paste0(method, '_copy_num')]
    present <- !is.na(col_a) & !is.na(col_b)
    cn_match <- present & round(col_a) == round(col_b)
    perc_msg('    Aggr. CN:  ', sum(cn_match), sum(present))
    if (!is.null(fun)) {
      cat(sprintf('        %s\n', fun(col_a, col_b)))
    }
    
    if (n_copies > 0) {
      present_all <- rep(T, nrow(df))
      par_cn_match <- rep(T, nrow(df))
      for (i in 1:n_copies) {
        col_a <- df[sprintf('%s_paralog%d', primary, i)]
        col_b <- df[sprintf('%s_paralog%d', method, i)]
        present <- !is.na(col_a) & !is.na(col_b)
        curr_match <- present & round(col_a) == round(col_b)
        if (write_all) {
          perc_msg(sprintf('    Copy %d:    ', i), sum(curr_match), sum(present))
        }
        if (!is.null(fun) && write_all) {
          cat(sprintf('        %s\n', fun(col_a, col_b)))
        }
        present_all <- present_all & present
        par_cn_match <- par_cn_match & curr_match
      }
      perc_msg('    Paralog CN:', sum(par_cn_match), sum(present_all))
    }
  }
}

# ====== SMN1 ======

gene <- 'SMN1'
df <- load(gene, 'mlpa', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(2)
df <- filter(df, population != 'BEB' & population != 'TSI')
combined_16 <- filter(df, method == 'mlpa_16') %>% extract_b('MLPA')
combined_16 <- left_join(combined_16, filter(df, method == 'mlpa_16') %>% extract_a)

combined_78 <- filter(df, method == 'mlpa_78') %>% extract_b('MLPA')
combined_78 <- left_join(combined_78, filter(df, method == 'mlpa_78') %>% extract_a)

df <- load(gene, 'caller', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(2)
combined_16 <- left_join(combined_16, filter(df, method == 'caller_16') %>% extract_b('CNC'))
combined_78 <- left_join(combined_78, filter(df, method == 'caller_78') %>% extract_b('CNC'))

qm2 <- load(gene, 'qm2', keep_paralog=T, keep_b=T) %>% extend_paralog(2)
combined_16 <- left_join(combined_16, filter(qm2, region == 'SMN1_exon2') %>%
                           extract_b('QM2'))
combined_78 <- left_join(combined_78, filter(qm2, region == 'SMN1_exon7') %>%
                           extract_b('QM2'))
compare_to_exper(combined_16)
compare_to_exper(combined_78)
filter(combined_16, round(CNC_copy_num) != round(MLPA_copy_num)
       & !is.na(CNC_copy_num) & !is.na(MLPA_copy_num))

# View(filter(combined_16, round(MLPA_copy_num) != round(QM2_copy_num)))
rm(combined_16, combined_78)

# ====== RHCE ======

gene <- 'RHCE'
df <- load(gene, 'comp', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(2)
combined <- filter(df, method == 'MIP') %>% extract_b('MIP')
combined <- left_join(combined, filter(df, method == 'WGS') %>% extract_b('WGS'))

qm2 <- load(gene, 'qm2', keep_paralog=T, keep_b=T) %>% extend_paralog(2)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)
compare_to_exper(filter(combined, !is.na(QM2_copy_num)))
rm(combined)

# View(filter(combined, round(MIP_copy_num) != round(QM2_copy_num)))

# ====== AMY1C ======

gene <- 'AMY1C'
df <- load(gene, 'qpcr', keep_paralog=F, keep_b=T, keep_qual=T)
combined <- filter(df, method == 'PRT') %>% extract_b('PRT')
combined <- left_join(combined, filter(df, method == 'g1k') %>% extract_b('WGS'))

qm2 <- load(gene, 'qm2', keep_paralog=F, keep_b=T)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined, function(a, b) round(cor(a, b, use = 'pairwise'), 3))
compare_to_exper(combined, function(a, b)
  round(sum(abs(a - b), na.rm=T) / sum(!is.na(a) & !is.na(b)), 3))
rm(combined)

# ======= SRGAP2 =======

gene <- 'SRGAP2'
df <- load(gene, 'comp', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(4)
combined <- filter(df, method == 'MIP') %>% extract_b('MIP')
combined <- left_join(combined, filter(df, method == 'WGS') %>% extract_b('WGS'))

qm2 <- load(gene, 'qm2', keep_paralog=T, keep_b=T) %>% extend_paralog(4)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)
rm(combined)

combined

# View(filter(df, round(paralog4) != round(b_paralog4)))

# ====== C4A ======

gene <- 'C4A'
df <- load(gene, 'comp', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(2)
combined <- filter(df, method == 'SB') %>% extract_b('SB')
combined <- left_join(combined, filter(df, method == 'PRT') %>% extract_b('PRT'))

qm2 <- load(gene, 'qm2', keep_paralog=T, keep_b=T) %>% extend_paralog(2)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)

# ====== NPY4R ======

gene <- 'NPY4R'
df <- load(gene, 'comp', keep_paralog=F, keep_b=T, keep_qual=T)
combined <- filter(df, method == 'ddPCR') %>% extract_b('ddPCR')
# combined <- left_join(combined, filter(df, method == 'ddPCR') %>% extract_b('ddPCR'))
combined <- left_join(combined, filter(df, method == 'FREEC') %>% extract_b('FREEC'))
combined <- left_join(combined, filter(df, method == 'CNVnator') %>% extract_b('CNVnator'))

qm2 <- load(gene, 'qm2', keep_paralog=F, keep_b=T)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)

# ====== FCGR3A ======

gene <- 'FCGR3A'
df <- load(gene, 'comp', keep_paralog=T, keep_b=T, keep_qual=T) %>% extend_paralog(2)
df$b_copy_num <- as.numeric(df$b_copy_num)
main <- 'TaqMan'
combined <- filter(df, method == main) %>% extract_b(main)
for (curr in unique(df$method)) {
  if (curr != main) {
    combined <- left_join(combined, filter(df, method == curr)
                          %>% extract_b(curr))
  }
}

qm2 <- load(gene, 'qm2', keep_paralog=T, keep_b=T) %>% extend_paralog(2)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)

combined_exper_match <- filter(combined, TaqMan_copy_num == PRT_REDVR_copy_num &
         TaqMan_copy_num == SYBR_Green_copy_num)
compare_to_exper(combined_exper_match)
filter(combined_exper_match, TaqMan_paralog2 == PRT_REDVR_paralog2 &
         TaqMan_paralog2 == STR_paralog2) %>% compare_to_exper()

# ====== PMS2 ======

rm(combined)
gene <- 'PMS2'
df <- load(gene, 'comp', keep_paralog=F, keep_b=T, keep_qual=T)
combined_15 <- filter(df, method == 'exon15') %>% extract_b('Cohort')
combined_14 <- filter(df, method == 'exon14') %>% extract_b('Cohort')

qm2_m <- load(gene, 'qm2_matrix', keep_paralog=F, keep_b=T)
combined_15 <- left_join(combined_15,
                         filter(qm2_m, region == 'exon15') %>% extract_b('QM2-matrix'))
combined_14 <- left_join(combined_14,
                         filter(qm2_m, region == 'exon14') %>% extract_b('QM2-matrix'))
qm2_b <- load(gene, 'qm2_bed', keep_paralog=F, keep_b=T)
combined_15 <- left_join(combined_15,
                         filter(qm2_b, region == 'exon15') %>% extract_b('QM2-bed'))
combined_14 <- left_join(combined_14,
                         filter(qm2_b, region == 'exon14') %>% extract_b('QM2-bed'))

combined_15 <- left_join(combined_15, filter(df, method == 'exon15') %>% extract_a)
# combined_14 <- left_join(combined_14, filter(df, method == 'exon14') %>% extract_a)

compare_to_exper(combined_15)
# compare_to_exper(combined_14)
compare_to_exper(filter(combined_14, round(Cohort_copy_num) != 4))

combined_14 <- left_join(combined_14, filter(df, method == 'exon14') %>% extract_a(qual=0))
compare_to_exper(filter(combined_14, round(Cohort_copy_num) != 4))

# ====== APOBEC3A ======

gene <- 'APOBEC3A'
df <- load(gene, 'comp', keep_paralog=F, keep_b=T, keep_qual=T)
combined <- extract_b(df, 'PCR')

qm2 <- load(gene, 'qm2', keep_paralog=F, keep_b=T)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)

# ====== HYDIN ======

gene <- 'HYDIN'
df <- load(gene, 'comp', keep_paralog=F, keep_b=T, keep_qual=T)
combined <- extract_b(df, 'FISH')

qm2 <- load(gene, 'qm2', keep_paralog=F, keep_b=T)
combined <- left_join(combined, extract_b(qm2, 'QM2'))
combined <- left_join(combined, extract_a(df))
compare_to_exper(combined)
