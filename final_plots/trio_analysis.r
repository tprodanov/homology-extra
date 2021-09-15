library(tidyverse)
library(ggplot2)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

summaries <- data.frame()
dir <- '~/Data/hg38/jvc/comparisons/children/trios/'
for (file in list.files(dir)) {
  tmp <- read_delim(file.path(dir, file), '\t',
                    col_types = cols(.default = 'c'),
                    col_names = F)
  names(tmp) <- c('chrom', 'start', 'end', 'locus',
                  'sample', 'other_regions', 'agCN_filter', 'agCN', 'agCN_qual',
                  'psCN_filter', 'psCN', 'psCN_qual', 'gene', 'exon')
  dataset <- strsplit(file, '.', fixed = T)[[1]][2]
  summaries <- rbind(summaries, add_column(tmp, dataset=dataset))
}
summaries <- mutate(summaries,
    start=as.numeric(start), end=as.numeric(end), agCN_qual=as.numeric(agCN_qual),
    agCN=as.numeric(agCN), exon=as.numeric(exon))

summaries$refCN <- with(summaries, ifelse(other_regions == '*', 2,
                       (str_count(other_regions, ',') + 2) * 2))
cns <- unique(summaries[c('gene', 'refCN')])

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
  
  df
}

pedigree <- read_delim('~/Data/hg38/data/g1k/g1k.ped', '\t')
names(pedigree) <- sub(' ', '.', names(pedigree))

QUAL <- 20

summaries$other_regions <- NULL
summaries$locus <- NULL
summaries$refCN <- NULL
summaries <- extend_paralog(summaries, 5)
summaries_fam <- combine_families(summaries, pedigree)
summaries_fam <- filter(summaries_fam,
      agCN_filter == 'PASS' & !is.na(agCN) & agCN_qual >= QUAL)
summaries_child <- filter(summaries_fam, status == 'child')
summaries_mother <- filter(summaries_fam, status == 'mother')
summaries_father <- filter(summaries_fam, status == 'father')

summaries_join1 <- left_join(summaries_child, summaries_mother,
    by=c('chrom', 'family', 'gene', 'exon'), suffix=c('.child', '.mother'))
summaries_join2 <- left_join(summaries_child, summaries_father,
    by=c('chrom', 'family', 'gene', 'exon'), suffix=c('.child', '.father'))
summaries_join <- full_join(summaries_join1, summaries_join2,
    by=names(summaries_join1)[!grepl('.mother', names(summaries_join1), fixed=T)])
summaries_join <- filter(summaries_join,
    !is.na(agCN.child) & !is.na(agCN.mother) & !is.na(agCN.father))
rm(summaries_join1, summaries_join2, summaries_fam,
   summaries_child, summaries_mother, summaries_father)

length(unique(summaries_join$sample.child))
unique(filter(summaries_join, dataset.child != 'children')$sample.child)
unique(filter(summaries_join, dataset.child != 'children')$dataset.child)

summaries_join <- left_join(summaries_join, cns, by='gene')
# write_delim(summaries_join, '~/Data/hg38/jvc/trios.csv', '\t')

# --- Analysis ---

# curr_gene <- 'SMN1'
# curr_exon <- 2
# curr_refCN <- filter(cns, gene == curr_gene)$refCN[1]
# curr_summaries <- filter(summaries_join, gene == curr_gene & exon == curr_exon)
# curr_summaries <- select(curr_summaries, !matches('psCN(|_qual)[345]'))
# curr_summaries <- filter(curr_summaries,
#   psCN_filter.child == 'PASS' & psCN_filter.father == 'PASS' & psCN_filter.mother == 'PASS'
#   & psCN_qual1.child >= QUAL & psCN_qual1.father >= QUAL & psCN_qual1.mother >= QUAL)
# 
# 
# ggplot(curr_summaries) +
#   geom_point(aes(agCN.father, agCN.mother), alpha=.1) +
#   facet_wrap(~ agCN.child)
# 
# ggplot(curr_summaries) +
#   geom_point(aes(psCN1.father, psCN1.mother), alpha=.1) +
#   facet_wrap(~ psCN1.child)
# 
# filter(curr_summaries, psCN1.child != 2) %>%
#   double_filter(psCN1.father != psCN1.child & psCN1.mother != psCN1.child)
# filter(curr_summaries, psCN1.child != 2 &
#       psCN1.father != psCN1.child & psCN1.mother != psCN1.child) %>%
#   select(matches('sample'), matches('psCN[12]')) %>% as.data.frame
# 
# # --- STRC ---
# 
# curr_gene <- 'STRC'
# curr_refCN <- filter(cns, gene == curr_gene)$refCN[1]
# curr_summaries <- filter(summaries_join, gene == curr_gene)
# curr_summaries <- select(curr_summaries, !matches('psCN(|_qual)[345]'))
# curr_summaries <- filter(curr_summaries,
#  psCN_filter.child == 'PASS' & psCN_filter.father == 'PASS' & psCN_filter.mother == 'PASS'
#  & psCN_qual1.child >= QUAL & psCN_qual1.father >= QUAL & psCN_qual1.mother >= QUAL)
# 
# ggplot(curr_summaries) +
#   geom_point(aes(agCN.father, agCN.mother), alpha=.1) +
#   facet_wrap(~ agCN.child) +
#   theme_light()
# 
# filter(curr_summaries, agCN.child == 5 & agCN.father == 4 & agCN.mother == 4)
# 
# ggplot(curr_summaries) +
#   geom_point(aes(psCN1.father, psCN1.mother), alpha=.1) +
#   facet_wrap(~ psCN1.child) +
#   theme_light()
# 
# filter(curr_summaries, psCN1.child == 1 & psCN1.father == 2 & psCN1.mother == 2) %>%
#   select(matches('sample|psCN1|psCN_qual1')) %>% as.data.frame
# 
# filter(curr_summaries, psCN1.child != 2) %>%
#   double_filter(psCN1.father != psCN1.child & psCN1.mother != psCN1.child)
# filter(curr_summaries, psCN1.child != 2 &
#          psCN1.father != psCN1.child & psCN1.mother != psCN1.child) %>%
#   select(matches('sample'), matches('psCN[12]')) %>% as.data.frame

# ========

find_discrep <- function(df, gene, paralog=0, exon=NULL, qual=20, verbose=2) {
  if (verbose == 2) {
    print_df <- function(msg, filt) {
      n <- sum(filt)
      if (n > 0) {
        cat('----------\n')
        cat(sprintf('%s:      %d\n', msg, n))
        print(select(df[filt,], matches(sprintf('sample\\.|%s\\.|%s\\.', colv, colq))) %>%
                as.data.frame %>% head)
      }
    }
  } else {
    print_df <- function(msg, filt) {}
  }

  if (is.null(exon)) {
    df <- df[df$gene == gene,]
  } else {
    df <- df[df$gene == gene & df$exon == exon,]
  }
  if (paralog == 0) {
    colf <- 'agCN_filter'
    colv <- 'agCN'
    colq <- 'agCN_qual'
    ref <- df[1,]$refCN
  } else {
    colf <- 'psCN_filter'
    colv <- sprintf('psCN%d', paralog)
    colq <- sprintf('psCN_qual%d', paralog)
    ref <- 2
  }
  
  if (verbose > 0) {
    cat(sprintf('%s: refCN = %d\n', gene, ref))
  }
  n1 <- nrow(df)
  for (indiv in c('child', 'mother', 'father')) {
    df <- df[df[[paste(colf, indiv, sep='.')]] == 'PASS' &
             df[[paste(colq, indiv, sep='.')]] >= qual &
             !is.na(df[[paste(colv, indiv, sep='.')]]),]
  }
  n2 <- nrow(df)
  if (verbose > 0) {
    perc_msg('High quality trios:', n2, n1)
  }
  if (n2 == 0) {
    return()
  }
  
  vc <- paste0(colv, '.child')
  vf <- paste0(colv, '.father')
  vm <- paste0(colv, '.mother')
  
  names <- df$sample.child
  cn_c <- df[[vc]]
  cn_f <- df[[vf]]
  cn_m <- df[[vm]]
  
  discordant <- c()
  unclear <- c()
  
  filt <- cn_m + cn_f < cn_c
  discordant <- union(discordant, names[filt])
  print_df('[disc] Indiv > Maternal + Paternal', filt)
  
  filt <- cn_c < ref & (cn_m >= ref & cn_f >= ref)
  discordant <- union(discordant, names[filt])
  print_df('[disc] Indiv < ref  &  min(Maternal, Paternal) >= ref', filt)
  
  filt <- cn_c < ref & (cn_m < cn_c & cn_f < cn_c)
  unclear <- union(unclear, names[filt])
  print_df('[uncl] Indiv < ref  &  max(Maternal, Paternal) < Indiv', filt)
  
  filt <- cn_c > ref & (cn_m < cn_c & cn_f < cn_c)
  discordant <- union(discordant, names[filt])
  print_df('[disc] Indiv > ref  &  max(Maternal, Paternal) < Indiv', filt)

  filt1 <- cn_c == ref & (cn_m <= cn_c - 2 & cn_f <= cn_c - 2)
  discordant <- union(discordant, names[filt1])
  print_df('[disc] Indiv = ref  &  max(Maternal, Paternal) <= ref - 2', filt1)

  filt2 <- cn_c == ref & (cn_m >= cn_c + 2 & cn_f >= cn_c + 2)
  discordant <- union(discordant, names[filt2])
  print_df('[disc] Indiv = ref  &  min(Maternal, Paternal) >= ref + 2', filt2)
  
  filt <- !filt1 & !filt2 & cn_c == ref & (cn_c != cn_f & cn_c != cn_m)
  unclear <- union(unclear, names[filt])
  print_df('[uncl] Indiv = ref  &  Maternal != Indiv & Paternal != Indiv', filt)

  
  concordant <- setdiff(names, union(unclear, discordant))
  if (verbose > 0) {
    cat('----------\n')
    perc_msg('Concordant:', concordant, names)
    perc_msg('Discordant:', discordant, names)
    perc_msg('Unclear:   ', unclear, names)
  }
  perc_msg('Concordant (out of clear):', concordant, c(concordant, discordant))
}

(filter(summaries_join, dataset.child == dataset.father &
         dataset.child == dataset.mother)$sample.child) %>% unique
(filter(summaries_join, dataset.child == dataset.father &
          dataset.child == dataset.mother)$dataset.child) %>% unique

find_discrep(summaries_join, 'CEL', paralog=1)

find_discrep(summaries_join, 'STRC', verbose=0)
find_discrep(summaries_join, 'STRC', paralog=1)
find_discrep(summaries_join, 'STRC', paralog=2)

find_discrep(summaries_join, 'SMN1', exon=2)
find_discrep(summaries_join, 'SMN1', exon=2, paralog=1)
find_discrep(summaries_join, 'SMN1', exon=2, paralog=2)

find_discrep(summaries_join, 'PMS2', exon=15)
find_discrep(summaries_join, 'PMS2', exon=14, qual=0)

find_discrep(summaries_join, 'NEB')
# find_discrep(summaries_join, 'SAA1', paralog=2)

find_discrep(summaries_join, 'TPSB2')

find_discrep(summaries_join, 'CFH')
find_discrep(summaries_join, 'CFH', paralog=1)
find_discrep(summaries_join, 'CFH', paralog=2)

find_discrep(summaries_join, 'ABCC6')
find_discrep(summaries_join, 'ABCC6', paralog=1)
find_discrep(summaries_join, 'ABCC6', paralog=2)

find_discrep(summaries_join, 'CEL')
find_discrep(summaries_join, 'CEL', paralog=1)
find_discrep(summaries_join, 'CEL', paralog=2)


ggplot(filter(summaries, gene == 'STRC')) +
  geom_histogram(aes(agCN_qual))
tmp <- filter(summaries, gene == 'STRC')
nrow(tmp)
sum(with(tmp, agCN_qual < 200))
quantile(tmp$agCN_qual, 0.01)


par_genes <- data.frame(gene = unique(summaries_join$gene),
                        trios = 0,
                        agcn_var = 0)
par_genes <- left_join(par_genes, cns)
for (i in 1:nrow(par_genes)) {
  curr_gene <- par_genes[i,]$gene
  df <- filter(summaries_join, gene == curr_gene)
  filt <- rep(T, nrow(df))
  for (suff in c('.child', '.father', '.mother')) {
    filt <- filt & !is.na(df[[paste0('psCN1', suff)]])
    filt <- filt & !is.na(df[[paste0('psCN2', suff)]])
    filt <- filt & df[[paste0('psCN_filter', suff)]] == 'PASS'
    filt <- filt & df[[paste0('psCN_qual1', suff)]] >= QUAL
    filt <- filt & df[[paste0('psCN_qual2', suff)]] >= QUAL
  }
  par_genes[i,]$trios <- sum(filt)
  par_genes[i,]$agcn_var <- var(df[filt,]$agCN.child)
}
filter(par_genes, agcn_var >= 0.05)
filter(par_genes, agcn_var >= 0.1)

discr_wr <- function(...) {
  for (par in 0:3) {
    cat(sprintf('Paralog %d:   ', par))
    find_discrep(summaries_join, ..., paralog=par, verbose=0)
  }
}

for (curr_gene in filter(par_genes, agcn_var >= 0.1 | refCN != 4)$gene) {
  cat(sprintf('%-10s   %d\n', curr_gene, filter(cns, gene == curr_gene)$refCN))
  discr_wr(curr_gene)
  cat('\n-------\n')
}
par_genes


centers <- read_delim('~/Data/hg38/jvc/regions/gene_centers2.bed', '\t', col_names=F)

names(centers) <- c('chrom', 'start', 'end', 'gene', 'exon')
centers
centers <- left_join(centers, cns)



library(tidyverse)
library(ggplot2)
library(ggrepel)

trios <- read_delim('~/Data/hg38/jvc/comparisons/children/trios/analysis/summary.csv',
           '\t', comment = '#') %>% filter(n_trios > 400) %>%
  mutate(perc_1 = ge_1 / n_trios * 100,
         perc_2 = ge_2 / n_trios * 100,
         perc_3 = ge_3 / n_trios * 100)
trios <- filter(trios,
                (gene != 'PMS2' | exon == 15))


scales::show_col(ggthemes::tableau_color_pal()(10))
colors <- ggthemes::tableau_color_pal()(10)[c(1, 3, 5)]

ggplot(trios) +
  geom_text_repel(aes(cn_std, perc_2, label = sprintf('%s (copy %d)', gene, paralog)),
                  data = filter(trios, cn_std > 0.75 | perc_2 < 99),
                  seed = 7, min.segment.length = 0) +
  geom_point(aes(cn_std, perc_2, color = factor(ref_cn))) +
  scale_x_continuous('psCN standard deviation',
                     breaks=seq(0, 1, 0.2)) +
  scale_y_continuous('Concordant trios (%)',
                     breaks=seq(0, 100, 2), minor_breaks=0:100) +
  scale_color_manual('refCN', values = colors) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_light() +
  theme(legend.position = c(0.032, 0.062),
        legend.justification = c('left', 'bottom'),
        legend.background = element_rect('transparent'),
        legend.key = element_rect(fill = 'transparent'),
        legend.key.height = unit(15, 'pt'))
ggsave('~/Tmp/1.png', width=10, height=6, scale=.6, dpi=450)
