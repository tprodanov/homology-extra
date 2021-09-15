library(ggplot2)
library(tidyverse)

filename <- '~/Data/hg38/jvc/depth/han_bgi.1603/window_depth.csv'
plots_dir <- '~/Data/hg38/jvc/plots/bg_depth/'

sample <- 'HG00475'
tech <- 'bgi'
borders <- c(22 - 0.5, 72 + 0.5)

# ===== Load
dir <- dirname(filename)
windows <- read_delim(file.path(dir, 'windows.bed'), '\t')
names(windows)[1] <- 'chrom'
windows$window_ix <- 1:nrow(windows) - 1

win_depth_all <- read_delim(filename, '\t', comment='#')
win_depth <- win_depth_all[win_depth_all$sample == sample,]
win_depth <- left_join(win_depth, windows, by='window_ix') %>% filter(use)

bg_depth_all <- read_delim(file.path(dir, 'depth.csv'), '\t', comment='#')
bg_depth <- bg_depth_all[bg_depth_all$sample == sample & bg_depth_all$read_end == 1,]

# =====

sep_column <- function(vec, ix=1, sep=',') {
  suppressWarnings(as.numeric(sapply(strsplit(vec, sep), `[`, ix)))
}

bg_depth_filt <- filter(bg_depth, mean_loess > 0 & !is.na(mean))
#bg_depth_filt$q25 <- qnbinom(0.25, bg_depth_filt$nbinom_n, bg_depth_filt$nbinom_p)
#bg_depth_filt$q50 <- qnbinom(0.50, bg_depth_filt$nbinom_n, bg_depth_filt$nbinom_p)
#bg_depth_filt$q75 <- qnbinom(0.75, bg_depth_filt$nbinom_n, bg_depth_filt$nbinom_p)

ggplot(win_depth) +
  annotate(geom='rect', xmin=-Inf, xmax=borders[1], ymin=-Inf, ymax=Inf, alpha=.1) +
  annotate(geom='rect', xmin=borders[2], xmax=Inf, ymin=-Inf, ymax=Inf, alpha=.1) +
  geom_boxplot(aes(gc_content, depth1, group=gc_content), outlier.alpha=0.2) +
  geom_line(aes(gc_content, mean_loess), color='red',
            data=bg_depth_filt) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  theme_bw()

ggplot(filter(bg_depth, !is.na(mean))) +
  annotate(geom='rect', xmin=-Inf, xmax=borders[1], ymin=-Inf, ymax=Inf, alpha=.1) +
  annotate(geom='rect', xmin=borders[2], xmax=Inf, ymin=-Inf, ymax=Inf, alpha=.1) +
  geom_line(aes(gc_content, mean_loess, color='Mean'), size=0.9,
            data=filter(bg_depth, gc_content >= borders[1] & gc_content <= borders[2])) +
  geom_line(aes(gc_content, var_loess, color='Variance'), size=0.9,
            data=filter(bg_depth, gc_content >= borders[1] & gc_content <= borders[2])) +
  geom_point(aes(gc_content, mean, color='Mean'), alpha=.9) +
  geom_point(aes(gc_content, var, color='Variance'), alpha=.9) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Value', limits=c(0, NA)) +
  scale_color_manual(NULL, values=ggthemes::tableau_color_pal()(2)) +
  guides(alpha=F) +
  theme_bw() +
  theme(legend.position=c(0.99, 0.99), legend.justification=c('right', 'top'),
        legend.background=element_blank(),
        legend.key=element_rect(fill='gray92'))
ggsave('~/Tmp/1.png', width=8, height=4)

# =======

GC_CONT <- c(35, 45)
win_depth_filt <- filter(win_depth, gc_content >= GC_CONT[1] & gc_content <= GC_CONT[2])

n_wind <- nrow(win_depth_filt)
x <- min(win_depth_filt$depth1):max(win_depth_filt$depth1)
m <- mean(win_depth_filt$depth1)
v <- var(win_depth_filt$depth1)
nb_r = m^2 / (v - m)
nb_p = m / v

pred <- rbind(
  data.frame(x=x, y=n_wind * dnorm(x, m, sqrt(v)), dist='Normal'),
  data.frame(x=x, y=n_wind * dpois(x, m), dist='Poisson'),
  data.frame(x=x, y=n_wind * dnbinom(x, nb_r, nb_p), dist='Neg.Binom.')
)

lik_pois <- sum(dpois(win_depth_filt$depth1, m, log=T))
lik_nbinom <- sum(dnbinom(win_depth_filt$depth1, nb_r, nb_p, log=T))
lik_normal <- sum(log(pnorm(win_depth_filt$depth1 + 0.5, m, sqrt(v)) -
  pnorm(win_depth_filt$depth1 - 0.5, m, sqrt(v))))
lik_msg <- sprintf('Negative Binomial: %s \nNormal: %s \nPoisson: %s ',
                   format(lik_nbinom, big.mark=',', digits=0, scientific=F),
                   format(lik_normal, big.mark=',', digits=0, scientific=F),
                   format(lik_pois, big.mark=',', digits=0, scientific=F))

ggplot(win_depth_filt) +
  geom_histogram(aes(depth1), binwidth=1, fill=ggthemes::tableau_color_pal()(1),
                 alpha=.5) +
  geom_line(aes(x, y, color=dist), data=pred, size=1) +
  annotate('text', x=Inf, y=Inf, hjust=1, vjust=1.1,
           label=lik_msg) +
  scale_x_continuous('Read depth') +
  scale_y_continuous('Number of 100bp windows') +
  scale_color_manual('Distribution',
                     values=RColorBrewer::brewer.pal(8, 'Dark2')) +
                     #values=ggthemes::tableau_color_pal()(10)[c(2, 3, 4)]) +
  theme_bw() +
  theme(legend.position=c(1, 0.5), legend.justification=c('right', 'center'),
        legend.background=element_blank())
ggsave('~/Tmp/1.png', width=8, height=4)

###############################




sample <- 'NA18582'
dir <- 'bgi_old'; tech <- 'BGI (old)'
dir <- 'bgi'; tech <- 'BGI'
dir <- 'g1k_old'; tech <- 'Illumina (old)'
dir <- 'g1k'; tech <- 'Illumina'
dir <- 'bgi.004'; tech <- 'BGI (004)'

title <- sprintf('%s: %s', sample, tech)
win_depth <- read_delim(sprintf('%s/window_depth.csv', dir), '\t', comment='#')
win_depth <- mutate(win_depth, total = depth1 + depth2) %>%
  mutate(many_low_mapq = low_mapq >= 0.1 * total,
         many_clipped = clipped >= 0.1 * total,
         many_unpaired = unpaired >= 0.1 * total) %>%
  mutate(good = !many_low_mapq & !many_clipped & !many_unpaired)
win_depth <- left_join(win_depth, windows, 'window_ix')

sum(win_depth$many_clipped) / nrow(win_depth) * 100
sum(win_depth$many_low_mapq) / nrow(win_depth) * 100
sum(win_depth$many_unpaired) / nrow(win_depth) * 100
mean(win_depth$depth1)
sum(win_depth$depth1) / 1e6

###############################

depth_subset <- filter(win_depth, good & 35 <= gc_count & gc_count <= 50)
m <- mean(depth_subset$depth1)
v <- var(depth_subset$depth1)
distr <- data.frame(depth1 = 0:max(depth_subset$depth1))
n <- nrow(depth_subset)

distr$poisson <- dpois(distr$depth1, m) * n
nb_r = m^2 / (v - m)
nb_p = m / v
distr$nbinom <- dnbinom(distr$depth1, nb_r, nb_p) * n
distr$gauss <- dnorm(distr$depth1, m, sqrt(v)) * n

distr$obs <- as.vector(table(depth_subset$depth1)[as.character(distr$depth1)])
distr$obs[is.na(distr$obs)] <- 0

err_nbinom <- sum(with(distr, (nbinom / n - obs / n)^2))
err_gauss <- sum(with(distr, (gauss / n - obs / n)^2))
err_poisson <- sum(with(distr, (poisson / n - obs / n)^2))

ggplot(distr, aes(depth1)) +
  geom_point(aes(y = obs), size=2) +
  geom_line(aes(y = nbinom, color = 'Neg. Binomial')) +
  geom_line(aes(y = poisson, color = 'Poisson')) +
  geom_line(aes(y = gauss, color = 'Gaussian')) +
  ggtitle(sprintf('%s. Background depth (GC-content 40-50)', title)) +
  scale_x_continuous('Read depth') +
  scale_y_continuous('Number of windows') +
  scale_color_discrete('Distribution') +
  annotate('text', x = Inf, y = Inf, vjust=1.1, hjust=1,
           label=sprintf('Errors:  \nNeg.Binomial: %.6f  \nGaussian: %.6f  \nPoisson: %.6f  ',
                         err_nbinom, err_gauss, err_poisson)) +
  theme_bw()
ggsave(sprintf('%s/distr_%s.png', plots_dir, dir), width=9, height=5)


##############################

# depth1 < 90 because it is an outlier that make infinite gaussian probabilities.
#good_depth <- filter(win_depth, good & depth1 < 90)
good_depth <- filter(win_depth, good)

prediction <- data.frame(gc_count = 0:100)
loess_fit <- loess(depth1 ~ gc_count, good_depth, se=F, span=.1)
prediction$loess_def <- predict(loess_fit, prediction$gc_count)
#prediction$loess_02 <- pmax(0.01, predict(loess_fit, prediction$gc_count))
#prediction$loess_05 <- pmax(0.01, predict(loess_fit, prediction$gc_count))
prediction$loess_01 <- pmax(0.01, predict(loess_fit, prediction$gc_count))
prediction$loess <- prediction$loess_01

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = loess, group = 1), data = prediction, color = 'red') +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  theme_bw()

nrow(good_depth)
nrow(filter(good_depth, gc_count >= 80 | gc_count <= 20))

##############################

create_quantiles <- function(FUN, params, prefix='') {
  quant <- apply(params, 1, function(row) {
    l <- c(list(c(0.25, 0.5, 0.75)), as.list(row))
    do.call(FUN, l)
  })
  res <- data.frame(gc_content = rep(0:100, 3))
  res$label <- rep(sprintf('%s%s', prefix, c(25, 50, 75)), each=101)
  res$value = as.vector(t(quant))
  res
}

calculate_likelihood <- function(read_depth, FUN, params) {
  counts <- count(read_depth, gc_count, depth1)
  lik <- 0
  for (i in 1:nrow(counts)) {
    gc_count <- counts$gc_count[i]
    depth <- counts$depth1[i]
    n <- counts$n[i]
    
    l <- c(list(depth), as.list(params[gc_count + 1,]), list(log = T))
    curr_lik <- do.call(FUN, l)
    if (is.infinite(curr_lik)) {
      cat(sprintf('GC-count %2d, depth %2d -> %.4f (*%d = %.4f)\n',
                  gc_count, depth, curr_lik, n, n * curr_lik))
    }
    lik <- lik + n * curr_lik
  }
  lik
}

###############################

# Simple prediction (same everywhere)
nb_r = m^2 / (v - m)
nb_p = m / v
params1 <- matrix(rep(c(nb_r, nb_p), 101), ncol=2, byrow=T)
quant1 <- create_quantiles(qnbinom, params1)
lik1 <- calculate_likelihood(good_depth, dnbinom, params1)

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = value, color = label), data = quant1, size = 1.) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  scale_color_discrete('Quantile') +
  annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
           label=sprintf('Likelihood: %.0f  ', lik1)) +
  ggtitle(sprintf('%s. Simple prediction', title)) +
  theme_bw()
ggsave(sprintf('%s/pred1_%s.png', plots_dir, dir), width=9, height=5)

##############################

# params2 <- matrix(c(prediction$loess^2 / (v - prediction$loess),
#                     rep(nb_p, 101)), ncol=2)
# quant2 <- create_quantiles(qnbinom, params2)
# lik2 <- calculate_likelihood(good_depth, dnbinom, params2)
# 
# ggplot(good_depth, aes(gc_count)) +
#   geom_boxplot(aes(y = depth1, group = gc_count)) +
#   geom_line(aes(y = value, color = label), data = quant2, size=1.) +
#   scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
#   scale_y_continuous('Read depth') +
#   scale_color_discrete('Quantile') +
#   annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
#            label=sprintf('Likelihood: %.0f  ', lik2)) +
#   ggtitle(sprintf('%s. Constant p parameter', title)) +
#   theme_bw()
# ggsave(sprintf('%s/pred2_%s.png', plots_dir, dir), width=9, height=5)

params3 <- matrix(c(prediction$loess^2 / (pmax(v, prediction$loess + 0.01) - prediction$loess),
                    prediction$loess / pmax(v, prediction$loess + 0.01)), ncol=2)
quant3 <- create_quantiles(qnbinom, params3)
lik3 <- calculate_likelihood(good_depth, dnbinom, params3)

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = value, color = label), data = quant3, size=1.) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  scale_color_discrete('Quantile') +
  annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
           label=sprintf('Likelihood: %.0f  ', lik3)) +
  ggtitle(sprintf('%s. Constant variance (loess span 0.5)', title)) +
  theme_bw()
ggsave(sprintf('%s/pred3c_%s.png', plots_dir, dir), width=9, height=5)

##############################

# Normal distribution

good_depth$norm_depth <- good_depth$depth1 / prediction$loess[good_depth$gc_count + 1]

ggplot(filter(good_depth, gc_count >= 25 & gc_count <= 75)) +
  geom_histogram(aes(norm_depth), binwidth=.1)
norm_v <- var(filter(good_depth, gc_count >= 25 & gc_count <= 75)$norm_depth)

params10 <- matrix(c(prediction$loess,
                     sqrt(prediction$loess^2 * norm_v)), ncol=2)
quant10 <- create_quantiles(qnorm, params10)

norm_integer_prob <- function(x, mean, sd, log = F) {
  a <- pnorm(x - 0.5, mean, sd, log.p = log)
  b <- pnorm(x + 0.5, mean, sd, log.p = log)
  if (log) {
    ifelse(is.infinite(a), b, base::log(exp(b) - exp(a)))
  } else {
    b - a
  }
}

lik10 <- calculate_likelihood(good_depth, norm_integer_prob, params10)

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = value, color = label), data = quant10, size=1.) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  scale_color_discrete('Quantile') +
  annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
           label=sprintf('Likelihood: %.0f  ', lik10)) +
  ggtitle(sprintf('%s. Gaussian', title)) +
  theme_bw()
ggsave(sprintf('%s/pred10_%s.png', plots_dir, dir), width=9, height=5)

##############################

vars <- prediction$loess^2 * norm_v
vars_nb <- pmax(prediction$loess + 0.01, vars)
params4 <- matrix(c(prediction$loess^2 / (vars_nb - prediction$loess),
                    prediction$loess / vars_nb), ncol=2)
quant4 <- create_quantiles(qnbinom, params4)
lik4 <- calculate_likelihood(good_depth, dnbinom, params4)

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = value, color = label), data = quant4, size=1.) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  scale_color_discrete('Quantile') +
  annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
           label=sprintf('Likelihood: %.0f  ', lik4)) +
  ggtitle(sprintf('%s. Neg.Binom. scaled variance', title)) +
  theme_bw()
ggsave(sprintf('%s/pred4_%s.png', plots_dir, dir), width=9, height=5)

##############################

max_var <- max(prediction$loess^2 * norm_v, na.rm = T)
params5 <- matrix(c(prediction$loess^2 / (max_var - prediction$loess),
                    prediction$loess / max_var), ncol=2)
quant5 <- create_quantiles(qnbinom, params5)
lik5 <- calculate_likelihood(good_depth, dnbinom, params5)

ggplot(good_depth, aes(gc_count)) +
  geom_boxplot(aes(y = depth1, group = gc_count)) +
  geom_line(aes(y = value, color = label), data = quant5, size=1.) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Read depth') +
  scale_color_discrete('Quantile') +
  annotate('text', x = Inf, y = Inf, vjust=1.3, hjust=1,
           label=sprintf('Likelihood: %.0f  ', lik5)) +
  ggtitle(sprintf('%s. Neg.Binom. const variance (bigger)', title)) +
  theme_bw()
ggsave(sprintf('%s/pred5_%s.png', plots_dir, dir), width=9, height=5)

##############################

values <- filter(win_depth, good)$depth1
#values <- filter(win_depth, good & gc_count >= 30 & gc_count <= 60)$depth1
#values <- filter(win_depth, good & gc_count >= 50 & gc_count <= 60)$depth1
#values <- win_depth$depth1
m <- mean(values)
v <- var(values)

good_win_depth <- filter(win_depth, good)
df <- aggregate(depth1 ~ gc_count, good_win_depth, mean) %>%
  rename(mean = depth1)
df$var <- aggregate(depth1 ~ gc_count, good_win_depth, var)$depth1

good_win_depth$norm_depth <- good_win_depth$depth1 /
  df$mean[match(good_win_depth$gc_count, df$gc_count)]
df$norm_mean <- aggregate(norm_depth ~ gc_count, good_win_depth, mean)$norm_depth
df$norm_var <- aggregate(norm_depth ~ gc_count, good_win_depth, var)$norm_depth

ggplot(df) +
  geom_line(aes(gc_count, norm_mean)) +
  geom_line(aes(gc_count, norm_var), color='red')

print(c(m, v))

nb_r = m^2 / (v - m)
nb_p = m / v
nb_df <- data.frame(x = 0:max(round(values)))
nb_df$d <- dnbinom(nb_df$x, nb_r, nb_p) * length(values)
nb_df$poisson <- dpois(nb_df$x, m) * length(values)

ggplot(win_depth) +
  geom_histogram(aes(depth1, fill=good), binwidth=1) +
  geom_line(aes(x, d, color='Neg. Binom'), data=nb_df) +
  geom_line(aes(x, poisson, color='Poisson'), data=nb_df) +
  scale_x_continuous('Read depth', limits=c(0, 80)) +
  scale_y_continuous('Number of windows') +
  scale_color_manual('Distribution', values=c('black', 'blue')) +
  scale_fill_discrete('Good windows') +
  ggtitle('NA18582 (BGI) - background read depth') +
  theme_bw()
ggsave(sprintf('%s/bgi_bg_depth.png', plots_dir), width=9, height=5)



depth_subset <- filter(win_depth, good & 40 <= gc_count & gc_count <= 50)
m <- mean(depth_subset$depth1)
v <- var(depth_subset$depth1)
distr <- data.frame(depth1 = 0:max(depth_subset$depth1))

distr$poisson <- dpois(distr$depth1, m) * nrow(depth_subset)
nb_r = m^2 / (v - m)
nb_p = m / v
distr$nbinom <- dnbinom(distr$depth1, nb_r, nb_p) * nrow(depth_subset)
distr$gauss <- dnorm(distr$depth1, m, sqrt(v)) * nrow(depth_subset)

distr$obs <- as.vector(table(depth_subset$depth1)[as.character(distr$depth1)])
distr$obs[is.na(distr$obs)] <- 0

ggplot(distr, aes(depth1)) +
  geom_point(aes(y = obs), size=2) +
  geom_line(aes(y = nbinom, color = 'Neg. Binomial')) +
  geom_line(aes(y = poisson, color = 'Poisson')) +
  geom_line(aes(y = gauss, color = 'Gaussian')) +
  ggtitle('NA18582 background depth (BGI)') +
  scale_x_continuous('Read depth') +
  scale_y_continuous('Number of windows') +
  scale_color_discrete('Distribution') +
  theme_bw()
ggsave(sprintf('%s/distr_bgi.png', plots_dir), width=9, height=5)


sum(with(distr, abs(nbinom - obs)))
sum(with(distr, abs(gauss - obs)))
sum(with(distr, abs(poisson - obs)))

#####






loess_depth <- read_delim('NA18582/depth_bgi/depth.csv', '\t')

win_depth$gc_factor <- factor(sprintf('%03d', win_depth$gc_count))
loess_depth$gc_factor <- factor(sprintf('%03d', loess_depth$gc_count))

ggplot(loess_depth, aes(gc_factor, depth1)) +
  geom_boxplot(data=filter(win_depth, good)) +
  geom_line(aes(group=1), color='red') +
  scale_x_discrete('GC-content',
    labels=function(x) {
      x <- as.character(x)
      ifelse(endsWith(x, '0'), gsub('^0{1,2}', '', x), '')
    }) +
  scale_y_continuous('Read depth') +
  ggtitle('NA18582 (BGI) - background read depth') +
  theme_bw()
ggsave(sprintf('%s/bgi_bg_gc_content.png', plots_dir), width=9, height=5)


depth_dupl <- read_delim('runs/001.han_bgi/FRMPD2/extra/depth.csv', '\t', comment='#')
depth_dupl <- filter(depth_dupl, sample == 'NA18582')
windows_dupl <- read_delim('runs/001.han_bgi/FRMPD2/extra/windows.bed', '\t')
names(windows_dupl)[1] <- 'chrom'

depth_dupl <- left_join(depth_dupl, windows_dupl, 'window_ix')
ggplot(depth_dupl) +
  geom_line(aes(window_ix, depth1, color=factor(ploidy))) +
  theme_bw()

ggplot(depth_dupl) +
  geom_line(aes(window_ix, depth1 / bg_depth1 * 2, color=factor(ploidy))) +
  theme_bw()

ggplot(depth_dupl) +
  geom_line(aes(window_ix, gc_count, color=factor(ploidy))) +
  theme_bw()

ggplot(depth_dupl) +
  geom_point(aes(gc_count, depth1, color=factor(ploidy)), alpha=.5) +
  scale_color_discrete('Copy number') +
  scale_x_continuous('GC-content') +
  scale_y_continuous('Read depth') +
  ggtitle('NA18582 (BGI): FRMPD2 duplication') +
  theme_bw()
ggsave(sprintf('%s/bgi_dupl.png', plots_dir), width=9, height=5)

ggplot(depth_dupl) +
  geom_line(aes(window_ix, depth1 / bg_depth1 * 2,
                color=factor(round(gc_count / 20)), group=1)) +
  theme_bw()

ggplot(depth_dupl) +
  geom_line(aes(window_ix, depth1 / bg_depth1 * 2)) +
  geom_line(aes(window_ix, gc_count / 50, color='GC-count')) +
  theme_bw()

depth_dupl <- mutate(depth_dupl, total = depth1 + depth2) %>%
  mutate(many_low_mapq = low_mapq >= 0.1 * total,
         many_clipped = clipped >= 0.1 * total,
         many_unpaired = unpaired >= 0.1 * total) %>%
  mutate(good = !many_low_mapq & !many_clipped & !many_unpaired)

ggplot(depth_dupl) +
  geom_histogram(aes(depth1, fill=!many_low_mapq & !many_unpaired), binwidth = 1) +
  facet_wrap(~ ploidy) +
  theme_bw()

depth_dupl2 <- filter(depth_dupl, ploidy == 2)
nb_df2 <- data.frame(x = 0:max(round(depth_dupl2$depth1)))
nb_df2$d <- dnbinom(nb_df2$x, nb_r, nb_p) * nrow(depth_dupl2)
ggplot(depth_dupl2) +
  geom_histogram(aes(depth1, fill=good), binwidth=1) +
  geom_line(aes(x, d), data=nb_df2) +
  theme_bw()

depth_dupl4 <- filter(depth_dupl, ploidy == 4)
nb_df4 <- data.frame(x = 0:max(round(depth_dupl4$depth1)))
nb_df4$d <- dnbinom(nb_df4$x, nb_r * 2, nb_p) * nrow(depth_dupl4)
ggplot(depth_dupl4) +
  geom_histogram(aes(depth1, fill=good), binwidth=1) +
  geom_line(aes(x, d), data=nb_df4) +
  theme_bw()


setwd('~/Data/hg38/jvc/depth/NA18582/bgi.005')
setwd('~/Data/hg38/jvc/depth/NA18582/g1k.005')
bg_depth <- read_delim('depth.csv', '\t', comment = '#')
windows <- read_delim('windows.bed', '\t')
names(windows)[1] <- 'chrom'
windows$window_ix <- 1:nrow(windows) - 1
depth <- read_delim('window_depth.csv', '\t', comment = '#')
depth <- left_join(depth, windows, 'window_ix')
rm(windows)

for (i in 1:5) {
  bg_depth[sprintf('Q%d', 25 * (i - 1))] <- as.numeric(
    sapply(strsplit(bg_depth$quartiles, ',', fixed=T), `[`, i))
}
bg_depth$read_end <- factor(bg_depth$read_end)

ggplot(bg_depth, aes(gc_content, color = read_end)) +
  geom_point(aes(y = mean)) +
  geom_line(aes(y = mean_loess)) +
  theme_bw()
ggplot(bg_depth, aes(gc_content, color = read_end)) +
  geom_point(aes(y = var)) +
  geom_line(aes(y = var_loess)) +
  theme_bw()

ggplot(filter(bg_depth, read_end == '1'), aes(gc_content)) +
  #geom_line(aes(y = Q0, color='Q0')) +
  geom_line(aes(y = Q25, color='Q25')) +
  geom_line(aes(y = Q50, color='Q50')) +
  geom_line(aes(y = Q75, color='Q75')) +
  #geom_line(aes(y = Q100, color='Q100')) +
  geom_line(aes(y = mean_loess, color='loess')) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  theme_bw()

params <- as.matrix(filter(bg_depth, read_end == '1') %>% select('nbinom_n', 'nbinom_p'))
colnames(params) <- NULL
quant <- create_quantiles(qnbinom, params)

ggplot(filter(bg_depth, read_end == '1'), aes(gc_content)) +
  geom_line(aes(y = Q25, color='25'), size=2) +
  geom_line(aes(y = Q50, color='50'), size=2) +
  geom_line(aes(y = Q75, color='75'), size=2) +
  geom_line(aes(y = value, color = label), data=quant, linetype='dashed') +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  theme_bw()

mean_var <- data.frame(gc_count = 0:100, count=NA, mean=NA, var=NA)
for (i in 1:nrow(mean_var)) {
  subs <- filter(depth, gc_count == mean_var$gc_count[i])
  mean_var$count[i] <- nrow(subs)
  mean_var$mean[i] <- mean(subs$depth1)
  mean_var$var[i] <- var(subs$depth1)
}

mean_var$perc1 <- mean_var$count >= 0.01 * sum(mean_var$count)
ggplot(mean_var) +
  geom_line(aes(gc_count, var, color = perc1, group=1)) +
  geom_line(aes(gc_count, mean, color = perc1, group = 1))

var_loess <- loess(var ~ gc_count, filter(mean_var, perc1))
mean_var$var_pred <- predict(var_loess, mean_var$gc_count)

ggplot(filter(mean_var, perc1)) +
  geom_line(aes(gc_count, var)) +
  geom_line(aes(gc_count, var_pred), color='blue')

ggplot(mean_var) +
  geom_line(aes(gc_count, count / sum(count)))
