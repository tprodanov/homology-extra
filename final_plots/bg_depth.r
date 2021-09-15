library(ggplot2)
library(tidyverse)
library(cowplot)

plots <- list()
plots_dir <- '~/Data/hg38/jvc/plots/bg_depth/'

borders <- c(22 - 0.5, 72 + 0.5)
filename <- '~/Data/hg38/jvc/depth/han_bgi.1603/window_depth.csv'; tech <- 'bgi'
filename <- '~/Data/hg38/jvc/depth/g1k_2500.2903/window_depth.subset.csv'; tech <- 'g1k'

sample <- 'HG00476'
name <- sprintf('%s_%s', sample, tech)

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

# ggplot(win_depth) +
#   annotate(geom='rect', xmin=-Inf, xmax=borders[1], ymin=-Inf, ymax=Inf, alpha=.1) +
#   annotate(geom='rect', xmin=borders[2], xmax=Inf, ymin=-Inf, ymax=Inf, alpha=.1) +
#   geom_boxplot(aes(gc_content, depth1, group=gc_content), outlier.alpha=0.2) +
#   geom_line(aes(gc_content, mean_loess), color='red',
#             data=bg_depth_filt) +
#   scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
#   scale_y_continuous('Read depth') +
#   theme_bw()

(plots[[sprintf('%s.a', name)]] <- ggplot(filter(bg_depth, !is.na(mean))) +
  annotate(geom='rect', xmin=-Inf, xmax=borders[1], ymin=-Inf, ymax=Inf, alpha=.1) +
  annotate(geom='rect', xmin=borders[2], xmax=Inf, ymin=-Inf, ymax=Inf, alpha=.1) +
  geom_line(aes(gc_content, mean_loess, color='Mean'), size=0.9,
            data=filter(bg_depth, gc_content >= borders[1] & gc_content <= borders[2])) +
  geom_line(aes(gc_content, var_loess, color='Variance'), size=0.9,
            data=filter(bg_depth, gc_content >= borders[1] & gc_content <= borders[2])) +
  geom_point(aes(gc_content, mean, color='Mean'), alpha=.7, size=1) +
  geom_point(aes(gc_content, var, color='Variance'), alpha=.7, size=1) +
  scale_x_continuous('GC-content', breaks=seq(0, 100, 10)) +
  scale_y_continuous('Value', limits=c(0, NA)) +
  scale_color_manual(NULL, values=ggthemes::tableau_color_pal()(2)) +
  guides(alpha = F, color = guide_legend(override.aes=list(alpha=1))) +
  theme_bw() +
  theme(legend.position=c(0.5, 1.05),
        legend.justification=c('center', 'top'),
        legend.background=element_blank(),
        legend.key = element_rect(fill='transparent')))
# ggsave(sprintf('%s/%s.a.png', plots_dir, name), width=8, height=4)

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

(plots[[sprintf('%s.b', name)]] <- ggplot(win_depth_filt) +
  geom_histogram(aes(depth1), binwidth=1, fill=ggthemes::tableau_color_pal()(1),
                 alpha=.5) +
  geom_line(aes(x, y, color=dist), data=pred, size=1) +
  annotate('text', x=Inf, y=Inf, hjust=1, vjust=1.1, size=3.5,
           label=lik_msg) +
  scale_x_continuous('Read depth') +
  scale_y_continuous('Number of 100bp windows') +
  scale_color_manual('Distribution',
                     values=RColorBrewer::brewer.pal(8, 'Dark2')) +
  theme_bw() +
  theme(legend.position=c(1, 0.4),
        legend.justification=c('right', 'center'),
        legend.background=element_blank(),
        legend.key = element_rect(fill='transparent')))
# ggsave(sprintf('%s/%s.b.png', plots_dir, name), width=8, height=4)

###############################

plot_grid(plots[[sprintf('%s_g1k.a', sample)]],
          plots[[sprintf('%s_bgi.a', sample)]],
          plots[[sprintf('%s_g1k.b', sample)]],
          plots[[sprintf('%s_bgi.b', sample)]],
          ncol=2, labels=LETTERS)
ggsave(sprintf('%s/%s.png', plots_dir, sample), width=14, height=10, scale=.6, dpi=450)
