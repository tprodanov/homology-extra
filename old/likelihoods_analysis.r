library(dplyr)
library(viridis)
library(ggplot2)

dir <- '~/Data/hg38/jvc/runs/021/'
Sys.glob('*/paralog_ploidy/*/em_likelihoods.csv')

likelihoods <- read.csv('FCGR2B/paralog_ploidy/02-01/em_likelihoods.csv',
         sep='\t', comment.char='#')
likelihoods <- likelihoods %>% group_by(cluster) %>% top_n(1, iteration)
which.max(likelihoods$likelihood)
likelihoods[best_lik,]$likelihood - likelihoods[1,]$likelihood

ggplot(likelihoods) +
  geom_point(aes(likelihood, reliable_info, color=n_reliable)) +
  scale_color_viridis()

for (filename in Sys.glob(file.path(dir, '*/paralog_ploidy/*/em_likelihoods.csv'))) {
  print(filename)
  likelihoods <- read.csv(filename, sep='\t', comment.char='#')
  likelihoods <- likelihoods %>% group_by(cluster) %>% top_n(1, iteration)
  print(likelihoods)
  best_lik <- which.max(likelihoods$likelihood)
  max_info <- max(likelihoods$reliable_info)
  err <- ''
  if (likelihoods[best_lik,]$likelihood - likelihoods[1,]$likelihood > 20) {
    err <- sprintf('First cluster is not optimal (%.5f - %.5f = %.5f)',
                   likelihoods[best_lik,]$likelihood, likelihoods[1,]$likelihood,
                   likelihoods[best_lik,]$likelihood - likelihoods[1,]$likelihood)
  }
  
  if (is.nan(max_info)) {
    max_info <- 100
    err <- 'There is a cluster with no reliable PSVs'
  }
  
  if (max_info < 0.8) {
    err <- sprintf('Low maximal information content: %.3f', max_info)
  }

  max_n_reliable <- max(subset(likelihoods, reliable_info > max_info - 0.1)$n_reliable)
  if (likelihoods[best_lik, 'n_reliable'] < max_n_reliable / 2) {
    err <- 'Few reliable PSVs'
  }
  
  if (err != '') {
    cat(sprintf('%s in %s\n', err, filename))
    print(likelihoods)
    cat('\n')
  }
}
