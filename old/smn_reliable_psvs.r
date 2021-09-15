# SMN_caller hg19
# chr5:70246320-70246320
# chr5:70246793-70246793
# chr5:70247219-70247219
# chr5:70247290-70247290
# chr5:70247724-70247724
# chr5:70247773-70247773
# chr5:70247921-70247921
# chr5:70248036-70248036

# SMN_caller hg38
# chr5:70950493-70950493
# chr5:70950966-70950966
# chr5:70951392-70951392
# chr5:70951463-70951463
# chr5:70951897-70951897
# chr5:70951946-70951946
# chr5:70952094-70952094
# chr5:70952209-70952209

psv_weights <- read.csv('~/Data/hg38/jvc/runs/021/SMN1/paralog_ploidy.002/region1/psv_weights.csv',
         sep='\t', comment.char='#')
psv_weights <- psv_weights[psv_weights$iteration == max(psv_weights$iteration),]
psv_weights <- subset(psv_weights, copy1 != 'NaN')

smn_caller <- c(70246320, 70246793, 70247219, 70247290, 70247724, 70247773, 70247921, 70248036)
smn_caller <- c(70950493, 70950966, 70951392, 70951463,
                70951897, 70951946, 70952094, 70952209)
psv_weights$reliable <- psv_weights$psv %in% sprintf('chr5:%d', smn_caller)

ggplot(psv_weights) +
  geom_point(aes(copy1, copy2, color=reliable), alpha=.7) +
  scale_x_continuous('Copy 1 reliability') +
  scale_y_continuous('Copy 2 reliability') +
  theme_bw()

psv_weights$w <- pmin(psv_weights$copy1, psv_weights$copy2)
ggplot(psv_weights) +
  geom_histogram(aes(w, fill=reliable), binwidth=.01) +
  scale_x_continuous('PSV weight') +
  scale_y_continuous('Count') +
  scale_fill_discrete('SMN caller', labels=c('Unreliable', 'Reliable')) +
  theme_bw()
sum(psv_weights$w >= 0.8)  
ggsave('~/Work/presentations/november/reliable_psvs.png', width=8, height=3)
