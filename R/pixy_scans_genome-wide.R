# Hirundo genomic divergence landscape analysis
# pixy_scans_genome-wide.R - plot Fst, dxy, and π, and recombination results as genome-wide scans

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----
fst <- read.table('./results/pixy.all.1mb-100kb.fst.txt',header=T)
dxy <- read.table('./results/pixy.all.1mb-100kb.dxy.txt',header=T)
pi <- read.table('./results/pixy.all.1mb-100kb.pi.txt',header=T)
rho <- read.table('../pyrho/results/recombination.rustica.rmap.1Mb-100kb.txt',header=T,stringsAsFactors=F)

## Parse species
fst.sav.rus <- fst %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
fst.tyt.rus <- fst %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
fst.rus.aet <- fst %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
fst.rus.smi <- fst %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
fst.rus.neo <- fst %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
fst.rus.jav <- fst %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
fst.rus.dim <- fst %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
fst.aet.smi <- fst %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
fst.neo.jav <- fst %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

dxy.sav.rus <- dxy %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
dxy.tyt.rus <- dxy %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
dxy.rus.aet <- dxy %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
dxy.rus.smi <- dxy %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
dxy.rus.neo <- dxy %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
dxy.rus.jav <- dxy %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
dxy.rus.dim <- dxy %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
dxy.aet.smi <- dxy %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
dxy.neo.jav <- dxy %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

pi.sav <- pi %>% filter(str_detect(pop, 'HRS'))
pi.rus <- pi %>% filter(str_detect(pop, 'HRR'))
pi.tyt <- pi %>% filter(str_detect(pop, 'HRT'))
pi.aet <- pi %>% filter(str_detect(pop, 'HA'))
pi.smi <- pi %>% filter(str_detect(pop, 'HS'))
pi.neo <- pi %>% filter(str_detect(pop, 'HN'))
pi.jav <- pi %>% filter(str_detect(pop, 'HT'))
pi.dim <- pi %>% filter(str_detect(pop, 'HD'))

### Plot scans----

## Fst - plot 5 then output at 10 x 6, then plot second set
par(mfrow=c(5,1))
plot(fst.sav.rus$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='maroon',lwd=1.5)
plot(fst.tyt.rus$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
plot(fst.rus.aet$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='plum2',lwd=1.5)
plot(fst.rus.smi$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='goldenrod3',lwd=1.5)
plot(fst.rus.neo$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='springgreen3',lwd=1.5)
plot(fst.rus.dim$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='deepskyblue3',lwd=1.5)
plot(fst.aet.smi$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='skyblue',lwd=1.5)
plot(fst.neo.jav$avg_wc_fst,type='l',ylim=c(0,1),ylab='Fst',xlab='Chromosome Position',col='mediumseagreen',lwd=1.5)
plot(rho$rate,type='l',ylab='Recombination Rate',xlab='Chromosome Position',col='slategray',lwd=1.5)
place <- 0
for (scaff in unique(rho$chrom)){
  num <- max(as.integer(row.names(rho %>% filter(str_detect(scaff,chrom)))))
  place <- place + num
  abline(v=place)
}

## dxy - plot 5 then output at 10 x 6, then plot second set
par(mfrow=c(5,1))
plot(dxy.sav.rus$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='maroon',lwd=1.5)
plot(dxy.tyt.rus$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
plot(dxy.rus.aet$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='plum2',lwd=1.5)
plot(dxy.rus.smi$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='goldenrod3',lwd=1.5)
plot(dxy.rus.neo$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='springgreen3',lwd=1.5)
plot(dxy.rus.dim$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='deepskyblue3',lwd=1.5)
plot(dxy.aet.smi$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='skyblue',lwd=1.5)
plot(dxy.neo.jav$avg_dxy,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='mediumseagreen',lwd=1.5)
plot(rho$rate,type='l',ylab='Recombination Rate',xlab='Chromosome Position',col='slategray',lwd=1.5)
place <- 0
for (scaff in unique(rho$chrom)){
  num <- max(as.integer(row.names(rho %>% filter(str_detect(scaff,chrom)))))
  place <- place + num
  abline(v=place)
}

## pi - plot 5 then output at 10 x 6, then plot second set
par(mfrow=c(5,1))
plot(pi.rus$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
plot(pi.aet$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='plum2',lwd=1.5)
plot(pi.smi$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='goldenrod3',lwd=1.5)
plot(pi.neo$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='springgreen3',lwd=1.5)
plot(pi.jav$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='mediumseagreen',lwd=1.5)
plot(pi.dim$avg_pi,type='l',ylim=c(0,0.025),ylab='dxy',xlab='Chromosome Position',col='deepskyblue3',lwd=1.5)
plot(rho$rate,type='l',ylab='Recombination Rate',xlab='Chromosome Position',col='slategray',lwd=1.5)
place <- 0
for (scaff in unique(rho$chrom)){
  num <- max(as.integer(row.names(rho %>% filter(str_detect(scaff,chrom)))))
  place <- place + num
  abline(v=place)
}

