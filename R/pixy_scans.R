# Hirundo genomic divergence landscape analysis
# pixy_scans.R - plot Fst, dxy, and π results from Pixy

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----
fst.1a <- read.table('./results/pixy.1mb-100kb.chr1A_fst.txt',header=T)
dxy.1a <- read.table('./results/pixy.1mb-100kb.chr1A_dxy.txt',header=T)
pi.1a <- read.table('./results/pixy.1mb-100kb.chr1A_pi.txt',header=T)

fst.4 <- read.table('./results/pixy.1mb-100kb.chr4_fst.txt',header=T)
dxy.4 <- read.table('./results/pixy.1mb-100kb.chr4_dxy.txt',header=T)
pi.4 <- read.table('./results/pixy.1mb-100kb.chr4_pi.txt',header=T)

fst.z <- read.table('./results/pixy.1mb-100kb.chrZ_fst.txt',header=T)
dxy.z <- read.table('./results/pixy.1mb-100kb.chrZ_dxy.txt',header=T)
pi.z <- read.table('./results/pixy.1mb-100kb.chrZ_pi.txt',header=T)

## Parse species
fst.1a.sav.rus <- fst.1a %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
fst.1a.tyt.rus <- fst.1a %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
fst.1a.rus.aet <- fst.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
fst.1a.rus.smi <- fst.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
fst.1a.rus.neo <- fst.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
fst.1a.rus.jav <- fst.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
fst.1a.rus.dim <- fst.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
fst.1a.aet.smi <- fst.1a %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
fst.1a.neo.jav <- fst.1a %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

dxy.1a.sav.rus <- dxy.1a %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
dxy.1a.tyt.rus <- dxy.1a %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
dxy.1a.rus.aet <- dxy.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
dxy.1a.rus.smi <- dxy.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
dxy.1a.rus.neo <- dxy.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
dxy.1a.rus.jav <- dxy.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
dxy.1a.rus.dim <- dxy.1a %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
dxy.1a.aet.smi <- dxy.1a %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
dxy.1a.neo.jav <- dxy.1a %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

pi.1a.sav <- pi.1a %>% filter(str_detect(pop, 'HRS'))
pi.1a.rus <- pi.1a %>% filter(str_detect(pop, 'HRR'))
pi.1a.tyt <- pi.1a %>% filter(str_detect(pop, 'HRT'))
pi.1a.aet <- pi.1a %>% filter(str_detect(pop, 'HA'))
pi.1a.smi <- pi.1a %>% filter(str_detect(pop, 'HS'))
pi.1a.neo <- pi.1a %>% filter(str_detect(pop, 'HN'))
pi.1a.jav <- pi.1a %>% filter(str_detect(pop, 'HT'))
pi.1a.dim <- pi.1a %>% filter(str_detect(pop, 'HD'))

fst.4.sav.rus <- fst.4 %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
fst.4.tyt.rus <- fst.4 %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
fst.4.rus.aet <- fst.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
fst.4.rus.smi <- fst.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
fst.4.rus.neo <- fst.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
fst.4.rus.jav <- fst.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
fst.4.rus.dim <- fst.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
fst.4.aet.smi <- fst.4 %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
fst.4.neo.jav <- fst.4 %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

dxy.4.sav.rus <- dxy.4 %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
dxy.4.tyt.rus <- dxy.4 %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
dxy.4.rus.aet <- dxy.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
dxy.4.rus.smi <- dxy.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
dxy.4.rus.neo <- dxy.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
dxy.4.rus.jav <- dxy.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
dxy.4.rus.dim <- dxy.4 %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
dxy.4.aet.smi <- dxy.4 %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
dxy.4.neo.jav <- dxy.4 %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

pi.4.sav <- pi.4 %>% filter(str_detect(pop, 'HRS'))
pi.4.rus <- pi.4 %>% filter(str_detect(pop, 'HRR'))
pi.4.tyt <- pi.4 %>% filter(str_detect(pop, 'HRT'))
pi.4.aet <- pi.4 %>% filter(str_detect(pop, 'HA'))
pi.4.smi <- pi.4 %>% filter(str_detect(pop, 'HS'))
pi.4.neo <- pi.4 %>% filter(str_detect(pop, 'HN'))
pi.4.jav <- pi.4 %>% filter(str_detect(pop, 'HT'))
pi.4.dim <- pi.4 %>% filter(str_detect(pop, 'HD'))

fst.z.sav.rus <- fst.z %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
fst.z.tyt.rus <- fst.z %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
fst.z.rus.aet <- fst.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
fst.z.rus.smi <- fst.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
fst.z.rus.neo <- fst.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
fst.z.rus.jav <- fst.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
fst.z.rus.dim <- fst.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
fst.z.aet.smi <- fst.z %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
fst.z.neo.jav <- fst.z %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

dxy.z.sav.rus <- dxy.z %>% filter(str_detect(pop1, 'HRS') & str_detect(pop2,'HRR'))
dxy.z.tyt.rus <- dxy.z %>% filter(str_detect(pop1, 'HRT') & str_detect(pop2,'HRR'))
dxy.z.rus.aet <- dxy.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HA'))
dxy.z.rus.smi <- dxy.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HS'))
dxy.z.rus.neo <- dxy.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HN'))
dxy.z.rus.jav <- dxy.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HT'))
dxy.z.rus.dim <- dxy.z %>% filter(str_detect(pop1, 'HRR') & str_detect(pop2,'HD'))
dxy.z.aet.smi <- dxy.z %>% filter(str_detect(pop1, 'HA') & str_detect(pop2,'HS'))
dxy.z.neo.jav <- dxy.z %>% filter(str_detect(pop1, 'HN') & str_detect(pop2,'HT'))

pi.z.sav <- pi.z %>% filter(str_detect(pop, 'HRS'))
pi.z.rus <- pi.z %>% filter(str_detect(pop, 'HRR'))
pi.z.tyt <- pi.z %>% filter(str_detect(pop, 'HRT'))
pi.z.aet <- pi.z %>% filter(str_detect(pop, 'HA'))
pi.z.smi <- pi.z %>% filter(str_detect(pop, 'HS'))
pi.z.neo <- pi.z %>% filter(str_detect(pop, 'HN'))
pi.z.jav <- pi.z %>% filter(str_detect(pop, 'HT'))
pi.z.dim <- pi.z %>% filter(str_detect(pop, 'HD'))

## Read in and parse recombination rate
rho <- read.table('/Users/drewschield/Library/CloudStorage/GoogleDrive-drew.schield@gmail.com/My Drive/projects/hirundo_divergence_landscape/analysis/pyrho/results/recombination.rustica.rmap.1Mb-100kb.txt',header=T,stringsAsFactors=F)
rho.1a <- rho %>% filter(str_detect(chrom, 'NC_053453.1'))
rho.4 <- rho %>% filter(str_detect(chrom, 'NC_053454.1'))
rho.z <- rho %>% filter(str_detect(chrom, 'NC_053488.1'))

### Plot scans on focal chromosomes----

## Scale x-axis to the longest of the three chromosomes

## Fst
par(mfrow=c(4,3))
plot(fst.1a.sav.rus$window_pos_1,fst.1a.sav.rus$avg_wc_fst,type='l',ylim=c(0,1),xlim=c(0,9.02e+07),ylab='Fst',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(fst.1a.tyt.rus$window_pos_1,fst.1a.tyt.rus$avg_wc_fst,col='darkorange3',lwd=1.5)
lines(fst.1a.rus.aet$window_pos_1,fst.1a.rus.aet$avg_wc_fst,col='plum2',lwd=1.5)
lines(fst.1a.rus.smi$window_pos_1,fst.1a.rus.smi$avg_wc_fst,col='goldenrod3',lwd=1.5)
lines(fst.1a.rus.neo$window_pos_1,fst.1a.rus.neo$avg_wc_fst,col='springgreen3',lwd=1.5)
#lines(fst.1a.rus.jav$window_pos_1,fst.1a.rus.jav$avg_wc_fst,col='darkgreen',lwd=1.5)
lines(fst.1a.rus.dim$window_pos_1,fst.1a.rus.dim$avg_wc_fst,col='deepskyblue3',lwd=1.5)
lines(fst.1a.aet.smi$window_pos_1,fst.1a.aet.smi$avg_wc_fst,col='skyblue',lwd=1.5)
lines(fst.1a.neo.jav$window_pos_1,fst.1a.neo.jav$avg_wc_fst,col='mediumseagreen',lwd=1.5)

plot(fst.4.sav.rus$window_pos_1,fst.4.sav.rus$avg_wc_fst,type='l',ylim=c(0,1),xlim=c(0,9.02e+07),ylab='Fst',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(fst.4.tyt.rus$window_pos_1,fst.4.tyt.rus$avg_wc_fst,col='darkorange3',lwd=1.5)
lines(fst.4.rus.aet$window_pos_1,fst.4.rus.aet$avg_wc_fst,col='plum2',lwd=1.5)
lines(fst.4.rus.smi$window_pos_1,fst.4.rus.smi$avg_wc_fst,col='goldenrod3',lwd=1.5)
lines(fst.4.rus.neo$window_pos_1,fst.4.rus.neo$avg_wc_fst,col='springgreen3',lwd=1.5)
#lines(fst.4.rus.jav$window_pos_1,fst.4.rus.jav$avg_wc_fst,col='darkgreen',lwd=1.5)
lines(fst.4.rus.dim$window_pos_1,fst.4.rus.dim$avg_wc_fst,col='deepskyblue3',lwd=1.5)
lines(fst.4.aet.smi$window_pos_1,fst.4.aet.smi$avg_wc_fst,col='skyblue',lwd=1.5)
lines(fst.4.neo.jav$window_pos_1,fst.4.neo.jav$avg_wc_fst,col='mediumseagreen',lwd=1.5)

plot(fst.z.sav.rus$window_pos_1,fst.z.sav.rus$avg_wc_fst,type='l',ylim=c(0,1),xlim=c(0,9.02e+07),ylab='Fst',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(fst.z.tyt.rus$window_pos_1,fst.z.tyt.rus$avg_wc_fst,col='darkorange3',lwd=1.5)
lines(fst.z.rus.aet$window_pos_1,fst.z.rus.aet$avg_wc_fst,col='plum2',lwd=1.5)
lines(fst.z.rus.smi$window_pos_1,fst.z.rus.smi$avg_wc_fst,col='goldenrod3',lwd=1.5)
lines(fst.z.rus.neo$window_pos_1,fst.z.rus.neo$avg_wc_fst,col='springgreen3',lwd=1.5)
#lines(fst.z.rus.jav$window_pos_1,fst.z.rus.jav$avg_wc_fst,col='darkgreen',lwd=1.5)
lines(fst.z.rus.dim$window_pos_1,fst.z.rus.dim$avg_wc_fst,col='deepskyblue3',lwd=1.5)
lines(fst.z.aet.smi$window_pos_1,fst.z.aet.smi$avg_wc_fst,col='skyblue',lwd=1.5)
lines(fst.z.neo.jav$window_pos_1,fst.z.neo.jav$avg_wc_fst,col='mediumseagreen',lwd=1.5)

## dxy
plot(dxy.1a.sav.rus$window_pos_1,dxy.1a.sav.rus$avg_dxy,type='l',ylim=c(0,0.0125),xlim=c(0,9.02e+07),ylab='dxy',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(dxy.1a.tyt.rus$window_pos_1,dxy.1a.tyt.rus$avg_dxy,col='darkorange3',lwd=1.5)
lines(dxy.1a.rus.aet$window_pos_1,dxy.1a.rus.aet$avg_dxy,col='plum2',lwd=1.5)
lines(dxy.1a.rus.smi$window_pos_1,dxy.1a.rus.smi$avg_dxy,col='goldenrod3',lwd=1.5)
lines(dxy.1a.rus.neo$window_pos_1,dxy.1a.rus.neo$avg_dxy,col='springgreen3',lwd=1.5)
#lines(dxy.1a.rus.jav$window_pos_1,dxy.1a.rus.jav$avg_dxy,col='darkgreen',lwd=1.5)
lines(dxy.1a.rus.dim$window_pos_1,dxy.1a.rus.dim$avg_dxy,col='deepskyblue3',lwd=1.5)
lines(dxy.1a.aet.smi$window_pos_1,dxy.1a.aet.smi$avg_dxy,col='skyblue',lwd=1.5)
lines(dxy.1a.neo.jav$window_pos_1,dxy.1a.neo.jav$avg_dxy,col='mediumseagreen',lwd=1.5)

plot(dxy.4.sav.rus$window_pos_1,dxy.4.sav.rus$avg_dxy,type='l',ylim=c(0,0.0125),xlim=c(0,9.02e+07),ylab='dxy',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(dxy.4.tyt.rus$window_pos_1,dxy.4.tyt.rus$avg_dxy,col='darkorange3',lwd=1.5)
lines(dxy.4.rus.aet$window_pos_1,dxy.4.rus.aet$avg_dxy,col='plum2',lwd=1.5)
lines(dxy.4.rus.smi$window_pos_1,dxy.4.rus.smi$avg_dxy,col='goldenrod3',lwd=1.5)
lines(dxy.4.rus.neo$window_pos_1,dxy.4.rus.neo$avg_dxy,col='springgreen3',lwd=1.5)
#lines(dxy.4.rus.jav$window_pos_1,dxy.4.rus.jav$avg_dxy,col='darkgreen',lwd=1.5)
lines(dxy.4.rus.dim$window_pos_1,dxy.4.rus.dim$avg_dxy,col='deepskyblue3',lwd=1.5)
lines(dxy.4.aet.smi$window_pos_1,dxy.4.aet.smi$avg_dxy,col='skyblue',lwd=1.5)
lines(dxy.4.neo.jav$window_pos_1,dxy.4.neo.jav$avg_dxy,col='mediumseagreen',lwd=1.5)

plot(dxy.z.sav.rus$window_pos_1,dxy.z.sav.rus$avg_dxy,type='l',ylim=c(0,0.0125),xlim=c(0,9.02e+07),ylab='dxy',xlab='Chromosome Position',col='maroon',lwd=1.5)
lines(dxy.z.tyt.rus$window_pos_1,dxy.z.tyt.rus$avg_dxy,col='darkorange3',lwd=1.5)
lines(dxy.z.rus.aet$window_pos_1,dxy.z.rus.aet$avg_dxy,col='plum2',lwd=1.5)
lines(dxy.z.rus.smi$window_pos_1,dxy.z.rus.smi$avg_dxy,col='goldenrod3',lwd=1.5)
lines(dxy.z.rus.neo$window_pos_1,dxy.z.rus.neo$avg_dxy,col='springgreen3',lwd=1.5)
#lines(dxy.z.rus.jav$window_pos_1,dxy.z.rus.jav$avg_dxy,col='darkgreen',lwd=1.5)
lines(dxy.z.rus.dim$window_pos_1,dxy.z.rus.dim$avg_dxy,col='deepskyblue3',lwd=1.5)
lines(dxy.z.aet.smi$window_pos_1,dxy.z.aet.smi$avg_dxy,col='skyblue',lwd=1.5)
lines(dxy.z.neo.jav$window_pos_1,dxy.z.neo.jav$avg_dxy,col='mediumseagreen',lwd=1.5)

## π
plot(pi.1a.rus$window_pos_1,pi.1a.rus$avg_pi,type='l',ylim=c(0,0.01),xlim=c(0,9.02e+07),ylab='π',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
lines(pi.1a.aet$window_pos_1,pi.1a.aet$avg_pi,col='plum2',lwd=1.5)
lines(pi.1a.smi$window_pos_1,pi.1a.smi$avg_pi,col='goldenrod3',lwd=1.5)
lines(pi.1a.neo$window_pos_1,pi.1a.neo$avg_pi,col='springgreen3',lwd=1.5)
lines(pi.1a.jav$window_pos_1,pi.1a.jav$avg_pi,col='darkslategray4',lwd=1.5)
lines(pi.1a.dim$window_pos_1,pi.1a.dim$avg_pi,col='deepskyblue3',lwd=1.5)

plot(pi.4.rus$window_pos_1,pi.4.rus$avg_pi,type='l',ylim=c(0,0.01),xlim=c(0,9.02e+07),ylab='π',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
lines(pi.4.aet$window_pos_1,pi.4.aet$avg_pi,col='plum2',lwd=1.5)
lines(pi.4.smi$window_pos_1,pi.4.smi$avg_pi,col='goldenrod3',lwd=1.5)
lines(pi.4.neo$window_pos_1,pi.4.neo$avg_pi,col='springgreen3',lwd=1.5)
lines(pi.4.jav$window_pos_1,pi.4.jav$avg_pi,col='darkslategray4',lwd=1.5)
lines(pi.4.dim$window_pos_1,pi.4.dim$avg_pi,col='deepskyblue3',lwd=1.5)

plot(pi.z.rus$window_pos_1,pi.z.rus$avg_pi,type='l',ylim=c(0,0.01),xlim=c(0,9.02e+07),ylab='π',xlab='Chromosome Position',col='darkorange3',lwd=1.5)
lines(pi.z.aet$window_pos_1,pi.z.aet$avg_pi,col='plum2',lwd=1.5)
lines(pi.z.smi$window_pos_1,pi.z.smi$avg_pi,col='goldenrod3',lwd=1.5)
lines(pi.z.neo$window_pos_1,pi.z.neo$avg_pi,col='springgreen3',lwd=1.5)
lines(pi.z.jav$window_pos_1,pi.z.jav$avg_pi,col='darkslategray4',lwd=1.5)
lines(pi.z.dim$window_pos_1,pi.z.dim$avg_pi,col='deepskyblue3',lwd=1.5)

## Recombination rate
plot(rho.1a$start,rho.1a$rate,type='l',col='slategray',xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',lwd=1.5)
plot(rho.4$start,rho.4$rate,type='l',col='slategray',xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',lwd=1.5)
plot(rho.z$start,rho.z$rate,type='l',col='slategray',xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',lwd=1.5)

## Output at 14 x 16 or similar aspect


