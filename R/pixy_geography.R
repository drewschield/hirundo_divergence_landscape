# Hirundo genomic divergence landscape analysis
# pixy_geography.R - compare relationships between Fst, dxy, and pi in differentiation islands for allopatric vs parapatric/partially sympatric species

# Note: we are not comparing relationships for species pairs where each n = 1, as we can't define Fst peaks between these.

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----
fst <- read.table('./results/pixy.all.1mb.fst.txt',header=T)
dxy <- read.table('./results/pixy.all.1mb.dxy.txt',header=T)
pi <- read.table('./results/pixy.all.1mb.pi.txt',header=T)

## Parse autosomes (we'll focus on autosomal regions and may also do focused Z chromosome stuff too)
fst <- fst %>% filter(chromosome != 'NC_053488.1')
dxy <- dxy %>% filter(chromosome != 'NC_053488.1')
pi <- pi %>% filter(chromosome != 'NC_053488.1')

### Parse species pairs----

fst.rus.aet <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.rus.ang <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HAN'))
fst.rus.nig <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HNI'))
fst.rus.smi <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
fst.rus.alb <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HAL'))
fst.rus.neo <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
fst.rus.jav <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
fst.rus.dim <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
fst.rus.atr <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HAT'))
fst.aet.ang <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HAN'))
fst.aet.nig <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HNI'))
fst.aet.smi <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
fst.aet.alb <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HAL'))
fst.aet.neo <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
fst.aet.jav <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
fst.aet.dim <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
fst.aet.atr <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HAT'))
#fst.ang.nig <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HNI'))
fst.ang.smi <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HS'))
#fst.ang.alb <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HAL'))
fst.ang.neo <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HN'))
fst.ang.jav <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HT'))
fst.ang.dim <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HD'))
#fst.ang.atr <- fst %>% filter((pop1 == 'HAN') & (pop2 == 'HAT'))
fst.nig.smi <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HS'))
#fst.nig.alb <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HAL'))
fst.nig.neo <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HN'))
fst.nig.jav <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HT'))
fst.nig.dim <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HD'))
#fst.nig.atr <- fst %>% filter((pop1 == 'HNI') & (pop2 == 'HAT'))
fst.smi.alb <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HAL'))
fst.smi.neo <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
fst.smi.jav <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
fst.smi.dim <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
fst.smi.atr <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HAT'))
fst.alb.neo <- fst %>% filter((pop1 == 'HAL') & (pop2 == 'HN'))
fst.alb.jav <- fst %>% filter((pop1 == 'HAL') & (pop2 == 'HT'))
fst.alb.dim <- fst %>% filter((pop1 == 'HAL') & (pop2 == 'HD'))
#fst.alb.atr <- fst %>% filter((pop1 == 'HAL') & (pop2 == 'HAT'))
fst.neo.jav <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.neo.dim <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
fst.neo.atr <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HAT'))
fst.jav.dim <- fst %>% filter((pop1 == 'HT') & (pop2 == 'HD'))
fst.jav.atr <- fst %>% filter((pop1 == 'HT') & (pop2 == 'HAT'))
fst.dim.atr <- fst %>% filter((pop1 == 'HD') & (pop2 == 'HAT'))

dxy.rus.aet <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
dxy.rus.ang <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HAN'))
dxy.rus.nig <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HNI'))
dxy.rus.smi <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
dxy.rus.alb <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HAL'))
dxy.rus.neo <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
dxy.rus.jav <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
dxy.rus.dim <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
dxy.rus.atr <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HAT'))
dxy.aet.ang <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HAN'))
dxy.aet.nig <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HNI'))
dxy.aet.smi <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
dxy.aet.alb <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HAL'))
dxy.aet.neo <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
dxy.aet.jav <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
dxy.aet.dim <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
dxy.aet.atr <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HAT'))
#dxy.ang.nig <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HNI'))
dxy.ang.smi <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HS'))
#dxy.ang.alb <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HAL'))
dxy.ang.neo <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HN'))
dxy.ang.jav <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HT'))
dxy.ang.dim <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HD'))
#dxy.ang.atr <- dxy %>% filter((pop1 == 'HAN') & (pop2 == 'HAT'))
dxy.nig.smi <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HS'))
#dxy.nig.alb <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HAL'))
dxy.nig.neo <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HN'))
dxy.nig.jav <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HT'))
dxy.nig.dim <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HD'))
#dxy.nig.atr <- dxy %>% filter((pop1 == 'HNI') & (pop2 == 'HAT'))
dxy.smi.alb <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HAL'))
dxy.smi.neo <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
dxy.smi.jav <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
dxy.smi.dim <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
dxy.smi.atr <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HAT'))
dxy.alb.neo <- dxy %>% filter((pop1 == 'HAL') & (pop2 == 'HN'))
dxy.alb.jav <- dxy %>% filter((pop1 == 'HAL') & (pop2 == 'HT'))
dxy.alb.dim <- dxy %>% filter((pop1 == 'HAL') & (pop2 == 'HD'))
#dxy.alb.atr <- dxy %>% filter((pop1 == 'HAL') & (pop2 == 'HAT'))
dxy.neo.jav <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
dxy.neo.dim <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
dxy.neo.atr <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HAT'))
dxy.jav.dim <- dxy %>% filter((pop1 == 'HT') & (pop2 == 'HD'))
dxy.jav.atr <- dxy %>% filter((pop1 == 'HT') & (pop2 == 'HAT'))
dxy.dim.atr <- dxy %>% filter((pop1 == 'HD') & (pop2 == 'HAT'))

pi.rus <- pi %>% filter(pop == 'HRR')
pi.aet <- pi %>% filter(pop == 'HA')
pi.ang <- pi %>% filter(pop == 'HAN')
pi.nig <- pi %>% filter(pop == 'HNI')
pi.alb <- pi %>% filter(pop == 'HAL')
pi.smi <- pi %>% filter(pop == 'HS')
pi.neo <- pi %>% filter(pop == 'HN')
pi.jav <- pi %>% filter(pop == 'HT')
pi.dim <- pi %>% filter(pop == 'HD')
pi.atr <- pi %>% filter(pop == 'HAT')

### Extract Fst peaks and backgrounds----
q.rus.aet <- quantile(fst.rus.aet$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.ang <- quantile(fst.rus.ang$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.nig <- quantile(fst.rus.nig$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.smi <- quantile(fst.rus.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.alb <- quantile(fst.rus.alb$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.neo <- quantile(fst.rus.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.jav <- quantile(fst.rus.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.dim <- quantile(fst.rus.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.rus.atr <- quantile(fst.rus.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.ang <- quantile(fst.aet.ang$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.nig <- quantile(fst.aet.nig$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.smi <- quantile(fst.aet.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.alb <- quantile(fst.aet.alb$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.neo <- quantile(fst.aet.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.jav <- quantile(fst.aet.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.dim <- quantile(fst.aet.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.aet.atr <- quantile(fst.aet.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.ang.smi <- quantile(fst.ang.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.ang.neo <- quantile(fst.ang.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.ang.jav <- quantile(fst.ang.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.ang.dim <- quantile(fst.ang.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.nig.smi <- quantile(fst.nig.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.nig.neo <- quantile(fst.nig.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.nig.jav <- quantile(fst.nig.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.nig.dim <- quantile(fst.nig.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.smi.alb <- quantile(fst.smi.alb$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.smi.neo <- quantile(fst.smi.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.smi.jav <- quantile(fst.smi.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.smi.dim <- quantile(fst.smi.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.smi.atr <- quantile(fst.smi.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.alb.neo <- quantile(fst.alb.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.alb.jav <- quantile(fst.alb.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.alb.dim <- quantile(fst.alb.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.neo.jav <- quantile(fst.neo.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.neo.dim <- quantile(fst.neo.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.neo.atr <- quantile(fst.neo.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.jav.dim <- quantile(fst.jav.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.jav.atr <- quantile(fst.jav.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)
q.dim.atr <- quantile(fst.dim.atr$avg_wc_fst,c(0.90),na.rm=T,names=F)

### Extract peak windows----
fst.rus.aet.is <- fst.rus.aet %>% filter(avg_wc_fst >= q.rus.aet)
fst.rus.ang.is <- fst.rus.ang %>% filter(avg_wc_fst >= q.rus.ang)
fst.rus.nig.is <- fst.rus.nig %>% filter(avg_wc_fst >= q.rus.nig)
fst.rus.smi.is <- fst.rus.smi %>% filter(avg_wc_fst >= q.rus.smi)
fst.rus.alb.is <- fst.rus.alb %>% filter(avg_wc_fst >= q.rus.alb)
fst.rus.neo.is <- fst.rus.neo %>% filter(avg_wc_fst >= q.rus.neo)
fst.rus.jav.is <- fst.rus.jav %>% filter(avg_wc_fst >= q.rus.jav)
fst.rus.dim.is <- fst.rus.dim %>% filter(avg_wc_fst >= q.rus.dim)
fst.rus.atr.is <- fst.rus.atr %>% filter(avg_wc_fst >= q.rus.atr)
fst.aet.ang.is <- fst.aet.ang %>% filter(avg_wc_fst >= q.aet.ang)
fst.aet.nig.is <- fst.aet.nig %>% filter(avg_wc_fst >= q.aet.nig)
fst.aet.smi.is <- fst.aet.smi %>% filter(avg_wc_fst >= q.aet.smi)
fst.aet.alb.is <- fst.aet.alb %>% filter(avg_wc_fst >= q.aet.alb)
fst.aet.neo.is <- fst.aet.neo %>% filter(avg_wc_fst >= q.aet.neo)
fst.aet.jav.is <- fst.aet.jav %>% filter(avg_wc_fst >= q.aet.jav)
fst.aet.dim.is <- fst.aet.dim %>% filter(avg_wc_fst >= q.aet.dim)
fst.aet.atr.is <- fst.aet.atr %>% filter(avg_wc_fst >= q.aet.atr)
fst.ang.smi.is <- fst.ang.smi %>% filter(avg_wc_fst >= q.ang.smi)
fst.ang.neo.is <- fst.ang.neo %>% filter(avg_wc_fst >= q.ang.neo)
fst.ang.jav.is <- fst.ang.jav %>% filter(avg_wc_fst >= q.ang.jav)
fst.ang.dim.is <- fst.ang.dim %>% filter(avg_wc_fst >= q.ang.dim)
fst.nig.smi.is <- fst.nig.smi %>% filter(avg_wc_fst >= q.nig.smi)
fst.nig.neo.is <- fst.nig.neo %>% filter(avg_wc_fst >= q.nig.neo)
fst.nig.jav.is <- fst.nig.jav %>% filter(avg_wc_fst >= q.nig.jav)
fst.nig.dim.is <- fst.nig.dim %>% filter(avg_wc_fst >= q.nig.dim)
fst.smi.alb.is <- fst.smi.alb %>% filter(avg_wc_fst >= q.smi.alb)
fst.smi.neo.is <- fst.smi.neo %>% filter(avg_wc_fst >= q.smi.neo)
fst.smi.jav.is <- fst.smi.jav %>% filter(avg_wc_fst >= q.smi.jav)
fst.smi.dim.is <- fst.smi.dim %>% filter(avg_wc_fst >= q.smi.dim)
fst.smi.atr.is <- fst.smi.atr %>% filter(avg_wc_fst >= q.smi.atr)
fst.alb.neo.is <- fst.alb.neo %>% filter(avg_wc_fst >= q.alb.neo)
fst.alb.jav.is <- fst.alb.jav %>% filter(avg_wc_fst >= q.alb.jav)
fst.alb.dim.is <- fst.alb.dim %>% filter(avg_wc_fst >= q.alb.dim)
fst.neo.jav.is <- fst.neo.jav %>% filter(avg_wc_fst >= q.neo.jav)
fst.neo.dim.is <- fst.neo.dim %>% filter(avg_wc_fst >= q.neo.dim)
fst.neo.atr.is <- fst.neo.atr %>% filter(avg_wc_fst >= q.neo.atr)
fst.jav.dim.is <- fst.jav.dim %>% filter(avg_wc_fst >= q.jav.dim)
fst.jav.atr.is <- fst.jav.atr %>% filter(avg_wc_fst >= q.jav.atr)
fst.dim.atr.is <- fst.dim.atr %>% filter(avg_wc_fst >= q.dim.atr)

dxy.rus.aet.is <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
dxy.rus.ang.is <- dxy.rus.ang[which(fst.rus.ang$avg_wc_fst>=q.rus.ang),]
dxy.rus.nig.is <- dxy.rus.nig[which(fst.rus.nig$avg_wc_fst>=q.rus.nig),]
dxy.rus.smi.is <- dxy.rus.smi[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
dxy.rus.alb.is <- dxy.rus.alb[which(fst.rus.alb$avg_wc_fst>=q.rus.alb),]
dxy.rus.neo.is <- dxy.rus.neo[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
dxy.rus.jav.is <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
dxy.rus.dim.is <- dxy.rus.dim[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
dxy.rus.atr.is <- dxy.rus.atr[which(fst.rus.atr$avg_wc_fst>=q.rus.atr),]
dxy.aet.ang.is <- dxy.aet.ang[which(fst.aet.ang$avg_wc_fst>=q.aet.ang),]
dxy.aet.nig.is <- dxy.aet.nig[which(fst.aet.nig$avg_wc_fst>=q.aet.nig),]
dxy.aet.smi.is <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
dxy.aet.alb.is <- dxy.aet.alb[which(fst.aet.alb$avg_wc_fst>=q.aet.alb),]
dxy.aet.neo.is <- dxy.aet.neo[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
dxy.aet.jav.is <- dxy.aet.jav[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
dxy.aet.dim.is <- dxy.aet.dim[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
dxy.aet.atr.is <- dxy.aet.atr[which(fst.aet.atr$avg_wc_fst>=q.aet.atr),]
dxy.ang.smi.is <- dxy.ang.smi[which(fst.ang.smi$avg_wc_fst>=q.ang.smi),]
dxy.ang.neo.is <- dxy.ang.neo[which(fst.ang.neo$avg_wc_fst>=q.ang.neo),]
dxy.ang.jav.is <- dxy.ang.jav[which(fst.ang.jav$avg_wc_fst>=q.ang.jav),]
dxy.ang.dim.is <- dxy.ang.dim[which(fst.ang.dim$avg_wc_fst>=q.ang.dim),]
dxy.nig.smi.is <- dxy.nig.smi[which(fst.nig.smi$avg_wc_fst>=q.nig.smi),]
dxy.nig.neo.is <- dxy.nig.neo[which(fst.nig.neo$avg_wc_fst>=q.nig.neo),]
dxy.nig.jav.is <- dxy.nig.jav[which(fst.nig.jav$avg_wc_fst>=q.nig.jav),]
dxy.nig.dim.is <- dxy.nig.dim[which(fst.nig.dim$avg_wc_fst>=q.nig.dim),]
dxy.smi.alb.is <- dxy.smi.alb[which(fst.smi.alb$avg_wc_fst>=q.smi.alb),]
dxy.smi.neo.is <- dxy.smi.neo[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
dxy.smi.jav.is <- dxy.smi.jav[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
dxy.smi.dim.is <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
dxy.smi.atr.is <- dxy.smi.atr[which(fst.smi.atr$avg_wc_fst>=q.smi.atr),]
dxy.alb.neo.is <- dxy.alb.neo[which(fst.alb.neo$avg_wc_fst>=q.alb.neo),]
dxy.alb.jav.is <- dxy.alb.jav[which(fst.alb.jav$avg_wc_fst>=q.alb.jav),]
dxy.alb.dim.is <- dxy.alb.dim[which(fst.alb.dim$avg_wc_fst>=q.alb.dim),]
dxy.neo.jav.is <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
dxy.neo.dim.is <- dxy.neo.dim[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
dxy.neo.atr.is <- dxy.neo.atr[which(fst.neo.atr$avg_wc_fst>=q.neo.atr),]
dxy.jav.dim.is <- dxy.jav.dim[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]
dxy.jav.atr.is <- dxy.jav.atr[which(fst.jav.atr$avg_wc_fst>=q.jav.atr),]
dxy.dim.atr.is <- dxy.dim.atr[which(fst.dim.atr$avg_wc_fst>=q.dim.atr),]

pi.rus.rus.aet.is <- pi.rus[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.rus.rus.ang.is <- pi.rus[which(fst.rus.ang$avg_wc_fst>=q.rus.ang),]
pi.rus.rus.nig.is <- pi.rus[which(fst.rus.nig$avg_wc_fst>=q.rus.nig),]
pi.rus.rus.smi.is <- pi.rus[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.rus.rus.alb.is <- pi.rus[which(fst.rus.alb$avg_wc_fst>=q.rus.alb),]
pi.rus.rus.neo.is <- pi.rus[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.rus.rus.jav.is <- pi.rus[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.rus.rus.dim.is <- pi.rus[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.rus.rus.atr.is <- pi.rus[which(fst.rus.atr$avg_wc_fst>=q.rus.atr),]
pi.aet.aet.ang.is <- pi.aet[which(fst.aet.ang$avg_wc_fst>=q.aet.ang),]
pi.aet.aet.nig.is <- pi.aet[which(fst.aet.nig$avg_wc_fst>=q.aet.nig),]
pi.aet.aet.smi.is <- pi.aet[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.aet.aet.alb.is <- pi.aet[which(fst.aet.alb$avg_wc_fst>=q.aet.alb),]
pi.aet.aet.neo.is <- pi.aet[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.aet.aet.jav.is <- pi.aet[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.aet.aet.dim.is <- pi.aet[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.aet.aet.atr.is <- pi.aet[which(fst.aet.atr$avg_wc_fst>=q.aet.atr),]
pi.ang.ang.smi.is <- pi.ang[which(fst.ang.smi$avg_wc_fst>=q.ang.smi),]
pi.ang.ang.neo.is <- pi.ang[which(fst.ang.neo$avg_wc_fst>=q.ang.neo),]
pi.ang.ang.jav.is <- pi.ang[which(fst.ang.jav$avg_wc_fst>=q.ang.jav),]
pi.ang.ang.dim.is <- pi.ang[which(fst.ang.dim$avg_wc_fst>=q.ang.dim),]
pi.nig.nig.smi.is <- pi.nig[which(fst.nig.smi$avg_wc_fst>=q.nig.smi),]
pi.nig.nig.neo.is <- pi.nig[which(fst.nig.neo$avg_wc_fst>=q.nig.neo),]
pi.nig.nig.jav.is <- pi.nig[which(fst.nig.jav$avg_wc_fst>=q.nig.jav),]
pi.nig.nig.dim.is <- pi.nig[which(fst.nig.dim$avg_wc_fst>=q.nig.dim),]
pi.smi.smi.alb.is <- pi.smi[which(fst.smi.alb$avg_wc_fst>=q.smi.alb),]
pi.smi.smi.neo.is <- pi.smi[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.smi.smi.jav.is <- pi.smi[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.smi.smi.dim.is <- pi.smi[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.smi.smi.atr.is <- pi.smi[which(fst.smi.atr$avg_wc_fst>=q.smi.atr),]
pi.alb.alb.neo.is <- pi.alb[which(fst.alb.neo$avg_wc_fst>=q.alb.neo),]
pi.alb.alb.jav.is <- pi.alb[which(fst.alb.jav$avg_wc_fst>=q.alb.jav),]
pi.alb.alb.dim.is <- pi.alb[which(fst.alb.dim$avg_wc_fst>=q.alb.dim),]
pi.neo.neo.jav.is <- pi.neo[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.neo.neo.dim.is <- pi.neo[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.neo.neo.atr.is <- pi.neo[which(fst.neo.atr$avg_wc_fst>=q.neo.atr),]
pi.jav.jav.dim.is <- pi.jav[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]
pi.jav.jav.atr.is <- pi.jav[which(fst.jav.atr$avg_wc_fst>=q.jav.atr),]
pi.dim.dim.atr.is <- pi.dim[which(fst.dim.atr$avg_wc_fst>=q.dim.atr),]

pi.aet.rus.aet.is <- pi.aet[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.ang.rus.ang.is <- pi.ang[which(fst.rus.ang$avg_wc_fst>=q.rus.ang),]
pi.nig.rus.nig.is <- pi.nig[which(fst.rus.nig$avg_wc_fst>=q.rus.nig),]
pi.smi.rus.smi.is <- pi.smi[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.alb.rus.alb.is <- pi.alb[which(fst.rus.alb$avg_wc_fst>=q.rus.alb),]
pi.neo.rus.neo.is <- pi.neo[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.jav.rus.jav.is <- pi.jav[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.dim.rus.dim.is <- pi.dim[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.atr.rus.atr.is <- pi.atr[which(fst.rus.atr$avg_wc_fst>=q.rus.atr),]
pi.ang.aet.ang.is <- pi.ang[which(fst.aet.ang$avg_wc_fst>=q.aet.ang),]
pi.nig.aet.nig.is <- pi.nig[which(fst.aet.nig$avg_wc_fst>=q.aet.nig),]
pi.smi.aet.smi.is <- pi.smi[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.alb.aet.alb.is <- pi.alb[which(fst.aet.alb$avg_wc_fst>=q.aet.alb),]
pi.neo.aet.neo.is <- pi.neo[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.jav.aet.jav.is <- pi.jav[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.dim.aet.dim.is <- pi.dim[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.atr.aet.atr.is <- pi.atr[which(fst.aet.atr$avg_wc_fst>=q.aet.atr),]
pi.smi.ang.smi.is <- pi.smi[which(fst.ang.smi$avg_wc_fst>=q.ang.smi),]
pi.neo.ang.neo.is <- pi.neo[which(fst.ang.neo$avg_wc_fst>=q.ang.neo),]
pi.jav.ang.jav.is <- pi.jav[which(fst.ang.jav$avg_wc_fst>=q.ang.jav),]
pi.dim.ang.dim.is <- pi.dim[which(fst.ang.dim$avg_wc_fst>=q.ang.dim),]
pi.smi.nig.smi.is <- pi.smi[which(fst.nig.smi$avg_wc_fst>=q.nig.smi),]
pi.neo.nig.neo.is <- pi.neo[which(fst.nig.neo$avg_wc_fst>=q.nig.neo),]
pi.jav.nig.jav.is <- pi.jav[which(fst.nig.jav$avg_wc_fst>=q.nig.jav),]
pi.dim.nig.dim.is <- pi.dim[which(fst.nig.dim$avg_wc_fst>=q.nig.dim),]
pi.alb.smi.alb.is <- pi.alb[which(fst.smi.alb$avg_wc_fst>=q.smi.alb),]
pi.neo.smi.neo.is <- pi.neo[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.jav.smi.jav.is <- pi.jav[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.dim.smi.dim.is <- pi.dim[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.atr.smi.atr.is <- pi.atr[which(fst.smi.atr$avg_wc_fst>=q.smi.atr),]
pi.neo.alb.neo.is <- pi.neo[which(fst.alb.neo$avg_wc_fst>=q.alb.neo),]
pi.jav.alb.jav.is <- pi.jav[which(fst.alb.jav$avg_wc_fst>=q.alb.jav),]
pi.dim.alb.dim.is <- pi.dim[which(fst.alb.dim$avg_wc_fst>=q.alb.dim),]
pi.jav.neo.jav.is <- pi.jav[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.dim.neo.dim.is <- pi.dim[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.atr.neo.atr.is <- pi.atr[which(fst.neo.atr$avg_wc_fst>=q.neo.atr),]
pi.dim.jav.dim.is <- pi.dim[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]
pi.atr.jav.atr.is <- pi.atr[which(fst.jav.atr$avg_wc_fst>=q.jav.atr),]
pi.atr.dim.atr.is <- pi.atr[which(fst.dim.atr$avg_wc_fst>=q.dim.atr),]

### Extract background windows----
fst.rus.aet.va <- fst.rus.aet %>% filter(avg_wc_fst < q.rus.aet)
fst.rus.ang.va <- fst.rus.ang %>% filter(avg_wc_fst < q.rus.ang)
fst.rus.nig.va <- fst.rus.nig %>% filter(avg_wc_fst < q.rus.nig)
fst.rus.smi.va <- fst.rus.smi %>% filter(avg_wc_fst < q.rus.smi)
fst.rus.alb.va <- fst.rus.alb %>% filter(avg_wc_fst < q.rus.alb)
fst.rus.neo.va <- fst.rus.neo %>% filter(avg_wc_fst < q.rus.neo)
fst.rus.jav.va <- fst.rus.jav %>% filter(avg_wc_fst < q.rus.jav)
fst.rus.dim.va <- fst.rus.dim %>% filter(avg_wc_fst < q.rus.dim)
fst.rus.atr.va <- fst.rus.atr %>% filter(avg_wc_fst < q.rus.atr)
fst.aet.ang.va <- fst.aet.ang %>% filter(avg_wc_fst < q.aet.ang)
fst.aet.nig.va <- fst.aet.nig %>% filter(avg_wc_fst < q.aet.nig)
fst.aet.smi.va <- fst.aet.smi %>% filter(avg_wc_fst < q.aet.smi)
fst.aet.alb.va <- fst.aet.alb %>% filter(avg_wc_fst < q.aet.alb)
fst.aet.neo.va <- fst.aet.neo %>% filter(avg_wc_fst < q.aet.neo)
fst.aet.jav.va <- fst.aet.jav %>% filter(avg_wc_fst < q.aet.jav)
fst.aet.dim.va <- fst.aet.dim %>% filter(avg_wc_fst < q.aet.dim)
fst.aet.atr.va <- fst.aet.atr %>% filter(avg_wc_fst < q.aet.atr)
fst.ang.smi.va <- fst.ang.smi %>% filter(avg_wc_fst < q.ang.smi)
fst.ang.neo.va <- fst.ang.neo %>% filter(avg_wc_fst < q.ang.neo)
fst.ang.jav.va <- fst.ang.jav %>% filter(avg_wc_fst < q.ang.jav)
fst.ang.dim.va <- fst.ang.dim %>% filter(avg_wc_fst < q.ang.dim)
fst.nig.smi.va <- fst.nig.smi %>% filter(avg_wc_fst < q.nig.smi)
fst.nig.neo.va <- fst.nig.neo %>% filter(avg_wc_fst < q.nig.neo)
fst.nig.jav.va <- fst.nig.jav %>% filter(avg_wc_fst < q.nig.jav)
fst.nig.dim.va <- fst.nig.dim %>% filter(avg_wc_fst < q.nig.dim)
fst.smi.alb.va <- fst.smi.alb %>% filter(avg_wc_fst < q.smi.alb)
fst.smi.neo.va <- fst.smi.neo %>% filter(avg_wc_fst < q.smi.neo)
fst.smi.jav.va <- fst.smi.jav %>% filter(avg_wc_fst < q.smi.jav)
fst.smi.dim.va <- fst.smi.dim %>% filter(avg_wc_fst < q.smi.dim)
fst.smi.atr.va <- fst.smi.atr %>% filter(avg_wc_fst < q.smi.atr)
fst.alb.neo.va <- fst.alb.neo %>% filter(avg_wc_fst < q.alb.neo)
fst.alb.jav.va <- fst.alb.jav %>% filter(avg_wc_fst < q.alb.jav)
fst.alb.dim.va <- fst.alb.dim %>% filter(avg_wc_fst < q.alb.dim)
fst.neo.jav.va <- fst.neo.jav %>% filter(avg_wc_fst < q.neo.jav)
fst.neo.dim.va <- fst.neo.dim %>% filter(avg_wc_fst < q.neo.dim)
fst.neo.atr.va <- fst.neo.atr %>% filter(avg_wc_fst < q.neo.atr)
fst.jav.dim.va <- fst.jav.dim %>% filter(avg_wc_fst < q.jav.dim)
fst.jav.atr.va <- fst.jav.atr %>% filter(avg_wc_fst < q.jav.atr)
fst.dim.atr.va <- fst.dim.atr %>% filter(avg_wc_fst < q.dim.atr)

dxy.rus.aet.va <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
dxy.rus.ang.va <- dxy.rus.ang[which(fst.rus.ang$avg_wc_fst<q.rus.ang),]
dxy.rus.nig.va <- dxy.rus.nig[which(fst.rus.nig$avg_wc_fst<q.rus.nig),]
dxy.rus.smi.va <- dxy.rus.smi[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
dxy.rus.alb.va <- dxy.rus.alb[which(fst.rus.alb$avg_wc_fst<q.rus.alb),]
dxy.rus.neo.va <- dxy.rus.neo[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
dxy.rus.jav.va <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
dxy.rus.dim.va <- dxy.rus.dim[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
dxy.rus.atr.va <- dxy.rus.atr[which(fst.rus.atr$avg_wc_fst<q.rus.atr),]
dxy.aet.ang.va <- dxy.aet.ang[which(fst.aet.ang$avg_wc_fst<q.aet.ang),]
dxy.aet.nig.va <- dxy.aet.nig[which(fst.aet.nig$avg_wc_fst<q.aet.nig),]
dxy.aet.smi.va <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
dxy.aet.alb.va <- dxy.aet.alb[which(fst.aet.alb$avg_wc_fst<q.aet.alb),]
dxy.aet.neo.va <- dxy.aet.neo[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
dxy.aet.jav.va <- dxy.aet.jav[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
dxy.aet.dim.va <- dxy.aet.dim[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
dxy.aet.atr.va <- dxy.aet.atr[which(fst.aet.atr$avg_wc_fst<q.aet.atr),]
dxy.ang.smi.va <- dxy.ang.smi[which(fst.ang.smi$avg_wc_fst<q.ang.smi),]
dxy.ang.neo.va <- dxy.ang.neo[which(fst.ang.neo$avg_wc_fst<q.ang.neo),]
dxy.ang.jav.va <- dxy.ang.jav[which(fst.ang.jav$avg_wc_fst<q.ang.jav),]
dxy.ang.dim.va <- dxy.ang.dim[which(fst.ang.dim$avg_wc_fst<q.ang.dim),]
dxy.nig.smi.va <- dxy.nig.smi[which(fst.nig.smi$avg_wc_fst<q.nig.smi),]
dxy.nig.neo.va <- dxy.nig.neo[which(fst.nig.neo$avg_wc_fst<q.nig.neo),]
dxy.nig.jav.va <- dxy.nig.jav[which(fst.nig.jav$avg_wc_fst<q.nig.jav),]
dxy.nig.dim.va <- dxy.nig.dim[which(fst.nig.dim$avg_wc_fst<q.nig.dim),]
dxy.smi.alb.va <- dxy.smi.alb[which(fst.smi.alb$avg_wc_fst<q.smi.alb),]
dxy.smi.neo.va <- dxy.smi.neo[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
dxy.smi.jav.va <- dxy.smi.jav[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
dxy.smi.dim.va <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
dxy.smi.atr.va <- dxy.smi.atr[which(fst.smi.atr$avg_wc_fst<q.smi.atr),]
dxy.alb.neo.va <- dxy.alb.neo[which(fst.alb.neo$avg_wc_fst<q.alb.neo),]
dxy.alb.jav.va <- dxy.alb.jav[which(fst.alb.jav$avg_wc_fst<q.alb.jav),]
dxy.alb.dim.va <- dxy.alb.dim[which(fst.alb.dim$avg_wc_fst<q.alb.dim),]
dxy.neo.jav.va <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
dxy.neo.dim.va <- dxy.neo.dim[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
dxy.neo.atr.va <- dxy.neo.atr[which(fst.neo.atr$avg_wc_fst<q.neo.atr),]
dxy.jav.dim.va <- dxy.jav.dim[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]
dxy.jav.atr.va <- dxy.jav.atr[which(fst.jav.atr$avg_wc_fst<q.jav.atr),]
dxy.dim.atr.va <- dxy.dim.atr[which(fst.dim.atr$avg_wc_fst<q.dim.atr),]

pi.rus.rus.aet.va <- pi.rus[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
pi.rus.rus.ang.va <- pi.rus[which(fst.rus.ang$avg_wc_fst<q.rus.ang),]
pi.rus.rus.nig.va <- pi.rus[which(fst.rus.nig$avg_wc_fst<q.rus.nig),]
pi.rus.rus.smi.va <- pi.rus[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
pi.rus.rus.alb.va <- pi.rus[which(fst.rus.alb$avg_wc_fst<q.rus.alb),]
pi.rus.rus.neo.va <- pi.rus[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
pi.rus.rus.jav.va <- pi.rus[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
pi.rus.rus.dim.va <- pi.rus[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
pi.rus.rus.atr.va <- pi.rus[which(fst.rus.atr$avg_wc_fst<q.rus.atr),]
pi.aet.aet.ang.va <- pi.aet[which(fst.aet.ang$avg_wc_fst<q.aet.ang),]
pi.aet.aet.nig.va <- pi.aet[which(fst.aet.nig$avg_wc_fst<q.aet.nig),]
pi.aet.aet.smi.va <- pi.aet[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
pi.aet.aet.alb.va <- pi.aet[which(fst.aet.alb$avg_wc_fst<q.aet.alb),]
pi.aet.aet.neo.va <- pi.aet[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
pi.aet.aet.jav.va <- pi.aet[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
pi.aet.aet.dim.va <- pi.aet[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
pi.aet.aet.atr.va <- pi.aet[which(fst.aet.atr$avg_wc_fst<q.aet.atr),]
pi.ang.ang.smi.va <- pi.ang[which(fst.ang.smi$avg_wc_fst<q.ang.smi),]
pi.ang.ang.neo.va <- pi.ang[which(fst.ang.neo$avg_wc_fst<q.ang.neo),]
pi.ang.ang.jav.va <- pi.ang[which(fst.ang.jav$avg_wc_fst<q.ang.jav),]
pi.ang.ang.dim.va <- pi.ang[which(fst.ang.dim$avg_wc_fst<q.ang.dim),]
pi.nig.nig.smi.va <- pi.nig[which(fst.nig.smi$avg_wc_fst<q.nig.smi),]
pi.nig.nig.neo.va <- pi.nig[which(fst.nig.neo$avg_wc_fst<q.nig.neo),]
pi.nig.nig.jav.va <- pi.nig[which(fst.nig.jav$avg_wc_fst<q.nig.jav),]
pi.nig.nig.dim.va <- pi.nig[which(fst.nig.dim$avg_wc_fst<q.nig.dim),]
pi.smi.smi.alb.va <- pi.smi[which(fst.smi.alb$avg_wc_fst<q.smi.alb),]
pi.smi.smi.neo.va <- pi.smi[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
pi.smi.smi.jav.va <- pi.smi[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
pi.smi.smi.dim.va <- pi.smi[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
pi.smi.smi.atr.va <- pi.smi[which(fst.smi.atr$avg_wc_fst<q.smi.atr),]
pi.alb.alb.neo.va <- pi.alb[which(fst.alb.neo$avg_wc_fst<q.alb.neo),]
pi.alb.alb.jav.va <- pi.alb[which(fst.alb.jav$avg_wc_fst<q.alb.jav),]
pi.alb.alb.dim.va <- pi.alb[which(fst.alb.dim$avg_wc_fst<q.alb.dim),]
pi.neo.neo.jav.va <- pi.neo[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
pi.neo.neo.dim.va <- pi.neo[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
pi.neo.neo.atr.va <- pi.neo[which(fst.neo.atr$avg_wc_fst<q.neo.atr),]
pi.jav.jav.dim.va <- pi.jav[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]
pi.jav.jav.atr.va <- pi.jav[which(fst.jav.atr$avg_wc_fst<q.jav.atr),]
pi.dim.dim.atr.va <- pi.dim[which(fst.dim.atr$avg_wc_fst<q.dim.atr),]

pi.aet.rus.aet.va <- pi.aet[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
pi.ang.rus.ang.va <- pi.ang[which(fst.rus.ang$avg_wc_fst<q.rus.ang),]
pi.nig.rus.nig.va <- pi.nig[which(fst.rus.nig$avg_wc_fst<q.rus.nig),]
pi.smi.rus.smi.va <- pi.smi[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
pi.alb.rus.alb.va <- pi.alb[which(fst.rus.alb$avg_wc_fst<q.rus.alb),]
pi.neo.rus.neo.va <- pi.neo[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
pi.jav.rus.jav.va <- pi.jav[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
pi.dim.rus.dim.va <- pi.dim[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
pi.atr.rus.atr.va <- pi.atr[which(fst.rus.atr$avg_wc_fst<q.rus.atr),]
pi.ang.aet.ang.va <- pi.ang[which(fst.aet.ang$avg_wc_fst<q.aet.ang),]
pi.nig.aet.nig.va <- pi.nig[which(fst.aet.nig$avg_wc_fst<q.aet.nig),]
pi.smi.aet.smi.va <- pi.smi[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
pi.alb.aet.alb.va <- pi.alb[which(fst.aet.alb$avg_wc_fst<q.aet.alb),]
pi.neo.aet.neo.va <- pi.neo[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
pi.jav.aet.jav.va <- pi.jav[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
pi.dim.aet.dim.va <- pi.dim[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
pi.atr.aet.atr.va <- pi.atr[which(fst.aet.atr$avg_wc_fst<q.aet.atr),]
pi.smi.ang.smi.va <- pi.smi[which(fst.ang.smi$avg_wc_fst<q.ang.smi),]
pi.neo.ang.neo.va <- pi.neo[which(fst.ang.neo$avg_wc_fst<q.ang.neo),]
pi.jav.ang.jav.va <- pi.jav[which(fst.ang.jav$avg_wc_fst<q.ang.jav),]
pi.dim.ang.dim.va <- pi.dim[which(fst.ang.dim$avg_wc_fst<q.ang.dim),]
pi.smi.nig.smi.va <- pi.smi[which(fst.nig.smi$avg_wc_fst<q.nig.smi),]
pi.neo.nig.neo.va <- pi.neo[which(fst.nig.neo$avg_wc_fst<q.nig.neo),]
pi.jav.nig.jav.va <- pi.jav[which(fst.nig.jav$avg_wc_fst<q.nig.jav),]
pi.dim.nig.dim.va <- pi.dim[which(fst.nig.dim$avg_wc_fst<q.nig.dim),]
pi.alb.smi.alb.va <- pi.alb[which(fst.smi.alb$avg_wc_fst<q.smi.alb),]
pi.neo.smi.neo.va <- pi.neo[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
pi.jav.smi.jav.va <- pi.jav[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
pi.dim.smi.dim.va <- pi.dim[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
pi.atr.smi.atr.va <- pi.atr[which(fst.smi.atr$avg_wc_fst<q.smi.atr),]
pi.neo.alb.neo.va <- pi.neo[which(fst.alb.neo$avg_wc_fst<q.alb.neo),]
pi.jav.alb.jav.va <- pi.jav[which(fst.alb.jav$avg_wc_fst<q.alb.jav),]
pi.dim.alb.dim.va <- pi.dim[which(fst.alb.dim$avg_wc_fst<q.alb.dim),]
pi.jav.neo.jav.va <- pi.jav[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
pi.dim.neo.dim.va <- pi.dim[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
pi.atr.neo.atr.va <- pi.atr[which(fst.neo.atr$avg_wc_fst<q.neo.atr),]
pi.dim.jav.dim.va <- pi.dim[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]
pi.atr.jav.atr.va <- pi.atr[which(fst.jav.atr$avg_wc_fst<q.jav.atr),]
pi.atr.dim.atr.va <- pi.atr[which(fst.dim.atr$avg_wc_fst<q.dim.atr),]

### Calculate mean pi for peaks and backgrounds----

pi.rus.aet.is <- (pi.rus.rus.aet.is$avg_pi + pi.aet.rus.aet.is$avg_pi)/2
pi.rus.ang.is <- (pi.rus.rus.ang.is$avg_pi + pi.ang.rus.ang.is$avg_pi)/2
pi.rus.nig.is <- (pi.rus.rus.nig.is$avg_pi + pi.nig.rus.nig.is$avg_pi)/2
pi.rus.smi.is <- (pi.rus.rus.smi.is$avg_pi + pi.smi.rus.smi.is$avg_pi)/2
pi.rus.alb.is <- (pi.rus.rus.alb.is$avg_pi + pi.alb.rus.alb.is$avg_pi)/2
pi.rus.neo.is <- (pi.rus.rus.neo.is$avg_pi + pi.neo.rus.neo.is$avg_pi)/2
pi.rus.jav.is <- (pi.rus.rus.jav.is$avg_pi + pi.jav.rus.jav.is$avg_pi)/2
pi.rus.dim.is <- (pi.rus.rus.dim.is$avg_pi + pi.dim.rus.dim.is$avg_pi)/2
pi.rus.atr.is <- (pi.rus.rus.atr.is$avg_pi + pi.atr.rus.atr.is$avg_pi)/2
pi.aet.ang.is <- (pi.aet.aet.ang.is$avg_pi + pi.ang.aet.ang.is$avg_pi)/2
pi.aet.nig.is <- (pi.aet.aet.nig.is$avg_pi + pi.nig.aet.nig.is$avg_pi)/2
pi.aet.smi.is <- (pi.aet.aet.smi.is$avg_pi + pi.smi.aet.smi.is$avg_pi)/2
pi.aet.alb.is <- (pi.aet.aet.alb.is$avg_pi + pi.alb.aet.alb.is$avg_pi)/2
pi.aet.neo.is <- (pi.aet.aet.neo.is$avg_pi + pi.neo.aet.neo.is$avg_pi)/2
pi.aet.jav.is <- (pi.aet.aet.jav.is$avg_pi + pi.jav.aet.jav.is$avg_pi)/2
pi.aet.dim.is <- (pi.aet.aet.dim.is$avg_pi + pi.dim.aet.dim.is$avg_pi)/2
pi.aet.atr.is <- (pi.aet.aet.atr.is$avg_pi + pi.atr.aet.atr.is$avg_pi)/2
pi.ang.smi.is <- (pi.ang.ang.smi.is$avg_pi + pi.smi.ang.smi.is$avg_pi)/2
pi.ang.neo.is <- (pi.ang.ang.neo.is$avg_pi + pi.neo.ang.neo.is$avg_pi)/2
pi.ang.jav.is <- (pi.ang.ang.jav.is$avg_pi + pi.jav.ang.jav.is$avg_pi)/2
pi.ang.dim.is <- (pi.ang.ang.dim.is$avg_pi + pi.dim.ang.dim.is$avg_pi)/2
pi.nig.smi.is <- (pi.nig.nig.smi.is$avg_pi + pi.smi.nig.smi.is$avg_pi)/2
pi.nig.neo.is <- (pi.nig.nig.neo.is$avg_pi + pi.neo.nig.neo.is$avg_pi)/2
pi.nig.jav.is <- (pi.nig.nig.jav.is$avg_pi + pi.jav.nig.jav.is$avg_pi)/2
pi.nig.dim.is <- (pi.nig.nig.dim.is$avg_pi + pi.dim.nig.dim.is$avg_pi)/2
pi.smi.alb.is <- (pi.smi.smi.alb.is$avg_pi + pi.alb.smi.alb.is$avg_pi)/2
pi.smi.neo.is <- (pi.smi.smi.neo.is$avg_pi + pi.neo.smi.neo.is$avg_pi)/2
pi.smi.jav.is <- (pi.smi.smi.jav.is$avg_pi + pi.jav.smi.jav.is$avg_pi)/2
pi.smi.dim.is <- (pi.smi.smi.dim.is$avg_pi + pi.dim.smi.dim.is$avg_pi)/2
pi.smi.atr.is <- (pi.smi.smi.atr.is$avg_pi + pi.atr.smi.atr.is$avg_pi)/2
pi.alb.neo.is <- (pi.alb.alb.neo.is$avg_pi + pi.neo.alb.neo.is$avg_pi)/2
pi.alb.jav.is <- (pi.alb.alb.jav.is$avg_pi + pi.jav.alb.jav.is$avg_pi)/2
pi.alb.dim.is <- (pi.alb.alb.dim.is$avg_pi + pi.dim.alb.dim.is$avg_pi)/2
pi.neo.jav.is <- (pi.neo.neo.jav.is$avg_pi + pi.jav.neo.jav.is$avg_pi)/2
pi.neo.dim.is <- (pi.neo.neo.dim.is$avg_pi + pi.dim.neo.dim.is$avg_pi)/2
pi.neo.atr.is <- (pi.neo.neo.atr.is$avg_pi + pi.atr.neo.atr.is$avg_pi)/2
pi.jav.dim.is <- (pi.jav.jav.dim.is$avg_pi + pi.dim.jav.dim.is$avg_pi)/2
pi.jav.atr.is <- (pi.jav.jav.atr.is$avg_pi + pi.atr.jav.atr.is$avg_pi)/2
pi.dim.atr.is <- (pi.dim.dim.atr.is$avg_pi + pi.atr.dim.atr.is$avg_pi)/2

pi.rus.aet.va <- (pi.rus.rus.aet.va$avg_pi + pi.aet.rus.aet.va$avg_pi)/2
pi.rus.ang.va <- (pi.rus.rus.ang.va$avg_pi + pi.ang.rus.ang.va$avg_pi)/2
pi.rus.nig.va <- (pi.rus.rus.nig.va$avg_pi + pi.nig.rus.nig.va$avg_pi)/2
pi.rus.smi.va <- (pi.rus.rus.smi.va$avg_pi + pi.smi.rus.smi.va$avg_pi)/2
pi.rus.alb.va <- (pi.rus.rus.alb.va$avg_pi + pi.alb.rus.alb.va$avg_pi)/2
pi.rus.neo.va <- (pi.rus.rus.neo.va$avg_pi + pi.neo.rus.neo.va$avg_pi)/2
pi.rus.jav.va <- (pi.rus.rus.jav.va$avg_pi + pi.jav.rus.jav.va$avg_pi)/2
pi.rus.dim.va <- (pi.rus.rus.dim.va$avg_pi + pi.dim.rus.dim.va$avg_pi)/2
pi.rus.atr.va <- (pi.rus.rus.atr.va$avg_pi + pi.atr.rus.atr.va$avg_pi)/2
pi.aet.ang.va <- (pi.aet.aet.ang.va$avg_pi + pi.ang.aet.ang.va$avg_pi)/2
pi.aet.nig.va <- (pi.aet.aet.nig.va$avg_pi + pi.nig.aet.nig.va$avg_pi)/2
pi.aet.smi.va <- (pi.aet.aet.smi.va$avg_pi + pi.smi.aet.smi.va$avg_pi)/2
pi.aet.alb.va <- (pi.aet.aet.alb.va$avg_pi + pi.alb.aet.alb.va$avg_pi)/2
pi.aet.neo.va <- (pi.aet.aet.neo.va$avg_pi + pi.neo.aet.neo.va$avg_pi)/2
pi.aet.jav.va <- (pi.aet.aet.jav.va$avg_pi + pi.jav.aet.jav.va$avg_pi)/2
pi.aet.dim.va <- (pi.aet.aet.dim.va$avg_pi + pi.dim.aet.dim.va$avg_pi)/2
pi.aet.atr.va <- (pi.aet.aet.atr.va$avg_pi + pi.atr.aet.atr.va$avg_pi)/2
pi.ang.smi.va <- (pi.ang.ang.smi.va$avg_pi + pi.smi.ang.smi.va$avg_pi)/2
pi.ang.neo.va <- (pi.ang.ang.neo.va$avg_pi + pi.neo.ang.neo.va$avg_pi)/2
pi.ang.jav.va <- (pi.ang.ang.jav.va$avg_pi + pi.jav.ang.jav.va$avg_pi)/2
pi.ang.dim.va <- (pi.ang.ang.dim.va$avg_pi + pi.dim.ang.dim.va$avg_pi)/2
pi.nig.smi.va <- (pi.nig.nig.smi.va$avg_pi + pi.smi.nig.smi.va$avg_pi)/2
pi.nig.neo.va <- (pi.nig.nig.neo.va$avg_pi + pi.neo.nig.neo.va$avg_pi)/2
pi.nig.jav.va <- (pi.nig.nig.jav.va$avg_pi + pi.jav.nig.jav.va$avg_pi)/2
pi.nig.dim.va <- (pi.nig.nig.dim.va$avg_pi + pi.dim.nig.dim.va$avg_pi)/2
pi.smi.alb.va <- (pi.smi.smi.alb.va$avg_pi + pi.alb.smi.alb.va$avg_pi)/2
pi.smi.neo.va <- (pi.smi.smi.neo.va$avg_pi + pi.neo.smi.neo.va$avg_pi)/2
pi.smi.jav.va <- (pi.smi.smi.jav.va$avg_pi + pi.jav.smi.jav.va$avg_pi)/2
pi.smi.dim.va <- (pi.smi.smi.dim.va$avg_pi + pi.dim.smi.dim.va$avg_pi)/2
pi.smi.atr.va <- (pi.smi.smi.atr.va$avg_pi + pi.atr.smi.atr.va$avg_pi)/2
pi.alb.neo.va <- (pi.alb.alb.neo.va$avg_pi + pi.neo.alb.neo.va$avg_pi)/2
pi.alb.jav.va <- (pi.alb.alb.jav.va$avg_pi + pi.jav.alb.jav.va$avg_pi)/2
pi.alb.dim.va <- (pi.alb.alb.dim.va$avg_pi + pi.dim.alb.dim.va$avg_pi)/2
pi.neo.jav.va <- (pi.neo.neo.jav.va$avg_pi + pi.jav.neo.jav.va$avg_pi)/2
pi.neo.dim.va <- (pi.neo.neo.dim.va$avg_pi + pi.dim.neo.dim.va$avg_pi)/2
pi.neo.atr.va <- (pi.neo.neo.atr.va$avg_pi + pi.atr.neo.atr.va$avg_pi)/2
pi.jav.dim.va <- (pi.jav.jav.dim.va$avg_pi + pi.dim.jav.dim.va$avg_pi)/2
pi.jav.atr.va <- (pi.jav.jav.atr.va$avg_pi + pi.atr.jav.atr.va$avg_pi)/2
pi.dim.atr.va <- (pi.dim.dim.atr.va$avg_pi + pi.atr.dim.atr.va$avg_pi)/2

### Assemble correlations for allopatric and parapatric/partially sympatric pairs

## Fst vs dxy
allo.fst.dxy.back <- c()
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.aet.va$avg_wc_fst,dxy.rus.aet.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.ang.va$avg_wc_fst,dxy.rus.ang.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.alb.va$avg_wc_fst,dxy.rus.alb.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.nig.va$avg_wc_fst,dxy.rus.nig.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.neo.va$avg_wc_fst,dxy.rus.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.dim.va$avg_wc_fst,dxy.rus.dim.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.rus.atr.va$avg_wc_fst,dxy.rus.atr.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.aet.alb.va$avg_wc_fst,dxy.aet.alb.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.aet.neo.va$avg_wc_fst,dxy.aet.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.aet.jav.va$avg_wc_fst,dxy.aet.jav.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.aet.dim.va$avg_wc_fst,dxy.aet.dim.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.aet.atr.va$avg_wc_fst,dxy.aet.atr.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.ang.alb.va$avg_wc_fst,dxy.ang.alb.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.ang.neo.va$avg_wc_fst,dxy.ang.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.ang.jav.va$avg_wc_fst,dxy.ang.jav.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.alb.nig.va$avg_wc_fst,dxy.alb.nig.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.alb.neo.va$avg_wc_fst,dxy.alb.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.alb.jav.va$avg_wc_fst,dxy.alb.jav.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.nig.neo.va$avg_wc_fst,dxy.nig.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.nig.jav.va$avg_wc_fst,dxy.nig.jav.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.nig.atr.va$avg_wc_fst,dxy.nig.atr.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.smi.neo.va$avg_wc_fst,dxy.smi.neo.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.neo.dim.va$avg_wc_fst,dxy.neo.dim.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.neo.atr.va$avg_wc_fst,dxy.neo.atr.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.jav.dim.va$avg_wc_fst,dxy.jav.dim.va$avg_dxy,method='spearman'))
allo.fst.dxy.back <- append(allo.fst.dxy.back,cor(fst.jav.atr.va$avg_wc_fst,dxy.jav.atr.va$avg_dxy,method='spearman'))

allo.fst.dxy.peak <- c()
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.aet.is$avg_wc_fst,dxy.rus.aet.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.ang.is$avg_wc_fst,dxy.rus.ang.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.alb.is$avg_wc_fst,dxy.rus.alb.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.nig.is$avg_wc_fst,dxy.rus.nig.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.neo.is$avg_wc_fst,dxy.rus.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.dim.is$avg_wc_fst,dxy.rus.dim.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.rus.atr.is$avg_wc_fst,dxy.rus.atr.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.aet.alb.is$avg_wc_fst,dxy.aet.alb.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.aet.neo.is$avg_wc_fst,dxy.aet.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.aet.jav.is$avg_wc_fst,dxy.aet.jav.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.aet.dim.is$avg_wc_fst,dxy.aet.dim.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.aet.atr.is$avg_wc_fst,dxy.aet.atr.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.ang.alb.is$avg_wc_fst,dxy.ang.alb.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.ang.neo.is$avg_wc_fst,dxy.ang.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.ang.jav.is$avg_wc_fst,dxy.ang.jav.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.alb.nig.is$avg_wc_fst,dxy.alb.nig.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.alb.neo.is$avg_wc_fst,dxy.alb.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.alb.jav.is$avg_wc_fst,dxy.alb.jav.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.nig.neo.is$avg_wc_fst,dxy.nig.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.nig.jav.is$avg_wc_fst,dxy.nig.jav.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.nig.atr.is$avg_wc_fst,dxy.nig.atr.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.smi.neo.is$avg_wc_fst,dxy.smi.neo.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.neo.dim.is$avg_wc_fst,dxy.neo.dim.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.neo.atr.is$avg_wc_fst,dxy.neo.atr.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.jav.dim.is$avg_wc_fst,dxy.jav.dim.is$avg_dxy,method='spearman'))
allo.fst.dxy.peak <- append(allo.fst.dxy.peak,cor(fst.jav.atr.is$avg_wc_fst,dxy.jav.atr.is$avg_dxy,method='spearman'))

para.fst.dxy.back <- c()
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.rus.smi.va$avg_wc_fst,dxy.rus.smi.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.rus.jav.va$avg_wc_fst,dxy.rus.jav.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.aet.ang.va$avg_wc_fst,dxy.aet.ang.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.aet.nig.va$avg_wc_fst,dxy.aet.nig.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.aet.smi.va$avg_wc_fst,dxy.aet.smi.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.ang.nig.va$avg_wc_fst,dxy.ang.nig.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.ang.smi.va$avg_wc_fst,dxy.ang.smi.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.ang.dim.va$avg_wc_fst,dxy.ang.dim.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.ang.atr.va$avg_wc_fst,dxy.ang.atr.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.alb.smi.va$avg_wc_fst,dxy.alb.smi.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.alb.dim.va$avg_wc_fst,dxy.alb.dim.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.alb.atr.va$avg_wc_fst,dxy.alb.atr.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.nig.smi.va$avg_wc_fst,dxy.nig.smi.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.nig.dim.va$avg_wc_fst,dxy.nig.dim.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.smi.jav.va$avg_wc_fst,dxy.smi.jav.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.smi.dim.va$avg_wc_fst,dxy.smi.dim.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.smi.atr.va$avg_wc_fst,dxy.smi.atr.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.neo.jav.va$avg_wc_fst,dxy.neo.jav.va$avg_dxy,method='spearman'))
para.fst.dxy.back <- append(para.fst.dxy.back,cor(fst.dim.atr.va$avg_wc_fst,dxy.dim.atr.va$avg_dxy,method='spearman'))

para.fst.dxy.peak <- c()
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.rus.smi.is$avg_wc_fst,dxy.rus.smi.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.rus.jav.is$avg_wc_fst,dxy.rus.jav.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.aet.ang.is$avg_wc_fst,dxy.aet.ang.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.aet.nig.is$avg_wc_fst,dxy.aet.nig.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.aet.smi.is$avg_wc_fst,dxy.aet.smi.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.ang.nig.is$avg_wc_fst,dxy.ang.nig.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.ang.smi.is$avg_wc_fst,dxy.ang.smi.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.ang.dim.is$avg_wc_fst,dxy.ang.dim.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.ang.atr.is$avg_wc_fst,dxy.ang.atr.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.alb.smi.is$avg_wc_fst,dxy.alb.smi.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.alb.dim.is$avg_wc_fst,dxy.alb.dim.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.alb.atr.is$avg_wc_fst,dxy.alb.atr.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.nig.smi.is$avg_wc_fst,dxy.nig.smi.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.nig.dim.is$avg_wc_fst,dxy.nig.dim.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.smi.jav.is$avg_wc_fst,dxy.smi.jav.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.smi.dim.is$avg_wc_fst,dxy.smi.dim.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.smi.atr.is$avg_wc_fst,dxy.smi.atr.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.neo.jav.is$avg_wc_fst,dxy.neo.jav.is$avg_dxy,method='spearman'))
para.fst.dxy.peak <- append(para.fst.dxy.peak,cor(fst.dim.atr.is$avg_wc_fst,dxy.dim.atr.is$avg_dxy,method='spearman'))

par(mfrow=c(1,2))
boxplot(allo.fst.dxy.back,allo.fst.dxy.peak)
boxplot(para.fst.dxy.back,para.fst.dxy.peak)

wilcox.test(allo.fst.dxy.back,allo.fst.dxy.peak)
wilcox.test(para.fst.dxy.back,para.fst.dxy.peak)

## dxy vs pi 
allo.dxy.pi.back <- c()
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.aet.va$avg_dxy,pi.rus.aet.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.ang.va$avg_dxy,pi.rus.ang.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.alb.va$avg_dxy,pi.rus.alb.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.nig.va$avg_dxy,pi.rus.nig.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.neo.va$avg_dxy,pi.rus.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.dim.va$avg_dxy,pi.rus.dim.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.rus.atr.va$avg_dxy,pi.rus.atr.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.aet.alb.va$avg_dxy,pi.aet.alb.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.aet.neo.va$avg_dxy,pi.aet.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.aet.jav.va$avg_dxy,pi.aet.jav.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.aet.dim.va$avg_dxy,pi.aet.dim.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.aet.atr.va$avg_dxy,pi.aet.atr.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.ang.alb.va$avg_dxy,pi.ang.alb.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.ang.neo.va$avg_dxy,pi.ang.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.ang.jav.va$avg_dxy,pi.ang.jav.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.alb.nig.va$avg_dxy,pi.alb.nig.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.alb.neo.va$avg_dxy,pi.alb.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.alb.jav.va$avg_dxy,pi.alb.jav.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.nig.neo.va$avg_dxy,pi.nig.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.nig.jav.va$avg_dxy,pi.nig.jav.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.nig.atr.va$avg_dxy,pi.nig.atr.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.smi.neo.va$avg_dxy,pi.smi.neo.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.neo.dim.va$avg_dxy,pi.neo.dim.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.neo.atr.va$avg_dxy,pi.neo.atr.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.jav.dim.va$avg_dxy,pi.jav.dim.va,method='spearman'))
allo.dxy.pi.back <- append(allo.dxy.pi.back,cor(dxy.jav.atr.va$avg_dxy,pi.jav.atr.va,method='spearman'))

allo.dxy.pi.peak <- c()
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.aet.is$avg_dxy,pi.rus.aet.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.ang.is$avg_dxy,pi.rus.ang.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.alb.is$avg_dxy,pi.rus.alb.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.nig.is$avg_dxy,pi.rus.nig.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.neo.is$avg_dxy,pi.rus.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.dim.is$avg_dxy,pi.rus.dim.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.rus.atr.is$avg_dxy,pi.rus.atr.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.aet.alb.is$avg_dxy,pi.aet.alb.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.aet.neo.is$avg_dxy,pi.aet.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.aet.jav.is$avg_dxy,pi.aet.jav.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.aet.dim.is$avg_dxy,pi.aet.dim.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.aet.atr.is$avg_dxy,pi.aet.atr.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.ang.alb.is$avg_dxy,pi.ang.alb.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.ang.neo.is$avg_dxy,pi.ang.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.ang.jav.is$avg_dxy,pi.ang.jav.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.alb.nig.is$avg_dxy,pi.alb.nig.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.alb.neo.is$avg_dxy,pi.alb.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.alb.jav.is$avg_dxy,pi.alb.jav.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.nig.neo.is$avg_dxy,pi.nig.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.nig.jav.is$avg_dxy,pi.nig.jav.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.nig.atr.is$avg_dxy,pi.nig.atr.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.smi.neo.is$avg_dxy,pi.smi.neo.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.neo.dim.is$avg_dxy,pi.neo.dim.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.neo.atr.is$avg_dxy,pi.neo.atr.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.jav.dim.is$avg_dxy,pi.jav.dim.is,method='spearman'))
allo.dxy.pi.peak <- append(allo.dxy.pi.peak,cor(dxy.jav.atr.is$avg_dxy,pi.jav.atr.is,method='spearman'))

para.dxy.pi.back <- c()
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.rus.smi.va$avg_dxy,pi.rus.smi.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.rus.jav.va$avg_dxy,pi.rus.jav.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.aet.ang.va$avg_dxy,pi.aet.ang.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.aet.nig.va$avg_dxy,pi.aet.nig.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.aet.smi.va$avg_dxy,pi.aet.smi.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.ang.nig.va$avg_dxy,pi.ang.nig.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.ang.smi.va$avg_dxy,pi.ang.smi.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.ang.dim.va$avg_dxy,pi.ang.dim.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.ang.atr.va$avg_dxy,pi.ang.atr.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.alb.smi.va$avg_dxy,pi.alb.smi.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.alb.dim.va$avg_dxy,pi.alb.dim.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.alb.atr.va$avg_dxy,pi.alb.atr.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.nig.smi.va$avg_dxy,pi.nig.smi.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.nig.dim.va$avg_dxy,pi.nig.dim.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.smi.jav.va$avg_dxy,pi.smi.jav.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.smi.dim.va$avg_dxy,pi.smi.dim.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.smi.atr.va$avg_dxy,pi.smi.atr.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.neo.jav.va$avg_dxy,pi.neo.jav.va,method='spearman'))
para.dxy.pi.back <- append(para.dxy.pi.back,cor(dxy.dim.atr.va$avg_dxy,pi.dim.atr.va,method='spearman'))

para.dxy.pi.peak <- c()
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.rus.smi.is$avg_dxy,pi.rus.smi.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.rus.jav.is$avg_dxy,pi.rus.jav.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.aet.ang.is$avg_dxy,pi.aet.ang.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.aet.nig.is$avg_dxy,pi.aet.nig.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.aet.smi.is$avg_dxy,pi.aet.smi.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.ang.nig.is$avg_dxy,pi.ang.nig.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.ang.smi.is$avg_dxy,pi.ang.smi.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.ang.dim.is$avg_dxy,pi.ang.dim.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.ang.atr.is$avg_dxy,pi.ang.atr.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.alb.smi.is$avg_dxy,pi.alb.smi.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.alb.dim.is$avg_dxy,pi.alb.dim.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.alb.atr.is$avg_dxy,pi.alb.atr.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.nig.smi.is$avg_dxy,pi.nig.smi.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.nig.dim.is$avg_dxy,pi.nig.dim.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.smi.jav.is$avg_dxy,pi.smi.jav.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.smi.dim.is$avg_dxy,pi.smi.dim.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.smi.atr.is$avg_dxy,pi.smi.atr.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.neo.jav.is$avg_dxy,pi.neo.jav.is,method='spearman'))
para.dxy.pi.peak <- append(para.dxy.pi.peak,cor(dxy.dim.atr.is$avg_dxy,pi.dim.atr.is,method='spearman'))

par(mfrow=c(2,1))
#boxplot(allo.fst.dxy.back,allo.fst.dxy.peak,outline=F,names=c('back','peak'),ylab='spearman rho Fst ~ dxy',main='Allopatry')
#boxplot(para.fst.dxy.back,para.fst.dxy.peak,outline=F,names=c('back','peak'),ylab='spearman rho Fst ~ dxy',main='Parapatry/Partial Sympatry')
boxplot(para.dxy.pi.back,para.dxy.pi.peak,outline=F,names=c('back','peak'),ylab='spearman rho dxy ~ pi',main='Parapatry/Partial Sympatry')
boxplot(allo.dxy.pi.back,allo.dxy.pi.peak,outline=F,names=c('back','peak'),ylab='spearman rho dxy ~ pi',main='Allopatry')

#wilcox.test(allo.fst.dxy.back,allo.fst.dxy.peak)
#wilcox.test(para.fst.dxy.back,para.fst.dxy.peak)

wilcox.test(para.dxy.pi.back,para.dxy.pi.peak)
wilcox.test(allo.dxy.pi.back,allo.dxy.pi.peak)

wilcox.test(allo.fst.dxy.peak,para.fst.dxy.peak)


