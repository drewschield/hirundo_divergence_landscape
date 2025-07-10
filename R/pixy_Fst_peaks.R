# Hirundo genomic divergence landscape analysis
# pixy_Fst_peaks.R - summarize distributions of dxy and pi in differentiation islands

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)
library(Rmisc)
library(rstatix)
library(ggpubr)
library(car)
library(multcomp)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----
fst <- read.table('./results/pixy.all.1mb.fst.txt',header=T)
dxy <- read.table('./results/pixy.all.1mb.dxy.txt',header=T)
pi <- read.table('./results/pixy.all.1mb.pi.txt',header=T)

## Parse autosomes (we'll focus on autosomal regions and may also do focused Z chromosome stuff too)
fst <- fst %>% filter(chromosome != 'NC_053488.1')
dxy <- dxy %>% filter(chromosome != 'NC_053488.1')
pi <- pi %>% filter(chromosome != 'NC_053488.1')

### Save landscapes with centromere regions parsed out----
fst.nc <- fst %>% filter(!(chromosome == 'NC_053451.1' & window_pos_1 >= 6500000 & window_pos_2 < 12000000) &
                           !(chromosome == 'NC_053453.1' & window_pos_1 >= 39000000 & window_pos_2 < 44000000) &
                           !(chromosome == 'NC_053450.1' & window_pos_1 >= 97500000 & window_pos_2 < 101250000) &
                           !(chromosome == 'NC_053452.1' & window_pos_1 >= 12500000 & window_pos_2 < 17000000) &
                           !(chromosome == 'NC_053454.1' & window_pos_1 >= 65000000 & window_pos_2 < 70000000) &
                           !(chromosome == 'NC_053470.1' & window_pos_1 >= 8400000 & window_pos_2 < 9617204) &
                           !(chromosome == 'NC_053455.1' & window_pos_1 >= 48000000 & window_pos_2 < 55000000) &
                           !(chromosome == 'NC_053457.1' & window_pos_1 >= 0 & window_pos_2 < 6000000) &
                           !(chromosome == 'NC_053456.1' & window_pos_1 >= 15000000 & window_pos_2 < 22000000) &
                           !(chromosome == 'NC_053458.1' & window_pos_1 >= 22500000 & window_pos_2 < 27500000))

dxy.nc <- dxy %>% filter(!(chromosome == 'NC_053451.1' & window_pos_1 >= 6500000 & window_pos_2 < 12000000) &
                           !(chromosome == 'NC_053453.1' & window_pos_1 >= 39000000 & window_pos_2 < 44000000) &
                           !(chromosome == 'NC_053450.1' & window_pos_1 >= 97500000 & window_pos_2 < 101250000) &
                           !(chromosome == 'NC_053452.1' & window_pos_1 >= 12500000 & window_pos_2 < 17000000) &
                           !(chromosome == 'NC_053454.1' & window_pos_1 >= 65000000 & window_pos_2 < 70000000) &
                           !(chromosome == 'NC_053470.1' & window_pos_1 >= 8400000 & window_pos_2 < 9617204) &
                           !(chromosome == 'NC_053455.1' & window_pos_1 >= 48000000 & window_pos_2 < 55000000) &
                           !(chromosome == 'NC_053457.1' & window_pos_1 >= 0 & window_pos_2 < 6000000) &
                           !(chromosome == 'NC_053456.1' & window_pos_1 >= 15000000 & window_pos_2 < 22000000) &
                           !(chromosome == 'NC_053458.1' & window_pos_1 >= 22500000 & window_pos_2 < 27500000))

pi.nc <- pi %>% filter(!(chromosome == 'NC_053451.1' & window_pos_1 >= 6500000 & window_pos_2 < 12000000) &
                           !(chromosome == 'NC_053453.1' & window_pos_1 >= 39000000 & window_pos_2 < 44000000) &
                           !(chromosome == 'NC_053450.1' & window_pos_1 >= 97500000 & window_pos_2 < 101250000) &
                           !(chromosome == 'NC_053452.1' & window_pos_1 >= 12500000 & window_pos_2 < 17000000) &
                           !(chromosome == 'NC_053454.1' & window_pos_1 >= 65000000 & window_pos_2 < 70000000) &
                           !(chromosome == 'NC_053470.1' & window_pos_1 >= 8400000 & window_pos_2 < 9617204) &
                           !(chromosome == 'NC_053455.1' & window_pos_1 >= 48000000 & window_pos_2 < 55000000) &
                           !(chromosome == 'NC_053457.1' & window_pos_1 >= 0 & window_pos_2 < 6000000) &
                           !(chromosome == 'NC_053456.1' & window_pos_1 >= 15000000 & window_pos_2 < 22000000) &
                           !(chromosome == 'NC_053458.1' & window_pos_1 >= 22500000 & window_pos_2 < 27500000))

### Parse species----

## All windows
fst.rus.aet <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.neo.jav <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.smi.dim <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
fst.rus.jav <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
fst.aet.smi <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HS'))

dxy.rus.aet <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
dxy.neo.jav <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
dxy.smi.dim <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
dxy.rus.jav <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
dxy.aet.smi <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HS'))

pi.rus <- pi %>% filter(pop == 'HRR')
pi.aet <- pi %>% filter(pop == 'HA')
pi.smi <- pi %>% filter(pop == 'HS')
pi.neo <- pi %>% filter(pop == 'HN')
pi.jav <- pi %>% filter(pop == 'HT')
pi.dim <- pi %>% filter(pop == 'HD')

## Non-centromere regions
fst.nc.rus.aet <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.nc.neo.jav <- fst.nc %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.nc.smi.dim <- fst.nc %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
fst.nc.rus.jav <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
fst.nc.aet.smi <- fst.nc %>% filter((pop1 == 'HA') & (pop2 == 'HS'))

dxy.nc.rus.aet <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
dxy.nc.neo.jav <- dxy.nc %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
dxy.nc.smi.dim <- dxy.nc %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
dxy.nc.rus.jav <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
dxy.nc.aet.smi <- dxy.nc %>% filter((pop1 == 'HA') & (pop2 == 'HS'))

pi.nc.rus <- pi.nc %>% filter(pop == 'HRR')
pi.nc.aet <- pi.nc %>% filter(pop == 'HA')
pi.nc.smi <- pi.nc %>% filter(pop == 'HS')
pi.nc.neo <- pi.nc %>% filter(pop == 'HN')
pi.nc.jav <- pi.nc %>% filter(pop == 'HT')
pi.nc.dim <- pi.nc %>% filter(pop == 'HD')

### Define and extract statistics in Fst islands and valleys----

## Read in null distributions based on permutations
perm.rus.aet <- read.table('./permutations/permutation.fst.1Mb.rus.aet.txt',header=T)
perm.neo.jav <- read.table('./permutations/permutation.fst.1Mb.neo.jav.txt',header=T)
perm.smi.dim <- read.table('./permutations/permutation.fst.1Mb.smi.dim.txt',header=T)
perm.rus.jav <- read.table('./permutations/permutation.fst.1Mb.rus.jav.txt',header=T)
perm.aet.smi <- read.table('./permutations/permutation.fst.1Mb.aet.smi.txt',header=T)

## Parse Fst for all regions >= null distribution & matching regions in dxy and pi
fst.rus.aet.is <- fst.rus.aet %>% filter(avg_wc_fst >= max(perm.rus.aet$Fst))
fst.neo.jav.is <- fst.neo.jav %>% filter(avg_wc_fst >= max(perm.neo.jav$Fst))
fst.smi.dim.is <- fst.smi.dim %>% filter(avg_wc_fst >= max(perm.smi.dim$Fst))
fst.rus.jav.is <- fst.rus.jav %>% filter(avg_wc_fst >= max(perm.rus.jav$Fst))
fst.aet.smi.is <- fst.aet.smi %>% filter(avg_wc_fst >= max(perm.aet.smi$Fst))

dxy.rus.aet.is <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
dxy.neo.jav.is <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
dxy.smi.dim.is <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
dxy.rus.jav.is <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
dxy.aet.smi.is <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]

pi.rus.is <- pi.rus[which(fst.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
pi.aet.is <- pi.aet[which(fst.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
pi.neo.is <- pi.neo[which(fst.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
pi.jav.is <- pi.jav[which(fst.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
pi.smi.is <- pi.smi[which(fst.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
pi.dim.is <- pi.dim[which(fst.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
pi.rus.rj.is <- pi.rus[which(fst.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
pi.jav.rj.is <- pi.jav[which(fst.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
pi.aet.as.is <- pi.aet[which(fst.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]
pi.smi.as.is <- pi.smi[which(fst.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]

## Parse non-centromere regions above thresholds
fst.nc.rus.aet.is <- fst.nc.rus.aet %>% filter(avg_wc_fst >= max(perm.rus.aet$Fst))
fst.nc.neo.jav.is <- fst.nc.neo.jav %>% filter(avg_wc_fst >= max(perm.neo.jav$Fst))
fst.nc.smi.dim.is <- fst.nc.smi.dim %>% filter(avg_wc_fst >= max(perm.smi.dim$Fst))
fst.nc.rus.jav.is <- fst.nc.rus.jav %>% filter(avg_wc_fst >= max(perm.rus.jav$Fst))
fst.nc.aet.smi.is <- fst.nc.aet.smi %>% filter(avg_wc_fst >= max(perm.aet.smi$Fst))

dxy.nc.rus.aet.is <- dxy.nc.rus.aet[which(fst.nc.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
dxy.nc.neo.jav.is <- dxy.nc.neo.jav[which(fst.nc.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
dxy.nc.smi.dim.is <- dxy.nc.smi.dim[which(fst.nc.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
dxy.nc.rus.jav.is <- dxy.nc.rus.jav[which(fst.nc.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
dxy.nc.aet.smi.is <- dxy.nc.aet.smi[which(fst.nc.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]

pi.nc.rus.is <- pi.nc.rus[which(fst.nc.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
pi.nc.aet.is <- pi.nc.aet[which(fst.nc.rus.aet$avg_wc_fst>=max(perm.rus.aet$Fst)),]
pi.nc.neo.is <- pi.nc.neo[which(fst.nc.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
pi.nc.jav.is <- pi.nc.jav[which(fst.nc.neo.jav$avg_wc_fst>=max(perm.neo.jav$Fst)),]
pi.nc.smi.is <- pi.nc.smi[which(fst.nc.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
pi.nc.dim.is <- pi.nc.dim[which(fst.nc.smi.dim$avg_wc_fst>=max(perm.smi.dim$Fst)),]
pi.nc.rus.rj.is <- pi.nc.rus[which(fst.nc.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
pi.nc.jav.rj.is <- pi.nc.jav[which(fst.nc.rus.jav$avg_wc_fst>=max(perm.rus.jav$Fst)),]
pi.nc.aet.as.is <- pi.nc.aet[which(fst.nc.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]
pi.nc.smi.as.is <- pi.nc.smi[which(fst.nc.aet.smi$avg_wc_fst>=max(perm.aet.smi$Fst)),]

## Extract regions outside islands
dxy.rus.aet.va <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst<max(perm.rus.aet$Fst)),]
dxy.neo.jav.va <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst<max(perm.neo.jav$Fst)),]
dxy.smi.dim.va <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst<max(perm.smi.dim$Fst)),]
dxy.rus.jav.va <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst<max(perm.rus.jav$Fst)),]
dxy.aet.smi.va <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst<max(perm.aet.smi$Fst)),]

pi.rus.va <- pi.rus[which(fst.rus.aet$avg_wc_fst<max(perm.rus.aet$Fst)),]
pi.aet.va <- pi.aet[which(fst.rus.aet$avg_wc_fst<max(perm.rus.aet$Fst)),]
pi.neo.va <- pi.neo[which(fst.neo.jav$avg_wc_fst<max(perm.neo.jav$Fst)),]
pi.jav.va <- pi.jav[which(fst.neo.jav$avg_wc_fst<max(perm.neo.jav$Fst)),]
pi.smi.va <- pi.smi[which(fst.smi.dim$avg_wc_fst<max(perm.smi.dim$Fst)),]
pi.dim.va <- pi.dim[which(fst.smi.dim$avg_wc_fst<max(perm.smi.dim$Fst)),]
pi.rus.rj.va <- pi.rus[which(fst.rus.jav$avg_wc_fst<max(perm.rus.jav$Fst)),]
pi.jav.rj.va <- pi.jav[which(fst.rus.jav$avg_wc_fst<max(perm.rus.jav$Fst)),]
pi.aet.as.va <- pi.aet[which(fst.aet.smi$avg_wc_fst<max(perm.aet.smi$Fst)),]
pi.smi.as.va <- pi.smi[which(fst.aet.smi$avg_wc_fst<max(perm.aet.smi$Fst)),]

### Plot boxplots to compare distributions (for main figure)----

## dxy
par(mfrow=c(2,3))
boxplot(dxy.rus.aet.va$avg_dxy,dxy.rus.aet.is$avg_dxy,dxy.nc.rus.aet.is$avg_dxy,outline=F,ylab='dxy',names=c('Back','Island','Is. No Cent.'),main='rustica-aethiopica')
boxplot(dxy.neo.jav.va$avg_dxy,dxy.neo.jav.is$avg_dxy,dxy.nc.neo.jav.is$avg_dxy,outline=F,ylab='dxy',names=c('Back','Island','Is. No Cent.'),main='neoxena-javanica')
boxplot(dxy.smi.dim.va$avg_dxy,dxy.smi.dim.is$avg_dxy,dxy.nc.smi.dim.is$avg_dxy,outline=F,ylab='dxy',names=c('Back','Island','Is. No Cent.'),main='smithii-dimidiata')

pi.rus.aet.va <- (pi.rus.va$avg_pi + pi.aet.va$avg_pi)/2
pi.neo.jav.va <- (pi.neo.va$avg_pi + pi.jav.va$avg_pi)/2
pi.smi.dim.va <- (pi.smi.va$avg_pi + pi.dim.va$avg_pi)/2
pi.rus.jav.va <- (pi.rus.rj.va$avg_pi + pi.jav.rj.va$avg_pi)/2
pi.aet.smi.va <- (pi.aet.as.va$avg_pi + pi.smi.as.va$avg_pi)/2
pi.rus.aet.is <- (pi.rus.is$avg_pi + pi.aet.is$avg_pi)/2
pi.neo.jav.is <- (pi.neo.is$avg_pi + pi.jav.is$avg_pi)/2
pi.smi.dim.is <- (pi.smi.is$avg_pi + pi.dim.is$avg_pi)/2
pi.rus.jav.is <- (pi.rus.rj.is$avg_pi + pi.jav.rj.is$avg_pi)/2
pi.aet.smi.is <- (pi.aet.as.is$avg_pi + pi.smi.as.is$avg_pi)/2
pi.nc.rus.aet.is <- (pi.nc.rus.is$avg_pi + pi.nc.aet.is$avg_pi)/2
pi.nc.neo.jav.is <- (pi.nc.neo.is$avg_pi + pi.nc.jav.is$avg_pi)/2
pi.nc.smi.dim.is <- (pi.nc.smi.is$avg_pi + pi.nc.dim.is$avg_pi)/2
pi.nc.rus.jav.is <- (pi.nc.rus.rj.is$avg_pi + pi.nc.jav.rj.is$avg_pi)/2
pi.nc.aet.smi.is <- (pi.nc.aet.as.is$avg_pi + pi.nc.smi.as.is$avg_pi)/2

boxplot(pi.rus.aet.va,pi.rus.aet.is,pi.nc.rus.aet.is,outline=F,ylab='Mean pi',names=c('Back','Island','Is. No Cent.'))
boxplot(pi.neo.jav.va,pi.neo.jav.is,pi.nc.neo.jav.is,outline=F,ylab='Mean pi',names=c('Back','Island','Is. No Cent.'))
boxplot(pi.smi.dim.va,pi.smi.dim.is,pi.nc.smi.dim.is,outline=F,ylab='Mean pi',names=c('Back','Island','Is. No Cent.'))

## Output at 5.5 x 6, then make second plot for additional comparisons

boxplot(dxy.rus.jav.va$avg_dxy,dxy.rus.jav.is$avg_dxy,dxy.nc.rus.jav.is$avg_dxy,outline=F,ylab='dxy',names=c('Back','Island','Is. No Cent.'),main='rustica-javanica')
boxplot(dxy.aet.smi.va$avg_dxy,dxy.aet.smi.is$avg_dxy,dxy.nc.aet.smi.is$avg_dxy,outline=F,ylab='dxy',names=c('Back','Island','Is. No Cent.'),main='aethiopica-smithii')
plot(1)
boxplot(pi.rus.jav.va,pi.rus.jav.is,pi.nc.rus.jav.is,outline=F,ylab='Mean pi',names=c('Back','Island','Is. No Cent.'))
boxplot(pi.aet.smi.va,pi.aet.smi.is,pi.nc.aet.smi.is,outline=F,ylab='Mean pi',names=c('Back','Island','Is. No Cent.'))
plot(1)

### Perform ANOVA and Tukey post-hoc tests on distributions----

## Set up data frames
x = c(rep("va",length(dxy.rus.aet.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.aet.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.aet.is$avg_dxy)))
dxy.rus.aet <- data.frame(
  dxy = c(dxy.rus.aet.va$avg_dxy,dxy.rus.aet.is$avg_dxy,dxy.nc.rus.aet.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.neo.jav.va$avg_dxy)))
y = c(rep("is",length(dxy.neo.jav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.neo.jav.is$avg_dxy)))
dxy.neo.jav <- data.frame(
  dxy = c(dxy.neo.jav.va$avg_dxy,dxy.neo.jav.is$avg_dxy,dxy.nc.neo.jav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.smi.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.smi.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.smi.dim.is$avg_dxy)))
dxy.smi.dim <- data.frame(
  dxy = c(dxy.smi.dim.va$avg_dxy,dxy.smi.dim.is$avg_dxy,dxy.nc.smi.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.aet.va)))
y = c(rep("is",length(pi.rus.aet.is)))
z = c(rep("nc",length(pi.nc.rus.aet.is)))
pi.rus.aet <- data.frame(
  pi = c(pi.rus.aet.va,pi.rus.aet.is,pi.nc.rus.aet.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.neo.jav.va)))
y = c(rep("is",length(pi.neo.jav.is)))
z = c(rep("nc",length(pi.nc.neo.jav.is)))
pi.neo.jav <- data.frame(
  pi = c(pi.neo.jav.va,pi.neo.jav.is,pi.nc.neo.jav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.smi.dim.va)))
y = c(rep("is",length(pi.smi.dim.is)))
z = c(rep("nc",length(pi.nc.smi.dim.is)))
pi.smi.dim <- data.frame(
  pi = c(pi.smi.dim.va,pi.smi.dim.is,pi.nc.smi.dim.is),
  loci = factor(c(x,y,z))
)

## Perform ANOVA
aov.dxy.rus.aet <- aov(dxy ~ loci, data = dxy.rus.aet)
summary(aov.dxy.rus.aet)

aov.dxy.neo.jav <- aov(dxy ~ loci, data = dxy.neo.jav)
summary(aov.dxy.neo.jav)

aov.dxy.smi.dim <- aov(dxy ~ loci, data = dxy.smi.dim)
summary(aov.dxy.smi.dim)

aov.pi.rus.aet <- aov(pi ~ loci, data = pi.rus.aet)
summary(aov.pi.rus.aet)

aov.pi.neo.jav <- aov(pi ~ loci, data = pi.neo.jav)
summary(aov.pi.neo.jav)

aov.pi.smi.dim <- aov(pi ~ loci, data = pi.smi.dim)
summary(aov.pi.smi.dim)

## These are consistently significant, warranting post-hoc tests correcting for mulitple tests

## Tukey's post-hoc tests correcting for multiple testing
post.dxy.rus.aet <- glht(aov.dxy.rus.aet, linfct = mcp(loci = "Tukey"))
summary(post.dxy.rus.aet)

post.dxy.neo.jav <- glht(aov.dxy.neo.jav, linfct = mcp(loci = "Tukey"))
summary(post.dxy.neo.jav)

post.dxy.smi.dim <- glht(aov.dxy.smi.dim, linfct = mcp(loci = "Tukey"))
summary(post.dxy.smi.dim)

post.pi.rus.aet <- glht(aov.pi.rus.aet, linfct = mcp(loci = "Tukey"))
summary(post.pi.rus.aet)

post.pi.neo.jav <- glht(aov.pi.neo.jav, linfct = mcp(loci = "Tukey"))
summary(post.pi.neo.jav)

post.pi.smi.dim <- glht(aov.pi.smi.dim, linfct = mcp(loci = "Tukey"))
summary(post.pi.smi.dim)

## These show consistent significant differences between islands (with or without centromeres) and the genome background.
## There is no significant differences between islands with/without centromeres.

