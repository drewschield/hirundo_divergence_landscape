# Hirundo genomic divergence landscape analysis
# pixy_Fst_peaks_all.R - summarize distributions of dxy and pi in differentiation islands for all pairs of species with n > 1

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

### Parse species pairs----

## All windows

fst.rus.sav <- fst %>% filter((pop1 == 'HRS') & (pop2 == 'HRR'))
fst.rus.tyt <- fst %>% filter((pop1 == 'HRT') & (pop2 == 'HRR'))
fst.rus.aet <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.rus.smi <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
fst.rus.neo <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
fst.rus.jav <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
fst.rus.dim <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
fst.aet.smi <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
fst.aet.neo <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
fst.aet.jav <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
fst.aet.dim <- fst %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
fst.smi.neo <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
fst.smi.jav <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
fst.smi.dim <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
fst.neo.jav <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.neo.dim <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
fst.jav.dim <- fst %>% filter((pop1 == 'HT') & (pop2 == 'HD'))

dxy.rus.sav <- dxy %>% filter((pop1 == 'HRS') & (pop2 == 'HRR'))
dxy.rus.tyt <- dxy %>% filter((pop1 == 'HRT') & (pop2 == 'HRR'))
dxy.rus.aet <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
dxy.rus.smi <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
dxy.rus.neo <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
dxy.rus.jav <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
dxy.rus.dim <- dxy %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
dxy.aet.smi <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
dxy.aet.neo <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
dxy.aet.jav <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
dxy.aet.dim <- dxy %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
dxy.smi.neo <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
dxy.smi.jav <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
dxy.smi.dim <- dxy %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
dxy.neo.jav <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
dxy.neo.dim <- dxy %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
dxy.jav.dim <- dxy %>% filter((pop1 == 'HT') & (pop2 == 'HD'))

pi.rus <- pi %>% filter(pop == 'HRR')
pi.sav <- pi %>% filter(pop == 'HRS')
pi.tyt <- pi %>% filter(pop == 'HRT')
pi.aet <- pi %>% filter(pop == 'HA')
pi.smi <- pi %>% filter(pop == 'HS')
pi.neo <- pi %>% filter(pop == 'HN')
pi.jav <- pi %>% filter(pop == 'HT')
pi.dim <- pi %>% filter(pop == 'HD')

## Outside centromeres
fst.nc.rus.sav <- fst.nc %>% filter((pop1 == 'HRS') & (pop2 == 'HRR'))
fst.nc.rus.tyt <- fst.nc %>% filter((pop1 == 'HRT') & (pop2 == 'HRR'))
fst.nc.rus.aet <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.nc.rus.smi <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
fst.nc.rus.neo <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
fst.nc.rus.jav <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
fst.nc.rus.dim <- fst.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
fst.nc.aet.smi <- fst.nc %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
fst.nc.aet.neo <- fst.nc %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
fst.nc.aet.jav <- fst.nc %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
fst.nc.aet.dim <- fst.nc %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
fst.nc.smi.neo <- fst.nc %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
fst.nc.smi.jav <- fst.nc %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
fst.nc.smi.dim <- fst.nc %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
fst.nc.neo.jav <- fst.nc %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.nc.neo.dim <- fst.nc %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
fst.nc.jav.dim <- fst.nc %>% filter((pop1 == 'HT') & (pop2 == 'HD'))

dxy.nc.rus.sav <- dxy.nc %>% filter((pop1 == 'HRS') & (pop2 == 'HRR'))
dxy.nc.rus.tyt <- dxy.nc %>% filter((pop1 == 'HRT') & (pop2 == 'HRR'))
dxy.nc.rus.aet <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
dxy.nc.rus.smi <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HS'))
dxy.nc.rus.neo <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HN'))
dxy.nc.rus.jav <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HT'))
dxy.nc.rus.dim <- dxy.nc %>% filter((pop1 == 'HRR') & (pop2 == 'HD'))
dxy.nc.aet.smi <- dxy.nc %>% filter((pop1 == 'HA') & (pop2 == 'HS'))
dxy.nc.aet.neo <- dxy.nc %>% filter((pop1 == 'HA') & (pop2 == 'HN'))
dxy.nc.aet.jav <- dxy.nc %>% filter((pop1 == 'HA') & (pop2 == 'HT'))
dxy.nc.aet.dim <- dxy.nc %>% filter((pop1 == 'HA') & (pop2 == 'HD'))
dxy.nc.smi.neo <- dxy.nc %>% filter((pop1 == 'HS') & (pop2 == 'HN'))
dxy.nc.smi.jav <- dxy.nc %>% filter((pop1 == 'HS') & (pop2 == 'HT'))
dxy.nc.smi.dim <- dxy.nc %>% filter((pop1 == 'HS') & (pop2 == 'HD'))
dxy.nc.neo.jav <- dxy.nc %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
dxy.nc.neo.dim <- dxy.nc %>% filter((pop1 == 'HN') & (pop2 == 'HD'))
dxy.nc.jav.dim <- dxy.nc %>% filter((pop1 == 'HT') & (pop2 == 'HD'))

pi.nc.rus <- pi.nc %>% filter(pop == 'HRR')
pi.nc.sav <- pi.nc %>% filter(pop == 'HRS')
pi.nc.tyt <- pi.nc %>% filter(pop == 'HRT')
pi.nc.aet <- pi.nc %>% filter(pop == 'HA')
pi.nc.smi <- pi.nc %>% filter(pop == 'HS')
pi.nc.neo <- pi.nc %>% filter(pop == 'HN')
pi.nc.jav <- pi.nc %>% filter(pop == 'HT')
pi.nc.dim <- pi.nc %>% filter(pop == 'HD')

### Extract summary statistics in all Fst peaks----

## Read in null distributions based on permutations
perm.rus.sav <- read.table('./permutations/permutation.fst.1Mb.rus.sav.txt',header=T)
perm.rus.tyt <- read.table('./permutations/permutation.fst.1Mb.rus.tyt.txt',header=T)
perm.rus.aet <- read.table('./permutations/permutation.fst.1Mb.rus.aet.txt',header=T)
perm.rus.smi <- read.table('./permutations/permutation.fst.1Mb.rus.smi.txt',header=T)
perm.rus.neo <- read.table('./permutations/permutation.fst.1Mb.rus.neo.txt',header=T)
perm.rus.jav <- read.table('./permutations/permutation.fst.1Mb.rus.jav.txt',header=T)
perm.rus.dim <- read.table('./permutations/permutation.fst.1Mb.rus.dim.txt',header=T)
perm.aet.smi <- read.table('./permutations/permutation.fst.1Mb.aet.smi.txt',header=T)
perm.aet.neo <- read.table('./permutations/permutation.fst.1Mb.aet.neo.txt',header=T)
perm.aet.jav <- read.table('./permutations/permutation.fst.1Mb.aet.jav.txt',header=T)
perm.aet.dim <- read.table('./permutations/permutation.fst.1Mb.aet.dim.txt',header=T)
perm.smi.neo <- read.table('./permutations/permutation.fst.1Mb.smi.neo.txt',header=T)
perm.smi.jav <- read.table('./permutations/permutation.fst.1Mb.smi.jav.txt',header=T)
perm.smi.dim <- read.table('./permutations/permutation.fst.1Mb.smi.dim.txt',header=T)
perm.neo.jav <- read.table('./permutations/permutation.fst.1Mb.neo.jav.txt',header=T)
perm.neo.dim <- read.table('./permutations/permutation.fst.1Mb.neo.dim.txt',header=T)
perm.jav.dim <- read.table('./permutations/permutation.fst.1Mb.jav.dim.txt',header=T)

## Define and extract Fst peaks
q.rus.sav <- max(perm.rus.sav$Fst)
q.rus.tyt <- max(perm.rus.tyt$Fst)
q.rus.aet <- max(perm.rus.aet$Fst)
q.rus.smi <- max(perm.rus.smi$Fst)
q.rus.neo <- max(perm.rus.neo$Fst)
q.rus.jav <- max(perm.rus.jav$Fst)
q.rus.dim <- max(perm.rus.dim$Fst)
q.aet.smi <- max(perm.aet.smi$Fst)
q.aet.neo <- max(perm.aet.neo$Fst)
q.aet.jav <- max(perm.aet.jav$Fst)
q.aet.dim <- max(perm.aet.dim$Fst)
q.smi.neo <- max(perm.smi.neo$Fst)
q.smi.jav <- max(perm.smi.jav$Fst)
q.smi.dim <- max(perm.smi.dim$Fst)
q.neo.jav <- max(perm.neo.jav$Fst)
q.neo.dim <- max(perm.neo.dim$Fst)
q.jav.dim <- max(perm.jav.dim$Fst)

#q.rus.sav <- quantile(fst.rus.sav$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.tyt <- quantile(fst.rus.tyt$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.aet <- quantile(fst.rus.aet$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.smi <- quantile(fst.rus.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.neo <- quantile(fst.rus.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.jav <- quantile(fst.rus.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.rus.dim <- quantile(fst.rus.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.aet.smi <- quantile(fst.aet.smi$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.aet.neo <- quantile(fst.aet.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.aet.jav <- quantile(fst.aet.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.aet.dim <- quantile(fst.aet.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.smi.neo <- quantile(fst.smi.neo$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.smi.jav <- quantile(fst.smi.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.smi.dim <- quantile(fst.smi.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.neo.jav <- quantile(fst.neo.jav$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.neo.dim <- quantile(fst.neo.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)
#q.jav.dim <- quantile(fst.jav.dim$avg_wc_fst,c(0.90),na.rm=T,names=F)

fst.rus.sav.is <- fst.rus.sav %>% filter(avg_wc_fst >= q.rus.sav)
fst.rus.tyt.is <- fst.rus.tyt %>% filter(avg_wc_fst >= q.rus.tyt)
fst.rus.aet.is <- fst.rus.aet %>% filter(avg_wc_fst >= q.rus.aet)
fst.rus.smi.is <- fst.rus.smi %>% filter(avg_wc_fst >= q.rus.smi)
fst.rus.neo.is <- fst.rus.neo %>% filter(avg_wc_fst >= q.rus.neo)
fst.rus.jav.is <- fst.rus.jav %>% filter(avg_wc_fst >= q.rus.jav)
fst.rus.dim.is <- fst.rus.dim %>% filter(avg_wc_fst >= q.rus.dim)
fst.aet.smi.is <- fst.aet.smi %>% filter(avg_wc_fst >= q.aet.smi)
fst.aet.neo.is <- fst.aet.neo %>% filter(avg_wc_fst >= q.aet.neo)
fst.aet.jav.is <- fst.aet.jav %>% filter(avg_wc_fst >= q.aet.jav)
fst.aet.dim.is <- fst.aet.dim %>% filter(avg_wc_fst >= q.aet.dim)
fst.smi.neo.is <- fst.smi.neo %>% filter(avg_wc_fst >= q.smi.neo)
fst.smi.jav.is <- fst.smi.jav %>% filter(avg_wc_fst >= q.smi.jav)
fst.smi.dim.is <- fst.smi.dim %>% filter(avg_wc_fst >= q.smi.dim)
fst.neo.jav.is <- fst.neo.jav %>% filter(avg_wc_fst >= q.neo.jav)
fst.neo.dim.is <- fst.neo.dim %>% filter(avg_wc_fst >= q.neo.dim)
fst.jav.dim.is <- fst.jav.dim %>% filter(avg_wc_fst >= q.jav.dim)

## Extract dxy and pi in peaks
dxy.rus.sav.is <- dxy.rus.sav[which(fst.rus.sav$avg_wc_fst>=q.rus.sav),]
dxy.rus.tyt.is <- dxy.rus.tyt[which(fst.rus.tyt$avg_wc_fst>=q.rus.tyt),]
dxy.rus.aet.is <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
dxy.rus.smi.is <- dxy.rus.smi[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
dxy.rus.neo.is <- dxy.rus.neo[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
dxy.rus.jav.is <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
dxy.rus.dim.is <- dxy.rus.dim[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
dxy.aet.smi.is <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
dxy.aet.neo.is <- dxy.aet.neo[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
dxy.aet.jav.is <- dxy.aet.jav[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
dxy.aet.dim.is <- dxy.aet.dim[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
dxy.smi.neo.is <- dxy.smi.neo[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
dxy.smi.jav.is <- dxy.smi.jav[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
dxy.smi.dim.is <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
dxy.neo.jav.is <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
dxy.neo.dim.is <- dxy.neo.dim[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
dxy.jav.dim.is <- dxy.jav.dim[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]

pi.rus.rus.sav.is <- pi.rus[which(fst.rus.sav$avg_wc_fst>=q.rus.sav),]
pi.rus.rus.tyt.is <- pi.rus[which(fst.rus.tyt$avg_wc_fst>=q.rus.tyt),]
pi.rus.rus.aet.is <- pi.rus[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.rus.rus.smi.is <- pi.rus[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.rus.rus.neo.is <- pi.rus[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.rus.rus.jav.is <- pi.rus[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.rus.rus.dim.is <- pi.rus[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.aet.aet.smi.is <- pi.aet[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.aet.aet.neo.is <- pi.aet[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.aet.aet.jav.is <- pi.aet[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.aet.aet.dim.is <- pi.aet[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.smi.smi.neo.is <- pi.smi[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.smi.smi.jav.is <- pi.smi[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.smi.smi.dim.is <- pi.smi[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.neo.neo.jav.is <- pi.neo[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.neo.neo.dim.is <- pi.neo[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.jav.jav.dim.is <- pi.jav[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]

pi.sav.rus.sav.is <- pi.sav[which(fst.rus.sav$avg_wc_fst>=q.rus.sav),]
pi.tyt.rus.tyt.is <- pi.tyt[which(fst.rus.tyt$avg_wc_fst>=q.rus.tyt),]
pi.aet.rus.aet.is <- pi.aet[which(fst.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.smi.rus.smi.is <- pi.smi[which(fst.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.neo.rus.neo.is <- pi.neo[which(fst.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.jav.rus.jav.is <- pi.jav[which(fst.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.dim.rus.dim.is <- pi.dim[which(fst.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.smi.aet.smi.is <- pi.smi[which(fst.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.neo.aet.neo.is <- pi.neo[which(fst.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.jav.aet.jav.is <- pi.jav[which(fst.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.dim.aet.dim.is <- pi.dim[which(fst.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.neo.smi.neo.is <- pi.neo[which(fst.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.jav.smi.jav.is <- pi.jav[which(fst.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.dim.smi.dim.is <- pi.dim[which(fst.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.jav.neo.jav.is <- pi.jav[which(fst.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.dim.neo.dim.is <- pi.dim[which(fst.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.dim.jav.dim.is <- pi.dim[which(fst.jav.dim$avg_wc_fst>=q.jav.dim),]

### Extract summary statistics in Fst peaks outside centromeres----

## Extract Fst peaks
fst.nc.rus.sav.is <- fst.nc.rus.sav %>% filter(avg_wc_fst >= q.rus.sav)
fst.nc.rus.tyt.is <- fst.nc.rus.tyt %>% filter(avg_wc_fst >= q.rus.tyt)
fst.nc.rus.aet.is <- fst.nc.rus.aet %>% filter(avg_wc_fst >= q.rus.aet)
fst.nc.rus.smi.is <- fst.nc.rus.smi %>% filter(avg_wc_fst >= q.rus.smi)
fst.nc.rus.neo.is <- fst.nc.rus.neo %>% filter(avg_wc_fst >= q.rus.neo)
fst.nc.rus.jav.is <- fst.nc.rus.jav %>% filter(avg_wc_fst >= q.rus.jav)
fst.nc.rus.dim.is <- fst.nc.rus.dim %>% filter(avg_wc_fst >= q.rus.dim)
fst.nc.aet.smi.is <- fst.nc.aet.smi %>% filter(avg_wc_fst >= q.aet.smi)
fst.nc.aet.neo.is <- fst.nc.aet.neo %>% filter(avg_wc_fst >= q.aet.neo)
fst.nc.aet.jav.is <- fst.nc.aet.jav %>% filter(avg_wc_fst >= q.aet.jav)
fst.nc.aet.dim.is <- fst.nc.aet.dim %>% filter(avg_wc_fst >= q.aet.dim)
fst.nc.smi.neo.is <- fst.nc.smi.neo %>% filter(avg_wc_fst >= q.smi.neo)
fst.nc.smi.jav.is <- fst.nc.smi.jav %>% filter(avg_wc_fst >= q.smi.jav)
fst.nc.smi.dim.is <- fst.nc.smi.dim %>% filter(avg_wc_fst >= q.smi.dim)
fst.nc.neo.jav.is <- fst.nc.neo.jav %>% filter(avg_wc_fst >= q.neo.jav)
fst.nc.neo.dim.is <- fst.nc.neo.dim %>% filter(avg_wc_fst >= q.neo.dim)
fst.nc.jav.dim.is <- fst.nc.jav.dim %>% filter(avg_wc_fst >= q.jav.dim)

## Extract dxy and pi in peaks
dxy.nc.rus.sav.is <- dxy.nc.rus.sav[which(fst.nc.rus.sav$avg_wc_fst>=q.rus.sav),]
dxy.nc.rus.tyt.is <- dxy.nc.rus.tyt[which(fst.nc.rus.tyt$avg_wc_fst>=q.rus.tyt),]
dxy.nc.rus.aet.is <- dxy.nc.rus.aet[which(fst.nc.rus.aet$avg_wc_fst>=q.rus.aet),]
dxy.nc.rus.smi.is <- dxy.nc.rus.smi[which(fst.nc.rus.smi$avg_wc_fst>=q.rus.smi),]
dxy.nc.rus.neo.is <- dxy.nc.rus.neo[which(fst.nc.rus.neo$avg_wc_fst>=q.rus.neo),]
dxy.nc.rus.jav.is <- dxy.nc.rus.jav[which(fst.nc.rus.jav$avg_wc_fst>=q.rus.jav),]
dxy.nc.rus.dim.is <- dxy.nc.rus.dim[which(fst.nc.rus.dim$avg_wc_fst>=q.rus.dim),]
dxy.nc.aet.smi.is <- dxy.nc.aet.smi[which(fst.nc.aet.smi$avg_wc_fst>=q.aet.smi),]
dxy.nc.aet.neo.is <- dxy.nc.aet.neo[which(fst.nc.aet.neo$avg_wc_fst>=q.aet.neo),]
dxy.nc.aet.jav.is <- dxy.nc.aet.jav[which(fst.nc.aet.jav$avg_wc_fst>=q.aet.jav),]
dxy.nc.aet.dim.is <- dxy.nc.aet.dim[which(fst.nc.aet.dim$avg_wc_fst>=q.aet.dim),]
dxy.nc.smi.neo.is <- dxy.nc.smi.neo[which(fst.nc.smi.neo$avg_wc_fst>=q.smi.neo),]
dxy.nc.smi.jav.is <- dxy.nc.smi.jav[which(fst.nc.smi.jav$avg_wc_fst>=q.smi.jav),]
dxy.nc.smi.dim.is <- dxy.nc.smi.dim[which(fst.nc.smi.dim$avg_wc_fst>=q.smi.dim),]
dxy.nc.neo.jav.is <- dxy.nc.neo.jav[which(fst.nc.neo.jav$avg_wc_fst>=q.neo.jav),]
dxy.nc.neo.dim.is <- dxy.nc.neo.dim[which(fst.nc.neo.dim$avg_wc_fst>=q.neo.dim),]
dxy.nc.jav.dim.is <- dxy.nc.jav.dim[which(fst.nc.jav.dim$avg_wc_fst>=q.jav.dim),]

pi.nc.rus.rus.sav.is <- pi.nc.rus[which(fst.nc.rus.sav$avg_wc_fst>=q.rus.sav),]
pi.nc.rus.rus.tyt.is <- pi.nc.rus[which(fst.nc.rus.tyt$avg_wc_fst>=q.rus.tyt),]
pi.nc.rus.rus.aet.is <- pi.nc.rus[which(fst.nc.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.nc.rus.rus.smi.is <- pi.nc.rus[which(fst.nc.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.nc.rus.rus.neo.is <- pi.nc.rus[which(fst.nc.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.nc.rus.rus.jav.is <- pi.nc.rus[which(fst.nc.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.nc.rus.rus.dim.is <- pi.nc.rus[which(fst.nc.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.nc.aet.aet.smi.is <- pi.nc.aet[which(fst.nc.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.nc.aet.aet.neo.is <- pi.nc.aet[which(fst.nc.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.nc.aet.aet.jav.is <- pi.nc.aet[which(fst.nc.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.nc.aet.aet.dim.is <- pi.nc.aet[which(fst.nc.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.nc.smi.smi.neo.is <- pi.nc.smi[which(fst.nc.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.nc.smi.smi.jav.is <- pi.nc.smi[which(fst.nc.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.nc.smi.smi.dim.is <- pi.nc.smi[which(fst.nc.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.nc.neo.neo.jav.is <- pi.nc.neo[which(fst.nc.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.nc.neo.neo.dim.is <- pi.nc.neo[which(fst.nc.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.nc.jav.jav.dim.is <- pi.nc.jav[which(fst.nc.jav.dim$avg_wc_fst>=q.jav.dim),]

pi.nc.sav.rus.sav.is <- pi.nc.sav[which(fst.nc.rus.sav$avg_wc_fst>=q.rus.sav),]
pi.nc.tyt.rus.tyt.is <- pi.nc.tyt[which(fst.nc.rus.tyt$avg_wc_fst>=q.rus.tyt),]
pi.nc.aet.rus.aet.is <- pi.nc.aet[which(fst.nc.rus.aet$avg_wc_fst>=q.rus.aet),]
pi.nc.smi.rus.smi.is <- pi.nc.smi[which(fst.nc.rus.smi$avg_wc_fst>=q.rus.smi),]
pi.nc.neo.rus.neo.is <- pi.nc.neo[which(fst.nc.rus.neo$avg_wc_fst>=q.rus.neo),]
pi.nc.jav.rus.jav.is <- pi.nc.jav[which(fst.nc.rus.jav$avg_wc_fst>=q.rus.jav),]
pi.nc.dim.rus.dim.is <- pi.nc.dim[which(fst.nc.rus.dim$avg_wc_fst>=q.rus.dim),]
pi.nc.smi.aet.smi.is <- pi.nc.smi[which(fst.nc.aet.smi$avg_wc_fst>=q.aet.smi),]
pi.nc.neo.aet.neo.is <- pi.nc.neo[which(fst.nc.aet.neo$avg_wc_fst>=q.aet.neo),]
pi.nc.jav.aet.jav.is <- pi.nc.jav[which(fst.nc.aet.jav$avg_wc_fst>=q.aet.jav),]
pi.nc.dim.aet.dim.is <- pi.nc.dim[which(fst.nc.aet.dim$avg_wc_fst>=q.aet.dim),]
pi.nc.neo.smi.neo.is <- pi.nc.neo[which(fst.nc.smi.neo$avg_wc_fst>=q.smi.neo),]
pi.nc.jav.smi.jav.is <- pi.nc.jav[which(fst.nc.smi.jav$avg_wc_fst>=q.smi.jav),]
pi.nc.dim.smi.dim.is <- pi.nc.dim[which(fst.nc.smi.dim$avg_wc_fst>=q.smi.dim),]
pi.nc.jav.neo.jav.is <- pi.nc.jav[which(fst.nc.neo.jav$avg_wc_fst>=q.neo.jav),]
pi.nc.dim.neo.dim.is <- pi.nc.dim[which(fst.nc.neo.dim$avg_wc_fst>=q.neo.dim),]
pi.nc.dim.jav.dim.is <- pi.nc.dim[which(fst.nc.jav.dim$avg_wc_fst>=q.jav.dim),]

### Extract summary statistics in genome backgrounds----

fst.rus.sav.va <- fst.rus.sav %>% filter(avg_wc_fst < q.rus.sav)
fst.rus.tyt.va <- fst.rus.tyt %>% filter(avg_wc_fst < q.rus.tyt)
fst.rus.aet.va <- fst.rus.aet %>% filter(avg_wc_fst < q.rus.aet)
fst.rus.smi.va <- fst.rus.smi %>% filter(avg_wc_fst < q.rus.smi)
fst.rus.neo.va <- fst.rus.neo %>% filter(avg_wc_fst < q.rus.neo)
fst.rus.jav.va <- fst.rus.jav %>% filter(avg_wc_fst < q.rus.jav)
fst.rus.dim.va <- fst.rus.dim %>% filter(avg_wc_fst < q.rus.dim)
fst.aet.smi.va <- fst.aet.smi %>% filter(avg_wc_fst < q.aet.smi)
fst.aet.neo.va <- fst.aet.neo %>% filter(avg_wc_fst < q.aet.neo)
fst.aet.jav.va <- fst.aet.jav %>% filter(avg_wc_fst < q.aet.jav)
fst.aet.dim.va <- fst.aet.dim %>% filter(avg_wc_fst < q.aet.dim)
fst.smi.neo.va <- fst.smi.neo %>% filter(avg_wc_fst < q.smi.neo)
fst.smi.jav.va <- fst.smi.jav %>% filter(avg_wc_fst < q.smi.jav)
fst.smi.dim.va <- fst.smi.dim %>% filter(avg_wc_fst < q.smi.dim)
fst.neo.jav.va <- fst.neo.jav %>% filter(avg_wc_fst < q.neo.jav)
fst.neo.dim.va <- fst.neo.dim %>% filter(avg_wc_fst < q.neo.dim)
fst.jav.dim.va <- fst.jav.dim %>% filter(avg_wc_fst < q.jav.dim)

## Extract dxy and pi in peaks

dxy.rus.sav.va <- dxy.rus.sav[which(fst.rus.sav$avg_wc_fst<q.rus.sav),]
dxy.rus.tyt.va <- dxy.rus.tyt[which(fst.rus.tyt$avg_wc_fst<q.rus.tyt),]
dxy.rus.aet.va <- dxy.rus.aet[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
dxy.rus.smi.va <- dxy.rus.smi[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
dxy.rus.neo.va <- dxy.rus.neo[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
dxy.rus.jav.va <- dxy.rus.jav[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
dxy.rus.dim.va <- dxy.rus.dim[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
dxy.aet.smi.va <- dxy.aet.smi[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
dxy.aet.neo.va <- dxy.aet.neo[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
dxy.aet.jav.va <- dxy.aet.jav[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
dxy.aet.dim.va <- dxy.aet.dim[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
dxy.smi.neo.va <- dxy.smi.neo[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
dxy.smi.jav.va <- dxy.smi.jav[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
dxy.smi.dim.va <- dxy.smi.dim[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
dxy.neo.jav.va <- dxy.neo.jav[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
dxy.neo.dim.va <- dxy.neo.dim[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
dxy.jav.dim.va <- dxy.jav.dim[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]

pi.rus.rus.sav.va <- pi.rus[which(fst.rus.sav$avg_wc_fst<q.rus.sav),]
pi.rus.rus.tyt.va <- pi.rus[which(fst.rus.tyt$avg_wc_fst<q.rus.tyt),]
pi.rus.rus.aet.va <- pi.rus[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
pi.rus.rus.smi.va <- pi.rus[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
pi.rus.rus.neo.va <- pi.rus[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
pi.rus.rus.jav.va <- pi.rus[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
pi.rus.rus.dim.va <- pi.rus[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
pi.aet.aet.smi.va <- pi.aet[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
pi.aet.aet.neo.va <- pi.aet[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
pi.aet.aet.jav.va <- pi.aet[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
pi.aet.aet.dim.va <- pi.aet[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
pi.smi.smi.neo.va <- pi.smi[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
pi.smi.smi.jav.va <- pi.smi[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
pi.smi.smi.dim.va <- pi.smi[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
pi.neo.neo.jav.va <- pi.neo[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
pi.neo.neo.dim.va <- pi.neo[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
pi.jav.jav.dim.va <- pi.jav[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]

pi.sav.rus.sav.va <- pi.sav[which(fst.rus.sav$avg_wc_fst<q.rus.sav),]
pi.tyt.rus.tyt.va <- pi.tyt[which(fst.rus.tyt$avg_wc_fst<q.rus.tyt),]
pi.aet.rus.aet.va <- pi.aet[which(fst.rus.aet$avg_wc_fst<q.rus.aet),]
pi.smi.rus.smi.va <- pi.smi[which(fst.rus.smi$avg_wc_fst<q.rus.smi),]
pi.neo.rus.neo.va <- pi.neo[which(fst.rus.neo$avg_wc_fst<q.rus.neo),]
pi.jav.rus.jav.va <- pi.jav[which(fst.rus.jav$avg_wc_fst<q.rus.jav),]
pi.dim.rus.dim.va <- pi.dim[which(fst.rus.dim$avg_wc_fst<q.rus.dim),]
pi.smi.aet.smi.va <- pi.smi[which(fst.aet.smi$avg_wc_fst<q.aet.smi),]
pi.neo.aet.neo.va <- pi.neo[which(fst.aet.neo$avg_wc_fst<q.aet.neo),]
pi.jav.aet.jav.va <- pi.jav[which(fst.aet.jav$avg_wc_fst<q.aet.jav),]
pi.dim.aet.dim.va <- pi.dim[which(fst.aet.dim$avg_wc_fst<q.aet.dim),]
pi.neo.smi.neo.va <- pi.neo[which(fst.smi.neo$avg_wc_fst<q.smi.neo),]
pi.jav.smi.jav.va <- pi.jav[which(fst.smi.jav$avg_wc_fst<q.smi.jav),]
pi.dim.smi.dim.va <- pi.dim[which(fst.smi.dim$avg_wc_fst<q.smi.dim),]
pi.jav.neo.jav.va <- pi.jav[which(fst.neo.jav$avg_wc_fst<q.neo.jav),]
pi.dim.neo.dim.va <- pi.dim[which(fst.neo.dim$avg_wc_fst<q.neo.dim),]
pi.dim.jav.dim.va <- pi.dim[which(fst.jav.dim$avg_wc_fst<q.jav.dim),]

## Test that this worked as expected

par(mfrow=c(1,3))
boxplot(dxy.rus.aet.va$avg_dxy,dxy.rus.aet.is$avg_dxy,dxy.nc.rus.aet.is$avg_dxy)
boxplot(dxy.neo.jav.va$avg_dxy,dxy.neo.jav.is$avg_dxy,dxy.nc.neo.jav.is$avg_dxy)
boxplot(dxy.smi.dim.va$avg_dxy,dxy.smi.dim.is$avg_dxy,dxy.nc.smi.dim.is$avg_dxy)

boxplot(dxy.rus.sav.va$avg_dxy,dxy.rus.sav.is$avg_dxy,dxy.nc.rus.sav.is$avg_dxy)
boxplot(dxy.rus.tyt.va$avg_dxy,dxy.rus.tyt.is$avg_dxy,dxy.nc.rus.tyt.is$avg_dxy)
## Yup.

### Perform ANOVA and Tukey post-hoc tests (dxy)----

## Set up data frames for ANOVA
x = c(rep("va",length(dxy.rus.sav.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.sav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.sav.is$avg_dxy)))
dxy.rus.sav <- data.frame(
  dxy = c(dxy.rus.sav.va$avg_dxy,dxy.rus.sav.is$avg_dxy,dxy.nc.rus.sav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.tyt.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.tyt.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.tyt.is$avg_dxy)))
dxy.rus.tyt <- data.frame(
  dxy = c(dxy.rus.tyt.va$avg_dxy,dxy.rus.tyt.is$avg_dxy,dxy.nc.rus.tyt.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.aet.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.aet.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.aet.is$avg_dxy)))
dxy.rus.aet <- data.frame(
  dxy = c(dxy.rus.aet.va$avg_dxy,dxy.rus.aet.is$avg_dxy,dxy.nc.rus.aet.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.smi.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.smi.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.smi.is$avg_dxy)))
dxy.rus.smi <- data.frame(
  dxy = c(dxy.rus.smi.va$avg_dxy,dxy.rus.smi.is$avg_dxy,dxy.nc.rus.smi.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.neo.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.neo.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.neo.is$avg_dxy)))
dxy.rus.neo <- data.frame(
  dxy = c(dxy.rus.neo.va$avg_dxy,dxy.rus.neo.is$avg_dxy,dxy.nc.rus.neo.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.jav.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.jav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.jav.is$avg_dxy)))
dxy.rus.jav <- data.frame(
  dxy = c(dxy.rus.jav.va$avg_dxy,dxy.rus.jav.is$avg_dxy,dxy.nc.rus.jav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.rus.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.rus.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.rus.dim.is$avg_dxy)))
dxy.rus.dim <- data.frame(
  dxy = c(dxy.rus.dim.va$avg_dxy,dxy.rus.dim.is$avg_dxy,dxy.nc.rus.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.aet.smi.va$avg_dxy)))
y = c(rep("is",length(dxy.aet.smi.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.aet.smi.is$avg_dxy)))
dxy.aet.smi <- data.frame(
  dxy = c(dxy.aet.smi.va$avg_dxy,dxy.aet.smi.is$avg_dxy,dxy.nc.aet.smi.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.aet.neo.va$avg_dxy)))
y = c(rep("is",length(dxy.aet.neo.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.aet.neo.is$avg_dxy)))
dxy.aet.neo <- data.frame(
  dxy = c(dxy.aet.neo.va$avg_dxy,dxy.aet.neo.is$avg_dxy,dxy.nc.aet.neo.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.aet.jav.va$avg_dxy)))
y = c(rep("is",length(dxy.aet.jav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.aet.jav.is$avg_dxy)))
dxy.aet.jav <- data.frame(
  dxy = c(dxy.aet.jav.va$avg_dxy,dxy.aet.jav.is$avg_dxy,dxy.nc.aet.jav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.aet.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.aet.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.aet.dim.is$avg_dxy)))
dxy.aet.dim <- data.frame(
  dxy = c(dxy.aet.dim.va$avg_dxy,dxy.aet.dim.is$avg_dxy,dxy.nc.aet.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.smi.neo.va$avg_dxy)))
y = c(rep("is",length(dxy.smi.neo.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.smi.neo.is$avg_dxy)))
dxy.smi.neo <- data.frame(
  dxy = c(dxy.smi.neo.va$avg_dxy,dxy.smi.neo.is$avg_dxy,dxy.nc.smi.neo.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.smi.jav.va$avg_dxy)))
y = c(rep("is",length(dxy.smi.jav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.smi.jav.is$avg_dxy)))
dxy.smi.jav <- data.frame(
  dxy = c(dxy.smi.jav.va$avg_dxy,dxy.smi.jav.is$avg_dxy,dxy.nc.smi.jav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.smi.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.smi.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.smi.dim.is$avg_dxy)))
dxy.smi.dim <- data.frame(
  dxy = c(dxy.smi.dim.va$avg_dxy,dxy.smi.dim.is$avg_dxy,dxy.nc.smi.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.neo.jav.va$avg_dxy)))
y = c(rep("is",length(dxy.neo.jav.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.neo.jav.is$avg_dxy)))
dxy.neo.jav <- data.frame(
  dxy = c(dxy.neo.jav.va$avg_dxy,dxy.neo.jav.is$avg_dxy,dxy.nc.neo.jav.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.neo.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.neo.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.neo.dim.is$avg_dxy)))
dxy.neo.dim <- data.frame(
  dxy = c(dxy.neo.dim.va$avg_dxy,dxy.neo.dim.is$avg_dxy,dxy.nc.neo.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(dxy.jav.dim.va$avg_dxy)))
y = c(rep("is",length(dxy.jav.dim.is$avg_dxy)))
z = c(rep("nc",length(dxy.nc.jav.dim.is$avg_dxy)))
dxy.jav.dim <- data.frame(
  dxy = c(dxy.jav.dim.va$avg_dxy,dxy.jav.dim.is$avg_dxy,dxy.nc.jav.dim.is$avg_dxy),
  loci = factor(c(x,y,z))
)

## Perform ANOVA
aov.dxy.rus.sav <- aov(dxy ~ loci, data = dxy.rus.sav)
aov.dxy.rus.tyt <- aov(dxy ~ loci, data = dxy.rus.tyt)
aov.dxy.rus.aet <- aov(dxy ~ loci, data = dxy.rus.aet)
aov.dxy.rus.smi <- aov(dxy ~ loci, data = dxy.rus.smi)
aov.dxy.rus.neo <- aov(dxy ~ loci, data = dxy.rus.neo)
aov.dxy.rus.jav <- aov(dxy ~ loci, data = dxy.rus.jav)
aov.dxy.rus.dim <- aov(dxy ~ loci, data = dxy.rus.dim)
aov.dxy.aet.smi <- aov(dxy ~ loci, data = dxy.aet.smi)
aov.dxy.aet.neo <- aov(dxy ~ loci, data = dxy.aet.neo)
aov.dxy.aet.jav <- aov(dxy ~ loci, data = dxy.aet.jav)
aov.dxy.aet.dim <- aov(dxy ~ loci, data = dxy.aet.dim)
aov.dxy.smi.neo <- aov(dxy ~ loci, data = dxy.smi.neo)
aov.dxy.smi.jav <- aov(dxy ~ loci, data = dxy.smi.jav)
aov.dxy.smi.dim <- aov(dxy ~ loci, data = dxy.smi.dim)
aov.dxy.neo.jav <- aov(dxy ~ loci, data = dxy.neo.jav)
aov.dxy.neo.dim <- aov(dxy ~ loci, data = dxy.neo.dim)
aov.dxy.jav.dim <- aov(dxy ~ loci, data = dxy.jav.dim)

summary(aov.dxy.rus.sav)
summary(aov.dxy.rus.tyt)
summary(aov.dxy.rus.aet)
summary(aov.dxy.rus.smi)
summary(aov.dxy.rus.neo)
summary(aov.dxy.rus.jav)
summary(aov.dxy.rus.dim)
summary(aov.dxy.aet.smi)
summary(aov.dxy.aet.neo)
summary(aov.dxy.aet.jav)
summary(aov.dxy.aet.dim)
summary(aov.dxy.smi.neo)
summary(aov.dxy.smi.jav)
summary(aov.dxy.smi.dim)
summary(aov.dxy.neo.jav)
summary(aov.dxy.neo.dim)
summary(aov.dxy.jav.dim)

## Perform Tukey post-hoc
post.dxy.rus.sav <- glht(aov.dxy.rus.sav, linfct = mcp(loci = "Tukey"))
post.dxy.rus.tyt <- glht(aov.dxy.rus.tyt, linfct = mcp(loci = "Tukey"))
post.dxy.rus.aet <- glht(aov.dxy.rus.aet, linfct = mcp(loci = "Tukey"))
post.dxy.rus.smi <- glht(aov.dxy.rus.smi, linfct = mcp(loci = "Tukey"))
post.dxy.rus.neo <- glht(aov.dxy.rus.neo, linfct = mcp(loci = "Tukey"))
post.dxy.rus.jav <- glht(aov.dxy.rus.jav, linfct = mcp(loci = "Tukey"))
post.dxy.rus.dim <- glht(aov.dxy.rus.dim, linfct = mcp(loci = "Tukey"))
post.dxy.aet.smi <- glht(aov.dxy.aet.smi, linfct = mcp(loci = "Tukey"))
post.dxy.aet.neo <- glht(aov.dxy.aet.neo, linfct = mcp(loci = "Tukey"))
post.dxy.aet.jav <- glht(aov.dxy.aet.jav, linfct = mcp(loci = "Tukey"))
post.dxy.aet.dim <- glht(aov.dxy.aet.dim, linfct = mcp(loci = "Tukey"))
post.dxy.smi.neo <- glht(aov.dxy.smi.neo, linfct = mcp(loci = "Tukey"))
post.dxy.smi.jav <- glht(aov.dxy.smi.jav, linfct = mcp(loci = "Tukey"))
post.dxy.smi.dim <- glht(aov.dxy.smi.dim, linfct = mcp(loci = "Tukey"))
post.dxy.neo.jav <- glht(aov.dxy.neo.jav, linfct = mcp(loci = "Tukey"))
post.dxy.neo.dim <- glht(aov.dxy.neo.dim, linfct = mcp(loci = "Tukey"))
post.dxy.jav.dim <- glht(aov.dxy.jav.dim, linfct = mcp(loci = "Tukey"))

summary(post.dxy.rus.sav)
summary(post.dxy.rus.tyt)
summary(post.dxy.rus.aet)
summary(post.dxy.rus.smi)
summary(post.dxy.rus.neo)
summary(post.dxy.rus.jav)
summary(post.dxy.rus.dim)
summary(post.dxy.aet.smi)
summary(post.dxy.aet.neo)
summary(post.dxy.aet.jav)
summary(post.dxy.aet.dim)
summary(post.dxy.smi.neo)
summary(post.dxy.smi.jav)
summary(post.dxy.smi.dim)
summary(post.dxy.neo.jav)
summary(post.dxy.neo.dim)
summary(post.dxy.jav.dim)

#### Perform ANOVA and Tukey post-hoc tests (pi)----

## Calculate mean pi
pi.rus.sav.is <- (pi.rus.rus.sav.is$avg_pi + pi.sav.rus.sav.is$avg_pi)/2
pi.rus.tyt.is <- (pi.rus.rus.tyt.is$avg_pi + pi.tyt.rus.tyt.is$avg_pi)/2
pi.rus.aet.is <- (pi.rus.rus.aet.is$avg_pi + pi.aet.rus.aet.is$avg_pi)/2
pi.rus.smi.is <- (pi.rus.rus.smi.is$avg_pi + pi.smi.rus.smi.is$avg_pi)/2
pi.rus.neo.is <- (pi.rus.rus.neo.is$avg_pi + pi.neo.rus.neo.is$avg_pi)/2
pi.rus.jav.is <- (pi.rus.rus.jav.is$avg_pi + pi.jav.rus.jav.is$avg_pi)/2
pi.rus.dim.is <- (pi.rus.rus.dim.is$avg_pi + pi.dim.rus.dim.is$avg_pi)/2
pi.aet.smi.is <- (pi.aet.aet.smi.is$avg_pi + pi.smi.aet.smi.is$avg_pi)/2
pi.aet.neo.is <- (pi.aet.aet.neo.is$avg_pi + pi.neo.aet.neo.is$avg_pi)/2
pi.aet.jav.is <- (pi.aet.aet.jav.is$avg_pi + pi.jav.aet.jav.is$avg_pi)/2
pi.aet.dim.is <- (pi.aet.aet.dim.is$avg_pi + pi.dim.aet.dim.is$avg_pi)/2
pi.smi.neo.is <- (pi.smi.smi.neo.is$avg_pi + pi.neo.smi.neo.is$avg_pi)/2
pi.smi.jav.is <- (pi.smi.smi.jav.is$avg_pi + pi.jav.smi.jav.is$avg_pi)/2
pi.smi.dim.is <- (pi.smi.smi.dim.is$avg_pi + pi.dim.smi.dim.is$avg_pi)/2
pi.neo.jav.is <- (pi.neo.neo.jav.is$avg_pi + pi.jav.neo.jav.is$avg_pi)/2
pi.neo.dim.is <- (pi.neo.neo.dim.is$avg_pi + pi.dim.neo.dim.is$avg_pi)/2
pi.jav.dim.is <- (pi.jav.jav.dim.is$avg_pi + pi.dim.jav.dim.is$avg_pi)/2

pi.nc.rus.sav.is <- (pi.nc.rus.rus.sav.is$avg_pi + pi.nc.sav.rus.sav.is$avg_pi)/2
pi.nc.rus.tyt.is <- (pi.nc.rus.rus.tyt.is$avg_pi + pi.nc.tyt.rus.tyt.is$avg_pi)/2
pi.nc.rus.aet.is <- (pi.nc.rus.rus.aet.is$avg_pi + pi.nc.aet.rus.aet.is$avg_pi)/2
pi.nc.rus.smi.is <- (pi.nc.rus.rus.smi.is$avg_pi + pi.nc.smi.rus.smi.is$avg_pi)/2
pi.nc.rus.neo.is <- (pi.nc.rus.rus.neo.is$avg_pi + pi.nc.neo.rus.neo.is$avg_pi)/2
pi.nc.rus.jav.is <- (pi.nc.rus.rus.jav.is$avg_pi + pi.nc.jav.rus.jav.is$avg_pi)/2
pi.nc.rus.dim.is <- (pi.nc.rus.rus.dim.is$avg_pi + pi.nc.dim.rus.dim.is$avg_pi)/2
pi.nc.aet.smi.is <- (pi.nc.aet.aet.smi.is$avg_pi + pi.nc.smi.aet.smi.is$avg_pi)/2
pi.nc.aet.neo.is <- (pi.nc.aet.aet.neo.is$avg_pi + pi.nc.neo.aet.neo.is$avg_pi)/2
pi.nc.aet.jav.is <- (pi.nc.aet.aet.jav.is$avg_pi + pi.nc.jav.aet.jav.is$avg_pi)/2
pi.nc.aet.dim.is <- (pi.nc.aet.aet.dim.is$avg_pi + pi.nc.dim.aet.dim.is$avg_pi)/2
pi.nc.smi.neo.is <- (pi.nc.smi.smi.neo.is$avg_pi + pi.nc.neo.smi.neo.is$avg_pi)/2
pi.nc.smi.jav.is <- (pi.nc.smi.smi.jav.is$avg_pi + pi.nc.jav.smi.jav.is$avg_pi)/2
pi.nc.smi.dim.is <- (pi.nc.smi.smi.dim.is$avg_pi + pi.nc.dim.smi.dim.is$avg_pi)/2
pi.nc.neo.jav.is <- (pi.nc.neo.neo.jav.is$avg_pi + pi.nc.jav.neo.jav.is$avg_pi)/2
pi.nc.neo.dim.is <- (pi.nc.neo.neo.dim.is$avg_pi + pi.nc.dim.neo.dim.is$avg_pi)/2
pi.nc.jav.dim.is <- (pi.nc.jav.jav.dim.is$avg_pi + pi.nc.dim.jav.dim.is$avg_pi)/2

pi.rus.sav.va <- (pi.rus.rus.sav.va$avg_pi + pi.sav.rus.sav.va$avg_pi)/2
pi.rus.tyt.va <- (pi.rus.rus.tyt.va$avg_pi + pi.tyt.rus.tyt.va$avg_pi)/2
pi.rus.aet.va <- (pi.rus.rus.aet.va$avg_pi + pi.aet.rus.aet.va$avg_pi)/2
pi.rus.smi.va <- (pi.rus.rus.smi.va$avg_pi + pi.smi.rus.smi.va$avg_pi)/2
pi.rus.neo.va <- (pi.rus.rus.neo.va$avg_pi + pi.neo.rus.neo.va$avg_pi)/2
pi.rus.jav.va <- (pi.rus.rus.jav.va$avg_pi + pi.jav.rus.jav.va$avg_pi)/2
pi.rus.dim.va <- (pi.rus.rus.dim.va$avg_pi + pi.dim.rus.dim.va$avg_pi)/2
pi.aet.smi.va <- (pi.aet.aet.smi.va$avg_pi + pi.smi.aet.smi.va$avg_pi)/2
pi.aet.neo.va <- (pi.aet.aet.neo.va$avg_pi + pi.neo.aet.neo.va$avg_pi)/2
pi.aet.jav.va <- (pi.aet.aet.jav.va$avg_pi + pi.jav.aet.jav.va$avg_pi)/2
pi.aet.dim.va <- (pi.aet.aet.dim.va$avg_pi + pi.dim.aet.dim.va$avg_pi)/2
pi.smi.neo.va <- (pi.smi.smi.neo.va$avg_pi + pi.neo.smi.neo.va$avg_pi)/2
pi.smi.jav.va <- (pi.smi.smi.jav.va$avg_pi + pi.jav.smi.jav.va$avg_pi)/2
pi.smi.dim.va <- (pi.smi.smi.dim.va$avg_pi + pi.dim.smi.dim.va$avg_pi)/2
pi.neo.jav.va <- (pi.neo.neo.jav.va$avg_pi + pi.jav.neo.jav.va$avg_pi)/2
pi.neo.dim.va <- (pi.neo.neo.dim.va$avg_pi + pi.dim.neo.dim.va$avg_pi)/2
pi.jav.dim.va <- (pi.jav.jav.dim.va$avg_pi + pi.dim.jav.dim.va$avg_pi)/2

## Set up data frames for ANOVA
x = c(rep("va",length(pi.rus.sav.va)))
y = c(rep("is",length(pi.rus.sav.is)))
z = c(rep("nc",length(pi.nc.rus.sav.is)))
pi.rus.sav <- data.frame(
  pi = c(pi.rus.sav.va,pi.rus.sav.is,pi.nc.rus.sav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.tyt.va)))
y = c(rep("is",length(pi.rus.tyt.is)))
z = c(rep("nc",length(pi.nc.rus.tyt.is)))
pi.rus.tyt <- data.frame(
  pi = c(pi.rus.tyt.va,pi.rus.tyt.is,pi.nc.rus.tyt.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.aet.va)))
y = c(rep("is",length(pi.rus.aet.is)))
z = c(rep("nc",length(pi.nc.rus.aet.is)))
pi.rus.aet <- data.frame(
  pi = c(pi.rus.aet.va,pi.rus.aet.is,pi.nc.rus.aet.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.smi.va)))
y = c(rep("is",length(pi.rus.smi.is)))
z = c(rep("nc",length(pi.nc.rus.smi.is)))
pi.rus.smi <- data.frame(
  pi = c(pi.rus.smi.va,pi.rus.smi.is,pi.nc.rus.smi.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.neo.va)))
y = c(rep("is",length(pi.rus.neo.is)))
z = c(rep("nc",length(pi.nc.rus.neo.is)))
pi.rus.neo <- data.frame(
  pi = c(pi.rus.neo.va,pi.rus.neo.is,pi.nc.rus.neo.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.jav.va)))
y = c(rep("is",length(pi.rus.jav.is)))
z = c(rep("nc",length(pi.nc.rus.jav.is)))
pi.rus.jav <- data.frame(
  pi = c(pi.rus.jav.va,pi.rus.jav.is,pi.nc.rus.jav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.rus.dim.va)))
y = c(rep("is",length(pi.rus.dim.is)))
z = c(rep("nc",length(pi.nc.rus.dim.is)))
pi.rus.dim <- data.frame(
  pi = c(pi.rus.dim.va,pi.rus.dim.is,pi.nc.rus.dim.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.aet.smi.va)))
y = c(rep("is",length(pi.aet.smi.is)))
z = c(rep("nc",length(pi.nc.aet.smi.is)))
pi.aet.smi <- data.frame(
  pi = c(pi.aet.smi.va,pi.aet.smi.is,pi.nc.aet.smi.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.aet.neo.va)))
y = c(rep("is",length(pi.aet.neo.is)))
z = c(rep("nc",length(pi.nc.aet.neo.is)))
pi.aet.neo <- data.frame(
  pi = c(pi.aet.neo.va,pi.aet.neo.is,pi.nc.aet.neo.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.aet.jav.va)))
y = c(rep("is",length(pi.aet.jav.is)))
z = c(rep("nc",length(pi.nc.aet.jav.is)))
pi.aet.jav <- data.frame(
  pi = c(pi.aet.jav.va,pi.aet.jav.is,pi.nc.aet.jav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.aet.dim.va)))
y = c(rep("is",length(pi.aet.dim.is)))
z = c(rep("nc",length(pi.nc.aet.dim.is)))
pi.aet.dim <- data.frame(
  pi = c(pi.aet.dim.va,pi.aet.dim.is,pi.nc.aet.dim.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.smi.neo.va)))
y = c(rep("is",length(pi.smi.neo.is)))
z = c(rep("nc",length(pi.nc.smi.neo.is)))
pi.smi.neo <- data.frame(
  pi = c(pi.smi.neo.va,pi.smi.neo.is,pi.nc.smi.neo.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.smi.jav.va)))
y = c(rep("is",length(pi.smi.jav.is)))
z = c(rep("nc",length(pi.nc.smi.jav.is)))
pi.smi.jav <- data.frame(
  pi = c(pi.smi.jav.va,pi.smi.jav.is,pi.nc.smi.jav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.smi.dim.va)))
y = c(rep("is",length(pi.smi.dim.is)))
z = c(rep("nc",length(pi.nc.smi.dim.is)))
pi.smi.dim <- data.frame(
  pi = c(pi.smi.dim.va,pi.smi.dim.is,pi.nc.smi.dim.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.neo.jav.va)))
y = c(rep("is",length(pi.neo.jav.is)))
z = c(rep("nc",length(pi.nc.neo.jav.is)))
pi.neo.jav <- data.frame(
  pi = c(pi.neo.jav.va,pi.neo.jav.is,pi.nc.neo.jav.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.neo.dim.va)))
y = c(rep("is",length(pi.neo.dim.is)))
z = c(rep("nc",length(pi.nc.neo.dim.is)))
pi.neo.dim <- data.frame(
  pi = c(pi.neo.dim.va,pi.neo.dim.is,pi.nc.neo.dim.is),
  loci = factor(c(x,y,z))
)

x = c(rep("va",length(pi.jav.dim.va)))
y = c(rep("is",length(pi.jav.dim.is)))
z = c(rep("nc",length(pi.nc.jav.dim.is)))
pi.jav.dim <- data.frame(
  pi = c(pi.jav.dim.va,pi.jav.dim.is,pi.nc.jav.dim.is),
  loci = factor(c(x,y,z))
)

## Perform ANOVA
aov.pi.rus.sav <- aov(pi ~ loci, data = pi.rus.sav)
aov.pi.rus.tyt <- aov(pi ~ loci, data = pi.rus.tyt)
aov.pi.rus.aet <- aov(pi ~ loci, data = pi.rus.aet)
aov.pi.rus.smi <- aov(pi ~ loci, data = pi.rus.smi)
aov.pi.rus.neo <- aov(pi ~ loci, data = pi.rus.neo)
aov.pi.rus.jav <- aov(pi ~ loci, data = pi.rus.jav)
aov.pi.rus.dim <- aov(pi ~ loci, data = pi.rus.dim)
aov.pi.aet.smi <- aov(pi ~ loci, data = pi.aet.smi)
aov.pi.aet.neo <- aov(pi ~ loci, data = pi.aet.neo)
aov.pi.aet.jav <- aov(pi ~ loci, data = pi.aet.jav)
aov.pi.aet.dim <- aov(pi ~ loci, data = pi.aet.dim)
aov.pi.smi.neo <- aov(pi ~ loci, data = pi.smi.neo)
aov.pi.smi.jav <- aov(pi ~ loci, data = pi.smi.jav)
aov.pi.smi.dim <- aov(pi ~ loci, data = pi.smi.dim)
aov.pi.neo.jav <- aov(pi ~ loci, data = pi.neo.jav)
aov.pi.neo.dim <- aov(pi ~ loci, data = pi.neo.dim)
aov.pi.jav.dim <- aov(pi ~ loci, data = pi.jav.dim)

summary(aov.pi.rus.sav)
summary(aov.pi.rus.tyt)
summary(aov.pi.rus.aet)
summary(aov.pi.rus.smi)
summary(aov.pi.rus.neo)
summary(aov.pi.rus.jav)
summary(aov.pi.rus.dim)
summary(aov.pi.aet.smi)
summary(aov.pi.aet.neo)
summary(aov.pi.aet.jav)
summary(aov.pi.aet.dim)
summary(aov.pi.smi.neo)
summary(aov.pi.smi.jav)
summary(aov.pi.smi.dim)
summary(aov.pi.neo.jav)
summary(aov.pi.neo.dim)
summary(aov.pi.jav.dim)

## Perform Tukey post-hoc
post.pi.rus.sav <- glht(aov.pi.rus.sav, linfct = mcp(loci = "Tukey"))
post.pi.rus.tyt <- glht(aov.pi.rus.tyt, linfct = mcp(loci = "Tukey"))
post.pi.rus.aet <- glht(aov.pi.rus.aet, linfct = mcp(loci = "Tukey"))
post.pi.rus.smi <- glht(aov.pi.rus.smi, linfct = mcp(loci = "Tukey"))
post.pi.rus.neo <- glht(aov.pi.rus.neo, linfct = mcp(loci = "Tukey"))
post.pi.rus.jav <- glht(aov.pi.rus.jav, linfct = mcp(loci = "Tukey"))
post.pi.rus.dim <- glht(aov.pi.rus.dim, linfct = mcp(loci = "Tukey"))
post.pi.aet.smi <- glht(aov.pi.aet.smi, linfct = mcp(loci = "Tukey"))
post.pi.aet.neo <- glht(aov.pi.aet.neo, linfct = mcp(loci = "Tukey"))
post.pi.aet.jav <- glht(aov.pi.aet.jav, linfct = mcp(loci = "Tukey"))
post.pi.aet.dim <- glht(aov.pi.aet.dim, linfct = mcp(loci = "Tukey"))
post.pi.smi.neo <- glht(aov.pi.smi.neo, linfct = mcp(loci = "Tukey"))
post.pi.smi.jav <- glht(aov.pi.smi.jav, linfct = mcp(loci = "Tukey"))
post.pi.smi.dim <- glht(aov.pi.smi.dim, linfct = mcp(loci = "Tukey"))
post.pi.neo.jav <- glht(aov.pi.neo.jav, linfct = mcp(loci = "Tukey"))
post.pi.neo.dim <- glht(aov.pi.neo.dim, linfct = mcp(loci = "Tukey"))
post.pi.jav.dim <- glht(aov.pi.jav.dim, linfct = mcp(loci = "Tukey"))

summary(post.pi.rus.sav)
summary(post.pi.rus.tyt)
summary(post.pi.rus.aet)
summary(post.pi.rus.smi)
summary(post.pi.rus.neo)
summary(post.pi.rus.jav)
summary(post.pi.rus.dim)
summary(post.pi.aet.smi)
summary(post.pi.aet.neo)
summary(post.pi.aet.jav)
summary(post.pi.aet.dim)
summary(post.pi.smi.neo)
summary(post.pi.smi.jav)
summary(post.pi.smi.dim)
summary(post.pi.neo.jav)
summary(post.pi.neo.dim)
summary(post.pi.jav.dim)

### Get mean and standard deviation Fst, dxy, and pi values----

fst.va.m <- c()
fst.va.m <- append(fst.va.m,mean(fst.rus.sav.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.tyt.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.aet.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.smi.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.neo.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.jav.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.rus.dim.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.aet.smi.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.aet.neo.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.aet.jav.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.aet.dim.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.smi.neo.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.smi.jav.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.smi.dim.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.neo.jav.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.neo.dim.va$avg_wc_fst))
fst.va.m <- append(fst.va.m,mean(fst.jav.dim.va$avg_wc_fst))

for(i in fst.va.m){print(i)}

fst.va.s <- c()
fst.va.s <- append(fst.va.s,sd(fst.rus.sav.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.tyt.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.aet.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.smi.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.neo.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.jav.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.rus.dim.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.aet.smi.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.aet.neo.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.aet.jav.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.aet.dim.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.smi.neo.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.smi.jav.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.smi.dim.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.neo.jav.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.neo.dim.va$avg_wc_fst))
fst.va.s <- append(fst.va.s,sd(fst.jav.dim.va$avg_wc_fst))

for(i in fst.va.s){print(i)}

dxy.va.m <- c()
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.sav.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.tyt.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.aet.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.smi.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.neo.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.jav.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.rus.dim.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.aet.smi.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.aet.neo.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.aet.jav.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.aet.dim.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.smi.neo.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.smi.jav.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.smi.dim.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.neo.jav.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.neo.dim.va$avg_dxy))
dxy.va.m <- append(dxy.va.m,mean(dxy.jav.dim.va$avg_dxy))

for(i in dxy.va.m){print(i)}

dxy.va.s <- c()
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.sav.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.tyt.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.aet.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.smi.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.neo.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.jav.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.rus.dim.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.aet.smi.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.aet.neo.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.aet.jav.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.aet.dim.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.smi.neo.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.smi.jav.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.smi.dim.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.neo.jav.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.neo.dim.va$avg_dxy))
dxy.va.s <- append(dxy.va.s,sd(dxy.jav.dim.va$avg_dxy))

for(i in dxy.va.s){print(i)}

pi.va.m <- c()
pi.va.m <- append(pi.va.m,mean(pi.rus.sav.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.tyt.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.aet.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.smi.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.neo.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.jav.va))
pi.va.m <- append(pi.va.m,mean(pi.rus.dim.va))
pi.va.m <- append(pi.va.m,mean(pi.aet.smi.va))
pi.va.m <- append(pi.va.m,mean(pi.aet.neo.va))
pi.va.m <- append(pi.va.m,mean(pi.aet.jav.va))
pi.va.m <- append(pi.va.m,mean(pi.aet.dim.va))
pi.va.m <- append(pi.va.m,mean(pi.smi.neo.va))
pi.va.m <- append(pi.va.m,mean(pi.smi.jav.va))
pi.va.m <- append(pi.va.m,mean(pi.smi.dim.va))
pi.va.m <- append(pi.va.m,mean(pi.neo.jav.va))
pi.va.m <- append(pi.va.m,mean(pi.neo.dim.va))
pi.va.m <- append(pi.va.m,mean(pi.jav.dim.va))

for(i in pi.va.m){print(i)}

pi.va.s <- c()
pi.va.s <- append(pi.va.s,sd(pi.rus.sav.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.tyt.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.aet.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.smi.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.neo.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.jav.va))
pi.va.s <- append(pi.va.s,sd(pi.rus.dim.va))
pi.va.s <- append(pi.va.s,sd(pi.aet.smi.va))
pi.va.s <- append(pi.va.s,sd(pi.aet.neo.va))
pi.va.s <- append(pi.va.s,sd(pi.aet.jav.va))
pi.va.s <- append(pi.va.s,sd(pi.aet.dim.va))
pi.va.s <- append(pi.va.s,sd(pi.smi.neo.va))
pi.va.s <- append(pi.va.s,sd(pi.smi.jav.va))
pi.va.s <- append(pi.va.s,sd(pi.smi.dim.va))
pi.va.s <- append(pi.va.s,sd(pi.neo.jav.va))
pi.va.s <- append(pi.va.s,sd(pi.neo.dim.va))
pi.va.s <- append(pi.va.s,sd(pi.jav.dim.va))

for(i in pi.va.s){print(i)}

fst.is.m <- c()
fst.is.m <- append(fst.is.m,mean(fst.rus.sav.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.tyt.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.aet.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.smi.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.neo.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.jav.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.rus.dim.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.aet.smi.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.aet.neo.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.aet.jav.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.aet.dim.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.smi.neo.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.smi.jav.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.smi.dim.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.neo.jav.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.neo.dim.is$avg_wc_fst))
fst.is.m <- append(fst.is.m,mean(fst.jav.dim.is$avg_wc_fst))

for(i in fst.is.m){print(i)}

fst.is.s <- c()
fst.is.s <- append(fst.is.s,sd(fst.rus.sav.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.tyt.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.aet.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.smi.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.neo.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.jav.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.rus.dim.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.aet.smi.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.aet.neo.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.aet.jav.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.aet.dim.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.smi.neo.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.smi.jav.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.smi.dim.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.neo.jav.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.neo.dim.is$avg_wc_fst))
fst.is.s <- append(fst.is.s,sd(fst.jav.dim.is$avg_wc_fst))

for(i in fst.is.s){print(i)}

dxy.is.m <- c()
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.sav.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.tyt.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.aet.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.smi.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.neo.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.jav.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.rus.dim.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.aet.smi.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.aet.neo.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.aet.jav.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.aet.dim.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.smi.neo.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.smi.jav.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.smi.dim.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.neo.jav.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.neo.dim.is$avg_dxy))
dxy.is.m <- append(dxy.is.m,mean(dxy.jav.dim.is$avg_dxy))

for(i in dxy.is.m){print(i)}

dxy.is.s <- c()
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.sav.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.tyt.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.aet.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.smi.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.neo.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.jav.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.rus.dim.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.aet.smi.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.aet.neo.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.aet.jav.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.aet.dim.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.smi.neo.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.smi.jav.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.smi.dim.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.neo.jav.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.neo.dim.is$avg_dxy))
dxy.is.s <- append(dxy.is.s,sd(dxy.jav.dim.is$avg_dxy))

for(i in dxy.is.s){print(i)}

pi.is.m <- c()
pi.is.m <- append(pi.is.m,mean(pi.rus.sav.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.tyt.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.aet.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.smi.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.neo.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.jav.is))
pi.is.m <- append(pi.is.m,mean(pi.rus.dim.is))
pi.is.m <- append(pi.is.m,mean(pi.aet.smi.is))
pi.is.m <- append(pi.is.m,mean(pi.aet.neo.is))
pi.is.m <- append(pi.is.m,mean(pi.aet.jav.is))
pi.is.m <- append(pi.is.m,mean(pi.aet.dim.is))
pi.is.m <- append(pi.is.m,mean(pi.smi.neo.is))
pi.is.m <- append(pi.is.m,mean(pi.smi.jav.is))
pi.is.m <- append(pi.is.m,mean(pi.smi.dim.is))
pi.is.m <- append(pi.is.m,mean(pi.neo.jav.is))
pi.is.m <- append(pi.is.m,mean(pi.neo.dim.is))
pi.is.m <- append(pi.is.m,mean(pi.jav.dim.is))

for(i in pi.is.m){print(i)}

pi.is.s <- c()
pi.is.s <- append(pi.is.s,sd(pi.rus.sav.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.tyt.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.aet.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.smi.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.neo.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.jav.is))
pi.is.s <- append(pi.is.s,sd(pi.rus.dim.is))
pi.is.s <- append(pi.is.s,sd(pi.aet.smi.is))
pi.is.s <- append(pi.is.s,sd(pi.aet.neo.is))
pi.is.s <- append(pi.is.s,sd(pi.aet.jav.is))
pi.is.s <- append(pi.is.s,sd(pi.aet.dim.is))
pi.is.s <- append(pi.is.s,sd(pi.smi.neo.is))
pi.is.s <- append(pi.is.s,sd(pi.smi.jav.is))
pi.is.s <- append(pi.is.s,sd(pi.smi.dim.is))
pi.is.s <- append(pi.is.s,sd(pi.neo.jav.is))
pi.is.s <- append(pi.is.s,sd(pi.neo.dim.is))
pi.is.s <- append(pi.is.s,sd(pi.jav.dim.is))

for(i in pi.is.s){print(i)}

fst.nc.is.m <- c()
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.sav.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.tyt.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.aet.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.smi.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.neo.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.jav.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.rus.dim.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.aet.smi.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.aet.neo.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.aet.jav.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.aet.dim.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.smi.neo.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.smi.jav.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.smi.dim.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.neo.jav.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.neo.dim.is$avg_wc_fst))
fst.nc.is.m <- append(fst.nc.is.m,mean(fst.nc.jav.dim.is$avg_wc_fst))

for(i in fst.nc.is.m){print(i)}

fst.nc.is.s <- c()
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.sav.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.tyt.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.aet.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.smi.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.neo.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.jav.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.rus.dim.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.aet.smi.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.aet.neo.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.aet.jav.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.aet.dim.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.smi.neo.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.smi.jav.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.smi.dim.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.neo.jav.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.neo.dim.is$avg_wc_fst))
fst.nc.is.s <- append(fst.nc.is.s,sd(fst.nc.jav.dim.is$avg_wc_fst))

for(i in fst.nc.is.s){print(i)}

dxy.nc.is.m <- c()
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.sav.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.tyt.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.aet.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.smi.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.neo.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.jav.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.rus.dim.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.aet.smi.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.aet.neo.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.aet.jav.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.aet.dim.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.smi.neo.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.smi.jav.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.smi.dim.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.neo.jav.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.neo.dim.is$avg_dxy))
dxy.nc.is.m <- append(dxy.nc.is.m,mean(dxy.nc.jav.dim.is$avg_dxy))

for(i in dxy.nc.is.m){print(i)}

dxy.nc.is.s <- c()
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.sav.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.tyt.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.aet.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.smi.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.neo.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.jav.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.rus.dim.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.aet.smi.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.aet.neo.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.aet.jav.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.aet.dim.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.smi.neo.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.smi.jav.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.smi.dim.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.neo.jav.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.neo.dim.is$avg_dxy))
dxy.nc.is.s <- append(dxy.nc.is.s,sd(dxy.nc.jav.dim.is$avg_dxy))

for(i in dxy.nc.is.s){print(i)}

pi.nc.is.m <- c()
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.sav.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.tyt.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.aet.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.smi.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.neo.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.jav.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.rus.dim.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.aet.smi.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.aet.neo.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.aet.jav.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.aet.dim.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.smi.neo.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.smi.jav.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.smi.dim.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.neo.jav.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.neo.dim.is))
pi.nc.is.m <- append(pi.nc.is.m,mean(pi.nc.jav.dim.is))

for(i in pi.nc.is.m){print(i)}

pi.nc.is.s <- c()
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.sav.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.tyt.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.aet.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.smi.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.neo.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.jav.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.rus.dim.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.aet.smi.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.aet.neo.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.aet.jav.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.aet.dim.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.smi.neo.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.smi.jav.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.smi.dim.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.neo.jav.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.neo.dim.is))
pi.nc.is.s <- append(pi.nc.is.s,sd(pi.nc.jav.dim.is))

for(i in pi.nc.is.s){print(i)}
