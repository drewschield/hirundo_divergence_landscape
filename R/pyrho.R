# Hirundo genomic divergence landscape analysis
# pyrho.R - plot & summarize recombination rate results from pyrho

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pyrho/')

### Read in data--------------------------

hrrho <- read.table('./results/recombination.rustica.rmap.1Mb-100kb.txt',header=T)
harho <- read.table('./results/recombination.aethiopica.rmap.1Mb-100kb.txt',header=T)
hsrho <- read.table('./results/recombination.smithii.rmap.1Mb-100kb.txt',header=T)
hnrho <- read.table('./results/recombination.neoxena.rmap.1Mb-100kb.txt',header=T)
htrho <- read.table('./results/recombination.tahitica.rmap.1Mb-100kb.txt',header=T)
hdrho <- read.table('./results/recombination.dimidiata.rmap.1Mb-100kb.txt',header=T)

## Subset chr 1A, chr 4, and Z chromosome

hrrho.1a <- hrrho %>% filter(str_detect(chrom, 'NC_053453.1'))
hrrho.4 <- hrrho %>% filter(str_detect(chrom, 'NC_053454.1'))
hrrho.z <- hrrho %>% filter(str_detect(chrom, 'NC_053488.1'))

harho.1a <- harho %>% filter(str_detect(chrom, 'NC_053453.1'))
harho.4 <- harho %>% filter(str_detect(chrom, 'NC_053454.1'))
harho.z <- harho %>% filter(str_detect(chrom, 'NC_053488.1'))

hsrho.1a <- hsrho %>% filter(str_detect(chrom, 'NC_053453.1'))
hsrho.4 <- hsrho %>% filter(str_detect(chrom, 'NC_053454.1'))
hsrho.z <- hsrho %>% filter(str_detect(chrom, 'NC_053488.1'))

hnrho.1a <- hnrho %>% filter(str_detect(chrom, 'NC_053453.1'))
hnrho.4 <- hnrho %>% filter(str_detect(chrom, 'NC_053454.1'))
hnrho.z <- hnrho %>% filter(str_detect(chrom, 'NC_053488.1'))

htrho.1a <- htrho %>% filter(str_detect(chrom, 'NC_053453.1'))
htrho.4 <- htrho %>% filter(str_detect(chrom, 'NC_053454.1'))
htrho.z <- htrho %>% filter(str_detect(chrom, 'NC_053488.1'))

hdrho.1a <- hdrho %>% filter(str_detect(chrom, 'NC_053453.1'))
hdrho.4 <- hdrho %>% filter(str_detect(chrom, 'NC_053454.1'))
hdrho.z <- hdrho %>% filter(str_detect(chrom, 'NC_053488.1'))

### Plot recombination rate scans----------------

## Scale x-axis equal to longest of the three chromosomes (Z)

par(mfrow=c(1,3))
plot(hrrho.1a$start,hrrho.1a$rate,type='l',col='maroon',lwd=2,xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',main='Chromsome 1A')
lines(harho.1a$start,harho.1a$rate,col='goldenrod3',lwd=2)
lines(hsrho.1a$start,hsrho.1a$rate,col='darkorange3',lwd=2)
lines(hnrho.1a$start,hnrho.1a$rate,col='aquamarine4',lwd=2)
lines(htrho.1a$start,htrho.1a$rate,col='darkgreen',lwd=2)
lines(hdrho.1a$start,hdrho.1a$rate,col='skyblue',lwd=2)

plot(hrrho.4$start,hrrho.4$rate,type='l',col='maroon',lwd=2,xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',main='Chromsome 4')
lines(harho.4$start,harho.4$rate,col='goldenrod3',lwd=2)
lines(hsrho.4$start,hsrho.4$rate,col='darkorange3',lwd=2)
lines(hnrho.4$start,hnrho.4$rate,col='aquamarine4',lwd=2)
lines(htrho.4$start,htrho.4$rate,col='darkgreen',lwd=2)
lines(hdrho.4$start,hdrho.4$rate,col='skyblue',lwd=2)

plot(hrrho.z$start,hrrho.z$rate,type='l',col='maroon',lwd=2,ylim=c(0,2.5e-08),xlim=c(0,9.02e+07),ylab='Recombination Rate',xlab='Chromosome Position',main='Z chromsome')
lines(harho.z$start,harho.z$rate,col='goldenrod3',lwd=2)
lines(hsrho.z$start,hsrho.z$rate,col='darkorange3',lwd=2)
lines(hnrho.z$start,hnrho.z$rate,col='aquamarine4',lwd=2)
lines(htrho.z$start,htrho.z$rate,col='darkgreen',lwd=2)
lines(hdrho.z$start,hdrho.z$rate,col='skyblue',lwd=2)

## Export at 4 x 14

### Correlations between recombination landscapes----

hrrho <- read.table('./results/recombination.rustica.rmap.1mb.txt',header=T)
harho <- read.table('./results/recombination.aethiopica.rmap.1Mb.txt',header=T)
hsrho <- read.table('./results/recombination.smithii.rmap.1Mb.txt',header=T)
hnrho <- read.table('./results/recombination.neoxena.rmap.1Mb.txt',header=T)
htrho <- read.table('./results/recombination.tahitica.rmap.1Mb.txt',header=T)
hdrho <- read.table('./results/recombination.dimidiata.rmap.1Mb.txt',header=T)

par(mfrow=c(3,5))
plot(hrrho$rate,harho$rate,pch=20,col='slategray')
plot(hrrho$rate,hsrho$rate,pch=20,col='slategray')
plot(hrrho$rate,hnrho$rate,pch=20,col='slategray')
plot(hrrho$rate,htrho$rate,pch=20,col='slategray')
plot(hrrho$rate,hdrho$rate,pch=20,col='slategray')
plot(harho$rate,hsrho$rate,pch=20,col='slategray')
plot(harho$rate,hnrho$rate,pch=20,col='slategray')
plot(harho$rate,htrho$rate,pch=20,col='slategray')
plot(harho$rate,hdrho$rate,pch=20,col='slategray')
plot(hsrho$rate,hnrho$rate,pch=20,col='slategray')
plot(hsrho$rate,htrho$rate,pch=20,col='slategray')
plot(hsrho$rate,hdrho$rate,pch=20,col='slategray')
plot(hnrho$rate,htrho$rate,pch=20,col='slategray')
plot(hnrho$rate,hdrho$rate,pch=20,col='slategray')
plot(htrho$rate,hdrho$rate,pch=20,col='slategray')

cor.test(as.numeric(hrrho$rate),as.numeric(harho$rate),method='spearman')
cor.test(as.numeric(hrrho$rate),as.numeric(hsrho$rate),method='spearman')
cor.test(as.numeric(hrrho$rate),as.numeric(hnrho$rate),method='spearman')
cor.test(as.numeric(hrrho$rate),as.numeric(htrho$rate),method='spearman')
cor.test(as.numeric(hrrho$rate),as.numeric(hdrho$rate),method='spearman')
cor.test(as.numeric(harho$rate),as.numeric(hsrho$rate),method='spearman')
cor.test(as.numeric(harho$rate),as.numeric(hnrho$rate),method='spearman')
cor.test(as.numeric(harho$rate),as.numeric(htrho$rate),method='spearman')
cor.test(as.numeric(harho$rate),as.numeric(hdrho$rate),method='spearman')
cor.test(as.numeric(hsrho$rate),as.numeric(hnrho$rate),method='spearman')
cor.test(as.numeric(hsrho$rate),as.numeric(htrho$rate),method='spearman')
cor.test(as.numeric(hsrho$rate),as.numeric(hdrho$rate),method='spearman')
cor.test(as.numeric(hnrho$rate),as.numeric(htrho$rate),method='spearman')
cor.test(as.numeric(hnrho$rate),as.numeric(hdrho$rate),method='spearman')
cor.test(as.numeric(htrho$rate),as.numeric(hdrho$rate),method='spearman')

## All are very strongly positively correlated - no surprise there



