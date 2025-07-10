# Hirundo genomic divergence landscape analysis
# pixy_Fst_peaks_TajimaD_FayWuH.R - summarize distributions of Tajima's D and Fay & Wu's H in differentiation peaks vs the genome background

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)
library(Rmisc)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----

## Fst
fst <- read.table('./results/pixy.all.50kb.fst.txt',header=T)
fst <- fst %>% filter(chromosome != 'NC_053488.1')
colnames(fst)[4] <- "start"
colnames(fst)[5] <- "end"

fst.rus.aet <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
fst.neo.jav <- fst %>% filter((pop1 == 'HN') & (pop2 == 'HT'))
fst.smi.dim <- fst %>% filter((pop1 == 'HS') & (pop2 == 'HD'))

## Tajima's D (read in and format)
taj.rus <- read.table('../tajima/results/tajima.rustica.50kb.txt',header=T)
colnames(taj.rus)[1] <- "chromosome"
colnames(taj.rus)[2] <- "start"
colnames(taj.rus)[3] <- "end"
taj.rus$start <- taj.rus$start+1

taj.aet <- read.table('../tajima/results/tajima.aethiopica.50kb.txt',header=T)
colnames(taj.aet)[1] <- "chromosome"
colnames(taj.aet)[2] <- "start"
colnames(taj.aet)[3] <- "end"
taj.aet$start <- taj.aet$start+1

taj.neo <- read.table('../tajima/results/tajima.neoxena.50kb.txt',header=T)
colnames(taj.neo)[1] <- "chromosome"
colnames(taj.neo)[2] <- "start"
colnames(taj.neo)[3] <- "end"
taj.neo$start <- taj.neo$start+1

taj.jav <- read.table('../tajima/results/tajima.javanica.50kb.txt',header=T)
colnames(taj.jav)[1] <- "chromosome"
colnames(taj.jav)[2] <- "start"
colnames(taj.jav)[3] <- "end"
taj.jav$start <- taj.jav$start+1

taj.smi <- read.table('../tajima/results/tajima.smithii.50kb.txt',header=T)
colnames(taj.smi)[1] <- "chromosome"
colnames(taj.smi)[2] <- "start"
colnames(taj.smi)[3] <- "end"
taj.smi$start <- taj.smi$start+1

taj.dim <- read.table('../tajima/results/tajima.dimidiata.50kb.txt',header=T)
colnames(taj.dim)[1] <- "chromosome"
colnames(taj.dim)[2] <- "start"
colnames(taj.dim)[3] <- "end"
taj.dim$start <- taj.dim$start+1

## Fay & Wu's H
fay.rus <- read.table('../faywu/sfs.rustica.50kb.txt',header=T)
colnames(fay.rus)[1] <- "chromosome"
colnames(fay.rus)[2] <- "start"
colnames(fay.rus)[3] <- "end"
fay.rus$start <- fay.rus$start+1

fay.aet <- read.table('../faywu/sfs.aethiopica.50kb.txt',header=T)
colnames(fay.aet)[1] <- "chromosome"
colnames(fay.aet)[2] <- "start"
colnames(fay.aet)[3] <- "end"
fay.aet$start <- fay.aet$start+1

fay.neo <- read.table('../faywu/sfs.neoxena.50kb.txt',header=T)
colnames(fay.neo)[1] <- "chromosome"
colnames(fay.neo)[2] <- "start"
colnames(fay.neo)[3] <- "end"
fay.neo$start <- fay.neo$start+1

fay.jav <- read.table('../faywu/sfs.javanica.50kb.txt',header=T)
colnames(fay.jav)[1] <- "chromosome"
colnames(fay.jav)[2] <- "start"
colnames(fay.jav)[3] <- "end"
fay.jav$start <- fay.jav$start+1

fay.smi <- read.table('../faywu/sfs.smithii.50kb.txt',header=T)
colnames(fay.smi)[1] <- "chromosome"
colnames(fay.smi)[2] <- "start"
colnames(fay.smi)[3] <- "end"
fay.smi$start <- fay.smi$start+1

fay.dim <- read.table('../faywu/sfs.dimidiata.50kb.txt',header=T)
colnames(fay.dim)[1] <- "chromosome"
colnames(fay.dim)[2] <- "start"
colnames(fay.dim)[3] <- "end"
fay.dim$start <- fay.dim$start+1

## Permutations: Tajima's D
perm.taj.rus <- read.table('../tajima/permutations/permutation.tajima.50kb.rustica.txt',header=T)
perm.taj.aet <- read.table('../tajima/permutations/permutation.tajima.50kb.aethiopica.txt',header=T)
perm.taj.neo <- read.table('../tajima/permutations/permutation.tajima.50kb.neoxena.txt',header=T)
perm.taj.jav <- read.table('../tajima/permutations/permutation.tajima.50kb.javanica.txt',header=T)
perm.taj.smi <- read.table('../tajima/permutations/permutation.tajima.50kb.smithii.txt',header=T)
perm.taj.dim <- read.table('../tajima/permutations/permutation.tajima.50kb.dimidiata.txt',header=T)

## Permutations: Fay & Wu's H
perm.fay.rus <- read.table('../faywu/permutations/permutation.faywu.50kb.rustica.txt',header=T)
perm.fay.aet <- read.table('../faywu/permutations/permutation.faywu.50kb.aethiopica.txt',header=T)
perm.fay.neo <- read.table('../faywu/permutations/permutation.faywu.50kb.neoxena.txt',header=T)
perm.fay.jav <- read.table('../faywu/permutations/permutation.faywu.50kb.javanica.txt',header=T)
perm.fay.smi <- read.table('../faywu/permutations/permutation.faywu.50kb.smithii.txt',header=T)
perm.fay.dim <- read.table('../faywu/permutations/permutation.faywu.50kb.dimidiata.txt',header=T)


### Merge Fst & Tajima's D & Fay & Wu's H data----
comp.rus <- merge(fst.rus.aet, taj.rus, by=c("chromosome","start"),sort=F)
comp.aet <- merge(fst.rus.aet, taj.aet, by=c("chromosome","start"),sort=F)
comp.neo <- merge(fst.neo.jav, taj.neo, by=c("chromosome","start"),sort=F)
comp.jav <- merge(fst.neo.jav, taj.jav, by=c("chromosome","start"),sort=F)
comp.smi <- merge(fst.smi.dim, taj.smi, by=c("chromosome","start"),sort=F)
comp.dim <- merge(fst.smi.dim, taj.dim, by=c("chromosome","start"),sort=F)

## Additional merge for Fay & Wu's H
comp.rus <- merge(comp.rus, fay.rus, by=c("chromosome","start"),sort=F)
comp.aet <- merge(comp.aet, fay.aet, by=c("chromosome","start"),sort=F)
comp.neo <- merge(comp.neo, fay.neo, by=c("chromosome","start"),sort=F)
comp.jav <- merge(comp.jav, fay.jav, by=c("chromosome","start"),sort=F)
comp.smi <- merge(comp.smi, fay.smi, by=c("chromosome","start"),sort=F)
comp.dim <- merge(comp.dim, fay.dim, by=c("chromosome","start"),sort=F)


### Define and extract Fst islands----

## Read in null distributions based on permutations
perm.rus.aet <- read.table('./permutations/permutation.fst.50kb.rus.aet.txt',header=T)
perm.neo.jav <- read.table('./permutations/permutation.fst.50kb.neo.jav.txt',header=T)
perm.smi.dim <- read.table('./permutations/permutation.fst.50kb.smi.dim.txt',header=T)

## Parse Fst for all regions >= null distribution 
fst.rus.aet.is <- fst.rus.aet %>% filter(avg_wc_fst >= max(perm.rus.aet$Fst))
fst.neo.jav.is <- fst.neo.jav %>% filter(avg_wc_fst >= max(perm.neo.jav$Fst))
fst.smi.dim.is <- fst.smi.dim %>% filter(avg_wc_fst >= max(perm.smi.dim$Fst))


### 1. Examine overlap between Fst islands and Tajima's D and Fay & Wu's H outliers 

## Extract Tajima's D + Fay & Wu's H outliers in Fst islands & background
comp.rus.fst.out <- comp.rus %>% filter(avg_wc_fst >= max(perm.rus.aet$Fst) & TajimaD <= min(perm.taj.rus$TajimaD) & FAY_WU_H <= min(perm.fay.rus$FAY_WU_H))
comp.aet.fst.out <- comp.aet %>% filter(avg_wc_fst >= max(perm.rus.aet$Fst) & TajimaD <= min(perm.taj.aet$TajimaD) & FAY_WU_H <= min(perm.fay.aet$FAY_WU_H))
comp.neo.fst.out <- comp.neo %>% filter(avg_wc_fst >= max(perm.neo.jav$Fst) & TajimaD <= min(perm.taj.neo$TajimaD) & FAY_WU_H <= min(perm.fay.neo$FAY_WU_H))
comp.jav.fst.out <- comp.jav %>% filter(avg_wc_fst >= max(perm.neo.jav$Fst) & TajimaD <= min(perm.taj.jav$TajimaD) & FAY_WU_H <= min(perm.fay.jav$FAY_WU_H))
comp.smi.fst.out <- comp.smi %>% filter(avg_wc_fst >= max(perm.smi.dim$Fst) & TajimaD <= min(perm.taj.smi$TajimaD) & FAY_WU_H <= min(perm.fay.smi$FAY_WU_H))
comp.dim.fst.out <- comp.dim %>% filter(avg_wc_fst >= max(perm.smi.dim$Fst) & TajimaD <= min(perm.taj.dim$TajimaD) & FAY_WU_H <= min(perm.fay.dim$FAY_WU_H))

## Extract all D + H outliers (without regard for Fst)
comp.rus.out <- comp.rus %>% filter(TajimaD <= min(perm.taj.rus$TajimaD) & FAY_WU_H <= min(perm.fay.rus$FAY_WU_H))
comp.aet.out <- comp.aet %>% filter(TajimaD <= min(perm.taj.aet$TajimaD) & FAY_WU_H <= min(perm.fay.aet$FAY_WU_H))
comp.neo.out <- comp.neo %>% filter(TajimaD <= min(perm.taj.neo$TajimaD) & FAY_WU_H <= min(perm.fay.neo$FAY_WU_H))
comp.jav.out <- comp.jav %>% filter(TajimaD <= min(perm.taj.jav$TajimaD) & FAY_WU_H <= min(perm.fay.jav$FAY_WU_H))
comp.smi.out <- comp.smi %>% filter(TajimaD <= min(perm.taj.smi$TajimaD) & FAY_WU_H <= min(perm.fay.smi$FAY_WU_H))
comp.dim.out <- comp.dim %>% filter(TajimaD <= min(perm.taj.dim$TajimaD) & FAY_WU_H <= min(perm.fay.dim$FAY_WU_H))

## Extract background (no D/H outliers)
comp.rus.back <- comp.rus %>% filter(TajimaD > min(perm.taj.rus$TajimaD) & FAY_WU_H > min(perm.fay.rus$FAY_WU_H))
comp.aet.back <- comp.aet %>% filter(TajimaD > min(perm.taj.aet$TajimaD) & FAY_WU_H > min(perm.fay.aet$FAY_WU_H))
comp.neo.back <- comp.neo %>% filter(TajimaD > min(perm.taj.neo$TajimaD) & FAY_WU_H > min(perm.fay.neo$FAY_WU_H))
comp.jav.back <- comp.jav %>% filter(TajimaD > min(perm.taj.jav$TajimaD) & FAY_WU_H > min(perm.fay.jav$FAY_WU_H))
comp.smi.back <- comp.smi %>% filter(TajimaD > min(perm.taj.smi$TajimaD) & FAY_WU_H > min(perm.fay.smi$FAY_WU_H))
comp.dim.back <- comp.dim %>% filter(TajimaD > min(perm.taj.dim$TajimaD) & FAY_WU_H > min(perm.fay.dim$FAY_WU_H))

## Get numbers of D/H outliers in Fst islands vs all D/H outliers; also get numbers of Fst islands
length(comp.rus.fst.out$TajimaD)
length(comp.aet.fst.out$TajimaD)
length(comp.neo.fst.out$TajimaD)
length(comp.jav.fst.out$TajimaD)
length(comp.smi.fst.out$TajimaD)
length(comp.dim.fst.out$TajimaD)

length(comp.rus.out$TajimaD)
length(comp.aet.out$TajimaD)
length(comp.neo.out$TajimaD)
length(comp.jav.out$TajimaD)
length(comp.smi.out$TajimaD)
length(comp.dim.out$TajimaD)

length(fst.rus.aet.is$avg_wc_fst)
length(fst.neo.jav.is$avg_wc_fst)
length(fst.smi.dim.is$avg_wc_fst)

## Proportions of Fst outliers with both D & H outliers
prop.out.rus <- length(comp.rus.fst.out$TajimaD)/length(fst.rus.aet.is$avg_wc_fst)
prop.out.aet <- length(comp.aet.fst.out$TajimaD)/length(fst.rus.aet.is$avg_wc_fst)
prop.out.neo <- length(comp.neo.fst.out$TajimaD)/length(fst.neo.jav.is$avg_wc_fst)
prop.out.jav <- length(comp.jav.fst.out$TajimaD)/length(fst.neo.jav.is$avg_wc_fst)
prop.out.smi <- length(comp.smi.fst.out$TajimaD)/length(fst.smi.dim.is$avg_wc_fst)
prop.out.dim <- length(comp.dim.fst.out$TajimaD)/length(fst.smi.dim.is$avg_wc_fst)

## Fisher's exact tests of enrichment of Fst islands for D + H outliers
## Formula for matrices: number of Fst islands with outliers, total Fst islands, total outliers, all windows

test.rus <- matrix(c(63,412,99,19143),nrow=2,ncol=2)
fisher.test(test.rus)

test.aet <- matrix(c(53,412,229,19143),nrow=2,ncol=2)
fisher.test(test.aet)

test.neo <- matrix(c(103,461,139,19143),nrow=2,ncol=2)
fisher.test(test.neo)

test.jav <- matrix(c(13,461,30,19143),nrow=2,ncol=2)
fisher.test(test.jav)

test.smi <- matrix(c(20,1053,271,19143),nrow=2,ncol=2)
fisher.test(test.smi)

test.dim <- matrix(c(201,1053,235,19143),nrow=2,ncol=2)
fisher.test(test.dim)

# All significant except for H. smithii (p = 0.23)


### 2. Compare Fst in Tajima's D and Fay & Wu's H outlier windows vs background

par(mfrow=c(3,2))
boxplot(comp.rus.back$avg_wc_fst,comp.rus.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))
boxplot(comp.aet.back$avg_wc_fst,comp.aet.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))
boxplot(comp.neo.back$avg_wc_fst,comp.neo.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))
boxplot(comp.jav.back$avg_wc_fst,comp.jav.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))
boxplot(comp.smi.back$avg_wc_fst,comp.smi.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))
boxplot(comp.dim.back$avg_wc_fst,comp.dim.out$avg_wc_fst,pch=20,names=c('Background','Outliers'),ylab='FST',ylim=c(0,1))

## Export at 3 x 6.5

wilcox.test(comp.rus.back$avg_wc_fst,comp.rus.out$avg_wc_fst)
wilcox.test(comp.aet.back$avg_wc_fst,comp.aet.out$avg_wc_fst)
wilcox.test(comp.neo.back$avg_wc_fst,comp.neo.out$avg_wc_fst)
wilcox.test(comp.jav.back$avg_wc_fst,comp.jav.out$avg_wc_fst)
wilcox.test(comp.smi.back$avg_wc_fst,comp.smi.out$avg_wc_fst)
wilcox.test(comp.dim.back$avg_wc_fst,comp.dim.out$avg_wc_fst)


### 3. Get numbers of overlapping outlier islands between species
comp.rus.fst.out %>% inner_join(y=comp.aet.fst.out, by=c('chromosome','start','end.x')) # 6
comp.rus.fst.out %>% inner_join(y=comp.neo.fst.out, by=c('chromosome','start','end.x')) # 4
comp.rus.fst.out %>% inner_join(y=comp.jav.fst.out, by=c('chromosome','start','end.x')) # 0
comp.rus.fst.out %>% inner_join(y=comp.smi.fst.out, by=c('chromosome','start','end.x')) # 0
comp.rus.fst.out %>% inner_join(y=comp.dim.fst.out, by=c('chromosome','start','end.x')) # 13

comp.aet.fst.out %>% inner_join(y=comp.rus.fst.out, by=c('chromosome','start','end.x')) # 6
comp.aet.fst.out %>% inner_join(y=comp.neo.fst.out, by=c('chromosome','start','end.x')) # 5
comp.aet.fst.out %>% inner_join(y=comp.jav.fst.out, by=c('chromosome','start','end.x')) # 2
comp.aet.fst.out %>% inner_join(y=comp.smi.fst.out, by=c('chromosome','start','end.x')) # 0
comp.aet.fst.out %>% inner_join(y=comp.dim.fst.out, by=c('chromosome','start','end.x')) # 20

comp.neo.fst.out %>% inner_join(y=comp.rus.fst.out, by=c('chromosome','start','end.x')) # 4
comp.neo.fst.out %>% inner_join(y=comp.aet.fst.out, by=c('chromosome','start','end.x')) # 5
comp.neo.fst.out %>% inner_join(y=comp.jav.fst.out, by=c('chromosome','start','end.x')) # 7
comp.neo.fst.out %>% inner_join(y=comp.smi.fst.out, by=c('chromosome','start','end.x')) # 0
comp.neo.fst.out %>% inner_join(y=comp.dim.fst.out, by=c('chromosome','start','end.x')) # 22

comp.jav.fst.out %>% inner_join(y=comp.rus.fst.out, by=c('chromosome','start','end.x')) # 0
comp.jav.fst.out %>% inner_join(y=comp.aet.fst.out, by=c('chromosome','start','end.x')) # 2
comp.jav.fst.out %>% inner_join(y=comp.neo.fst.out, by=c('chromosome','start','end.x')) # 7
comp.jav.fst.out %>% inner_join(y=comp.smi.fst.out, by=c('chromosome','start','end.x')) # 0
comp.jav.fst.out %>% inner_join(y=comp.dim.fst.out, by=c('chromosome','start','end.x')) # 5

comp.smi.fst.out %>% inner_join(y=comp.rus.fst.out, by=c('chromosome','start','end.x')) # 0
comp.smi.fst.out %>% inner_join(y=comp.aet.fst.out, by=c('chromosome','start','end.x')) # 0
comp.smi.fst.out %>% inner_join(y=comp.neo.fst.out, by=c('chromosome','start','end.x')) # 0
comp.smi.fst.out %>% inner_join(y=comp.jav.fst.out, by=c('chromosome','start','end.x')) # 0
comp.smi.fst.out %>% inner_join(y=comp.dim.fst.out, by=c('chromosome','start','end.x')) # 4

comp.dim.fst.out %>% inner_join(y=comp.rus.fst.out, by=c('chromosome','start','end.x')) # 13
comp.dim.fst.out %>% inner_join(y=comp.aet.fst.out, by=c('chromosome','start','end.x')) # 20
comp.dim.fst.out %>% inner_join(y=comp.neo.fst.out, by=c('chromosome','start','end.x')) # 22
comp.dim.fst.out %>% inner_join(y=comp.jav.fst.out, by=c('chromosome','start','end.x')) # 5
comp.dim.fst.out %>% inner_join(y=comp.smi.fst.out, by=c('chromosome','start','end.x')) # 4


### Extract Tajima's D in peaks vs background----

## Set outliers
q.rus.aet <- max(perm.rus.aet$Fst)
q.neo.jav <- max(perm.neo.jav$Fst)
q.smi.dim <- max(perm.smi.dim$Fst)

## Extract peak and background regions
comp.rus.is <- comp.rus %>% filter(avg_wc_fst >= q.rus.aet)
comp.rus.va <- comp.rus %>% filter(avg_wc_fst < q.rus.aet)
comp.aet.is <- comp.aet %>% filter(avg_wc_fst >= q.rus.aet)
comp.aet.va <- comp.aet %>% filter(avg_wc_fst < q.rus.aet)
comp.neo.is <- comp.neo %>% filter(avg_wc_fst >= q.neo.jav)
comp.neo.va <- comp.neo %>% filter(avg_wc_fst < q.neo.jav)
comp.jav.is <- comp.jav %>% filter(avg_wc_fst >= q.neo.jav)
comp.jav.va <- comp.jav %>% filter(avg_wc_fst < q.neo.jav)
comp.smi.is <- comp.smi %>% filter(avg_wc_fst >= q.smi.dim)
comp.smi.va <- comp.smi %>% filter(avg_wc_fst < q.smi.dim)
comp.dim.is <- comp.dim %>% filter(avg_wc_fst >= q.smi.dim)
comp.dim.va <- comp.dim %>% filter(avg_wc_fst < q.smi.dim)

### Plot Tajima's D distributions----
par(mfrow=c(3,2))
boxplot(comp.rus.va$TajimaD,comp.rus.is$TajimaD,outline=F)
boxplot(comp.aet.va$TajimaD,comp.aet.is$TajimaD,outline=F)
boxplot(comp.neo.va$TajimaD,comp.neo.is$TajimaD,outline=F)
boxplot(comp.jav.va$TajimaD,comp.jav.is$TajimaD,outline=F)
boxplot(comp.smi.va$TajimaD,comp.smi.is$TajimaD,outline=F)
boxplot(comp.dim.va$TajimaD,comp.dim.is$TajimaD,outline=F)

## Comparison with null distributions
#par(mfrow=c(3,2))
#boxplot(perm.taj.rus$TajimaD,comp.rus.is$TajimaD,outline=F)
#boxplot(perm.taj.aet$TajimaD,comp.aet.is$TajimaD,outline=F)
#boxplot(perm.taj.neo$TajimaD,comp.neo.is$TajimaD,outline=F)
#boxplot(perm.taj.jav$TajimaD,comp.jav.is$TajimaD,outline=F)
#boxplot(perm.taj.smi$TajimaD,comp.smi.is$TajimaD,outline=F)
#boxplot(perm.taj.dim$TajimaD,comp.dim.is$TajimaD,outline=F)

plot(density(comp.rus.is$TajimaD),ylim=c(0,2.2),xlim=c(-3,2),col='seagreen',main='rustica',xlab="Tajima's D")
lines(density(comp.rus.va$TajimaD))
abline(v=quantile(comp.rus$TajimaD,c(0.05)),col='red')

plot(density(comp.aet.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='aethiopica',xlab="Tajima's D")
lines(density(comp.aet.va$TajimaD))
abline(v=quantile(comp.aet$TajimaD,c(0.05)),col='red')

plot(density(comp.neo.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='neoxena',xlab="Tajima's D")
lines(density(comp.neo.va$TajimaD))
abline(v=quantile(comp.neo$TajimaD,c(0.05)),col='red')

plot(density(comp.jav.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='javanica',xlab="Tajima's D")
lines(density(comp.jav.va$TajimaD))
abline(v=quantile(comp.jav$TajimaD,c(0.05)),col='red')

plot(density(comp.smi.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='smithii',xlab="Tajima's D")
lines(density(comp.smi.va$TajimaD))
abline(v=quantile(comp.smi$TajimaD,c(0.05)),col='red')

plot(density(comp.dim.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='dimidiata',xlab="Tajima's D")
lines(density(comp.dim.va$TajimaD))
abline(v=quantile(comp.dim$TajimaD,c(0.05)),col='red')

## Comparison with null distributions
#plot(density(comp.rus.is$TajimaD),ylim=c(0,2.2),xlim=c(-3,2),col='seagreen',main='rustica',xlab="Tajima's D")
#lines(density(perm.taj.rus$TajimaD))
#abline(v=quantile(comp.rus$TajimaD,c(0.05)),col='red')
#
#plot(density(comp.aet.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='aethiopica',xlab="Tajima's D")
#lines(density(perm.taj.aet$TajimaD))
#abline(v=quantile(comp.aet$TajimaD,c(0.05)),col='red')
#
#plot(density(comp.neo.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='neoxena',xlab="Tajima's D")
#lines(density(perm.taj.neo$TajimaD))
#abline(v=quantile(comp.neo$TajimaD,c(0.05)),col='red')
#
#plot(density(comp.jav.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='javanica',xlab="Tajima's D")
#lines(density(perm.taj.jav$TajimaD))
#abline(v=quantile(comp.jav$TajimaD,c(0.05)),col='red')
#
#plot(density(comp.smi.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='smithii',xlab="Tajima's D")
#lines(density(perm.taj.smi$TajimaD))
#abline(v=quantile(comp.smi$TajimaD,c(0.05)),col='red')
#
#plot(density(comp.dim.is$TajimaD),ylim=c(0,2.2),col='seagreen',main='dimidiata',xlab="Tajima's D")
#lines(density(perm.taj.dim$TajimaD))
#abline(v=quantile(comp.dim$TajimaD,c(0.05)),col='red')

## Export at 8 x 6.5

### Fisher's Exact Tests of enrichment of Fst peaks for low Tajima's D values ----

## rustica
low.is.rus <- as.numeric(length(row.names(comp.rus.is %>% filter(TajimaD <= quantile(comp.rus$TajimaD,c(0.05))))))
low.va.rus <- as.numeric(length(row.names(comp.rus.va %>% filter(TajimaD <= quantile(comp.rus$TajimaD,c(0.05))))))
all.is.rus <- as.numeric(length(rownames(comp.rus.is))) - low.is.rus
all.va.rus <- as.numeric(length(rownames(comp.rus.va))) - low.va.rus

test.rus <- matrix(c(low.is.rus,all.is.rus,low.va.rus,all.va.rus),nrow=2,ncol=2)
fisher.test(test.rus)

## aethiopica
low.is.aet <- as.numeric(length(row.names(comp.aet.is %>% filter(TajimaD <= quantile(comp.aet$TajimaD,c(0.05))))))
low.va.aet <- as.numeric(length(row.names(comp.aet.va %>% filter(TajimaD <= quantile(comp.aet$TajimaD,c(0.05))))))
all.is.aet <- as.numeric(length(rownames(comp.aet.is))) - low.is.aet
all.va.aet <- as.numeric(length(rownames(comp.aet.va))) - low.va.aet

test.aet <- matrix(c(low.is.aet,all.is.aet,low.va.aet,all.va.aet),nrow=2,ncol=2)
fisher.test(test.aet)

## neoxena
low.is.neo <- as.numeric(length(row.names(comp.neo.is %>% filter(TajimaD <= quantile(comp.neo$TajimaD,c(0.05))))))
low.va.neo <- as.numeric(length(row.names(comp.neo.va %>% filter(TajimaD <= quantile(comp.neo$TajimaD,c(0.05))))))
all.is.neo <- as.numeric(length(rownames(comp.neo.is))) - low.is.neo
all.va.neo <- as.numeric(length(rownames(comp.neo.va))) - low.va.neo

test.neo <- matrix(c(low.is.neo,all.is.neo,low.va.neo,all.va.neo),nrow=2,ncol=2)
fisher.test(test.neo)

## javanica
low.is.jav <- as.numeric(length(row.names(comp.jav.is %>% filter(TajimaD <= quantile(comp.jav$TajimaD,c(0.05))))))
low.va.jav <- as.numeric(length(row.names(comp.jav.va %>% filter(TajimaD <= quantile(comp.jav$TajimaD,c(0.05))))))
all.is.jav <- as.numeric(length(rownames(comp.jav.is))) - low.is.jav
all.va.jav <- as.numeric(length(rownames(comp.jav.va))) - low.va.jav

test.jav <- matrix(c(low.is.jav,all.is.jav,low.va.jav,all.va.jav),nrow=2,ncol=2)
fisher.test(test.jav)

## smithii
low.is.smi <- as.numeric(length(row.names(comp.smi.is %>% filter(TajimaD <= quantile(comp.smi$TajimaD,c(0.05))))))
low.va.smi <- as.numeric(length(row.names(comp.smi.va %>% filter(TajimaD <= quantile(comp.smi$TajimaD,c(0.05))))))
all.is.smi <- as.numeric(length(rownames(comp.smi.is))) - low.is.smi
all.va.smi <- as.numeric(length(rownames(comp.smi.va))) - low.va.smi

test.smi <- matrix(c(low.is.smi,all.is.smi,low.va.smi,all.va.smi),nrow=2,ncol=2)
fisher.test(test.smi)

## dimidiata
low.is.dim <- as.numeric(length(row.names(comp.dim.is %>% filter(TajimaD <= quantile(comp.dim$TajimaD,c(0.05))))))
low.va.dim <- as.numeric(length(row.names(comp.dim.va %>% filter(TajimaD <= quantile(comp.dim$TajimaD,c(0.05))))))
all.is.dim <- as.numeric(length(rownames(comp.dim.is))) - low.is.dim
all.va.dim <- as.numeric(length(rownames(comp.dim.va))) - low.va.dim

test.dim <- matrix(c(low.is.dim,all.is.dim,low.va.dim,all.va.dim),nrow=2,ncol=2)
fisher.test(test.dim)

### Get quick proportions of low Tajima's D values in peaks vs backgrounds ----

low.is.rus/as.numeric(length(rownames(comp.rus.is)))
low.va.rus/as.numeric(length(rownames(comp.rus.va)))

low.is.aet/as.numeric(length(rownames(comp.aet.is)))
low.va.aet/as.numeric(length(rownames(comp.aet.va)))

low.is.neo/as.numeric(length(rownames(comp.neo.is)))
low.va.neo/as.numeric(length(rownames(comp.neo.va)))

low.is.jav/as.numeric(length(rownames(comp.jav.is)))
low.va.jav/as.numeric(length(rownames(comp.jav.va)))

low.is.smi/as.numeric(length(rownames(comp.smi.is)))
low.va.smi/as.numeric(length(rownames(comp.smi.va)))

low.is.dim/as.numeric(length(rownames(comp.dim.is)))
low.va.dim/as.numeric(length(rownames(comp.dim.va)))

### Plot Fay & Wu's H distributions----
par(mfrow=c(3,2))
boxplot(comp.rus.va$FAY_WU_H,comp.rus.is$FAY_WU_H,outline=F)
boxplot(comp.aet.va$FAY_WU_H,comp.aet.is$FAY_WU_H,outline=F)
boxplot(comp.neo.va$FAY_WU_H,comp.neo.is$FAY_WU_H,outline=F)
boxplot(comp.jav.va$FAY_WU_H,comp.jav.is$FAY_WU_H,outline=F)
boxplot(comp.smi.va$FAY_WU_H,comp.smi.is$FAY_WU_H,outline=F)
boxplot(comp.dim.va$FAY_WU_H,comp.dim.is$FAY_WU_H,outline=F)

## Comparison with null distributions
#boxplot(perm.fay.rus$FAY_WU_H,comp.rus.is$FAY_WU_H,outline=F)
#boxplot(perm.fay.aet$FAY_WU_H,comp.aet.is$FAY_WU_H,outline=F)
#boxplot(perm.fay.neo$FAY_WU_H,comp.neo.is$FAY_WU_H,outline=F)
#boxplot(perm.fay.jav$FAY_WU_H,comp.jav.is$FAY_WU_H,outline=F)
#boxplot(perm.fay.smi$FAY_WU_H,comp.smi.is$FAY_WU_H,outline=F)
#boxplot(perm.fay.dim$FAY_WU_H,comp.dim.is$FAY_WU_H,outline=F)

plot(density(comp.rus.is$FAY_WU_H),ylim=c(0,3.2),xlim=c(-3,2),col='seagreen',main='rustica',xlab="Fay & Wu's H")
lines(density(comp.rus.va$FAY_WU_H))
abline(v=quantile(comp.rus$FAY_WU_H,c(0.05)),col='red')

plot(density(comp.aet.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='aethiopica',xlab="Fay & Wu's H")
lines(density(comp.aet.va$FAY_WU_H))
abline(v=quantile(comp.aet$FAY_WU_H,c(0.05)),col='red')

plot(density(comp.neo.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='neoxena',xlab="Fay & Wu's H")
lines(density(comp.neo.va$FAY_WU_H))
abline(v=quantile(comp.neo$FAY_WU_H,c(0.05)),col='red')

plot(density(comp.jav.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='javanica',xlab="Fay & Wu's H")
lines(density(comp.jav.va$FAY_WU_H))
abline(v=quantile(comp.jav$FAY_WU_H,c(0.05)),col='red')

plot(density(comp.smi.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='smithii',xlab="Fay & Wu's H")
lines(density(comp.smi.va$FAY_WU_H))
abline(v=quantile(comp.smi$FAY_WU_H,c(0.05)),col='red')

plot(density(comp.dim.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='dimidiata',xlab="Fay & Wu's H")
lines(density(comp.dim.va$FAY_WU_H))
abline(v=quantile(comp.dim$FAY_WU_H,c(0.05)),col='red')

## Comparison with null distributions
#plot(density(comp.rus.is$FAY_WU_H),ylim=c(0,3.2),xlim=c(-3,2),col='seagreen',main='rustica',xlab="Fay & Wu's H")
#lines(density(perm.fay.rus$FAY_WU_H))
#abline(v=quantile(comp.rus$FAY_WU_H,c(0.05)),col='red')
#
#plot(density(comp.aet.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='aethiopica',xlab="Fay & Wu's H")
#lines(density(perm.fay.aet$FAY_WU_H))
#abline(v=quantile(comp.aet$FAY_WU_H,c(0.05)),col='red')
#
#plot(density(comp.neo.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='neoxena',xlab="Fay & Wu's H")
#lines(density(perm.fay.neo$FAY_WU_H))
#abline(v=quantile(comp.neo$FAY_WU_H,c(0.05)),col='red')
#
#plot(density(comp.jav.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='javanica',xlab="Fay & Wu's H")
#lines(density(perm.fay.jav$FAY_WU_H))
#abline(v=quantile(comp.jav$FAY_WU_H,c(0.05)),col='red')
#
#plot(density(comp.smi.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='smithii',xlab="Fay & Wu's H")
#lines(density(perm.fay.smi$FAY_WU_H))
#abline(v=quantile(comp.smi$FAY_WU_H,c(0.05)),col='red')
#
#plot(density(comp.dim.is$FAY_WU_H),ylim=c(0,2.2),col='seagreen',main='dimidiata',xlab="Fay & Wu's H")
#lines(density(perm.fay.dim$FAY_WU_H))
#abline(v=quantile(comp.dim$FAY_WU_H,c(0.05)),col='red')

## Export at 8 x 6.5

### Fisher's Exact Tests of enrichment of Fst peaks for low Fay & Wu's H values ----

## rustica
low.is.rus <- as.numeric(length(row.names(comp.rus.is %>% filter(FAY_WU_H <= quantile(comp.rus$FAY_WU_H,c(0.05))))))
low.va.rus <- as.numeric(length(row.names(comp.rus.va %>% filter(FAY_WU_H <= quantile(comp.rus$FAY_WU_H,c(0.05))))))
all.is.rus <- as.numeric(length(rownames(comp.rus.is))) - low.is.rus
all.va.rus <- as.numeric(length(rownames(comp.rus.va))) - low.va.rus

test.rus <- matrix(c(low.is.rus,all.is.rus,low.va.rus,all.va.rus),nrow=2,ncol=2)
fisher.test(test.rus)

## aethiopica
low.is.aet <- as.numeric(length(row.names(comp.aet.is %>% filter(FAY_WU_H <= quantile(comp.aet$FAY_WU_H,c(0.05))))))
low.va.aet <- as.numeric(length(row.names(comp.aet.va %>% filter(FAY_WU_H <= quantile(comp.aet$FAY_WU_H,c(0.05))))))
all.is.aet <- as.numeric(length(rownames(comp.aet.is))) - low.is.aet
all.va.aet <- as.numeric(length(rownames(comp.aet.va))) - low.va.aet

test.aet <- matrix(c(low.is.aet,all.is.aet,low.va.aet,all.va.aet),nrow=2,ncol=2)
fisher.test(test.aet)

## neoxena
low.is.neo <- as.numeric(length(row.names(comp.neo.is %>% filter(FAY_WU_H <= quantile(comp.neo$FAY_WU_H,c(0.05))))))
low.va.neo <- as.numeric(length(row.names(comp.neo.va %>% filter(FAY_WU_H <= quantile(comp.neo$FAY_WU_H,c(0.05))))))
all.is.neo <- as.numeric(length(rownames(comp.neo.is))) - low.is.neo
all.va.neo <- as.numeric(length(rownames(comp.neo.va))) - low.va.neo

test.neo <- matrix(c(low.is.neo,all.is.neo,low.va.neo,all.va.neo),nrow=2,ncol=2)
fisher.test(test.neo)

## javanica
low.is.jav <- as.numeric(length(row.names(comp.jav.is %>% filter(FAY_WU_H <= quantile(comp.jav$FAY_WU_H,c(0.05))))))
low.va.jav <- as.numeric(length(row.names(comp.jav.va %>% filter(FAY_WU_H <= quantile(comp.jav$FAY_WU_H,c(0.05))))))
all.is.jav <- as.numeric(length(rownames(comp.jav.is))) - low.is.jav
all.va.jav <- as.numeric(length(rownames(comp.jav.va))) - low.va.jav

test.jav <- matrix(c(low.is.jav,all.is.jav,low.va.jav,all.va.jav),nrow=2,ncol=2)
fisher.test(test.jav)

## smithii
low.is.smi <- as.numeric(length(row.names(comp.smi.is %>% filter(FAY_WU_H <= quantile(comp.smi$FAY_WU_H,c(0.05))))))
low.va.smi <- as.numeric(length(row.names(comp.smi.va %>% filter(FAY_WU_H <= quantile(comp.smi$FAY_WU_H,c(0.05))))))
all.is.smi <- as.numeric(length(rownames(comp.smi.is))) - low.is.smi
all.va.smi <- as.numeric(length(rownames(comp.smi.va))) - low.va.smi

test.smi <- matrix(c(low.is.smi,all.is.smi,low.va.smi,all.va.smi),nrow=2,ncol=2)
fisher.test(test.smi)

## dimidiata
low.is.dim <- as.numeric(length(row.names(comp.dim.is %>% filter(FAY_WU_H <= quantile(comp.dim$FAY_WU_H,c(0.05))))))
low.va.dim <- as.numeric(length(row.names(comp.dim.va %>% filter(FAY_WU_H <= quantile(comp.dim$FAY_WU_H,c(0.05))))))
all.is.dim <- as.numeric(length(rownames(comp.dim.is))) - low.is.dim
all.va.dim <- as.numeric(length(rownames(comp.dim.va))) - low.va.dim

test.dim <- matrix(c(low.is.dim,all.is.dim,low.va.dim,all.va.dim),nrow=2,ncol=2)
fisher.test(test.dim)

### Get quick proportions of low Fay & Wu's H values in peaks vs backgrounds ----

low.is.rus/as.numeric(length(rownames(comp.rus.is)))
low.va.rus/as.numeric(length(rownames(comp.rus.va)))

low.is.aet/as.numeric(length(rownames(comp.aet.is)))
low.va.aet/as.numeric(length(rownames(comp.aet.va)))

low.is.neo/as.numeric(length(rownames(comp.neo.is)))
low.va.neo/as.numeric(length(rownames(comp.neo.va)))

low.is.jav/as.numeric(length(rownames(comp.jav.is)))
low.va.jav/as.numeric(length(rownames(comp.jav.va)))

low.is.smi/as.numeric(length(rownames(comp.smi.is)))
low.va.smi/as.numeric(length(rownames(comp.smi.va)))

low.is.dim/as.numeric(length(rownames(comp.dim.is)))
low.va.dim/as.numeric(length(rownames(comp.dim.va)))
