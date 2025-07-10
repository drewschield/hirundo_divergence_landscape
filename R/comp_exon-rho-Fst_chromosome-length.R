# Hirundo genomic divergence landscape analysis
# comp_exon-rho-Fst_chromosome-length.R - compare correlations between exon density and recombination rate and Fst in 50kb windows as a function of chromosome length

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/pixy')

### Read in data----
## Note: opened recombination map in BBedit and changed '.\n' to 'NA\n' to get rid of stupid missing data 'characters' that goof up import
rho <- read.table('../pyrho/results/recombination.rustica.rmap.50kb.txt',header=T,stringsAsFactors=F)
rho$start <- as.integer(rho$start+1)
colnames(rho)[1] <- "chromosome"
#rho$rate[rho$rate == '.'] <- 'NA'

## Read in and process exon density
exon <- read.table('../exon/exon.density.50kb.txt',header=T,stringsAsFactors=F)
exon$start <- exon$start+1
exon$density <- as.integer(exon$bp)/50000
colnames(exon)[1] <- "chromosome"

## Read in an parse Fst for rustica-aethiopica
fst <- read.table('./results/pixy.all.50kb.fst.txt',header=T)
fst <- fst %>% filter((pop1 == 'HRR') & (pop2 == 'HA'))
colnames(fst)[4] <- "start"
colnames(fst)[5] <- "end"

## Find overlap in chromosome and start position and merge data frames
merge(x, y, by=c("k1","k2")) # NA's match

foo <- merge(fst, exon, by=c("chromosome","start"),sort=F)
comp <- merge(foo, rho, by=c("chromosome","start"),sort=F)

### Set chromosome length vector (in logical order by normal ordering scheme)----
chrom.length <- c(119023421,76187387,156035725,116801625,73257097,9617204,63258489,36085389,38459648,31262510,25880253,20272128,21491857,20890524,
                  18810845,16541138,13985943,11194341,12073725,11382101,15277844,7507825,5297670,6778862,7098401,2102120,6843954,5236451,5553549,
                  1648998,1590086,784579,437724,606149,523230,338027,276370,90132487)

### Parse chromosomes----

comp.1 <- comp %>% filter(chromosome == 'NC_053451.1')
comp.1a <- comp %>% filter(chromosome == 'NC_053453.1')
comp.2 <- comp %>% filter(chromosome == 'NC_053450.1')
comp.3 <- comp %>% filter(chromosome == 'NC_053452.1')
comp.4 <- comp %>% filter(chromosome == 'NC_053454.1')
comp.4a <- comp %>% filter(chromosome == 'NC_053470.1')
comp.5 <- comp %>% filter(chromosome == 'NC_053455.1')
comp.6 <- comp %>% filter(chromosome == 'NC_053457.1')
comp.7 <- comp %>% filter(chromosome == 'NC_053456.1')
comp.8 <- comp %>% filter(chromosome == 'NC_053458.1')
comp.9 <- comp %>% filter(chromosome == 'NC_053459.1')
comp.10 <- comp %>% filter(chromosome == 'NC_053462.1')
comp.11 <- comp %>% filter(chromosome == 'NC_053460.1')
comp.12 <- comp %>% filter(chromosome == 'NC_053461.1')
comp.13 <- comp %>% filter(chromosome == 'NC_053463.1')
comp.14 <- comp %>% filter(chromosome == 'NC_053464.1')
comp.15 <- comp %>% filter(chromosome == 'NC_053466.1')
comp.17 <- comp %>% filter(chromosome == 'NC_053469.1')
comp.18 <- comp %>% filter(chromosome == 'NC_053467.1')
comp.19 <- comp %>% filter(chromosome == 'NC_053468.1')
comp.20 <- comp %>% filter(chromosome == 'NC_053465.1')
comp.21 <- comp %>% filter(chromosome == 'NC_053471.1')
comp.22 <- comp %>% filter(chromosome == 'NC_053477.1')
comp.23 <- comp %>% filter(chromosome == 'NC_053474.1')
comp.24 <- comp %>% filter(chromosome == 'NC_053472.1')
comp.25 <- comp %>% filter(chromosome == 'NC_053478.1')
comp.26 <- comp %>% filter(chromosome == 'NC_053473.1')
comp.27 <- comp %>% filter(chromosome == 'NC_053476.1')
comp.28 <- comp %>% filter(chromosome == 'NC_053475.1')
comp.29 <- comp %>% filter(chromosome == 'NC_053479.1')
comp.31 <- comp %>% filter(chromosome == 'NC_053480.1')
comp.32 <- comp %>% filter(chromosome == 'NC_053481.1')
comp.33 <- comp %>% filter(chromosome == 'NC_053482.1')
comp.34 <- comp %>% filter(chromosome == 'NC_053483.1')
comp.35 <- comp %>% filter(chromosome == 'NC_053484.1')
comp.36 <- comp %>% filter(chromosome == 'NC_053485.1')
comp.37 <- comp %>% filter(chromosome == 'NC_053486.1')
comp.Z <- comp %>% filter(chromosome == 'NC_053488.1')

### Calculate correlations and append to vectors----

## Exon density vs recombination rate

cor.exon.rho <- c()
cor.exon.rho <- append(cor.exon.rho,cor(comp.1$density,comp.1$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.1a$density,comp.1a$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.2$density,comp.2$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.3$density,comp.3$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.4$density,comp.4$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.4a$density,comp.4a$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.5$density,comp.5$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.6$density,comp.6$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.7$density,comp.7$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.8$density,comp.8$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.9$density,comp.9$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.10$density,comp.10$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.11$density,comp.11$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.12$density,comp.12$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.13$density,comp.13$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.14$density,comp.14$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.15$density,comp.15$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.17$density,comp.17$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.18$density,comp.18$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.19$density,comp.19$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.20$density,comp.20$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.21$density,comp.21$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.22$density,comp.22$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.23$density,comp.23$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.24$density,comp.24$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.25$density,comp.25$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.26$density,comp.26$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.27$density,comp.27$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.28$density,comp.28$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.29$density,comp.29$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.31$density,comp.31$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.32$density,comp.32$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.33$density,comp.33$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.34$density,comp.34$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.35$density,comp.35$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.36$density,comp.36$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.37$density,comp.37$rate,method='spearman',use='complete.obs'))
cor.exon.rho <- append(cor.exon.rho,cor(comp.Z$density,comp.Z$rate,method='spearman',use='complete.obs'))

for (i in cor.exon.rho){print(i)}


cor.exon.fst <- c()
cor.exon.fst <- append(cor.exon.fst,cor(comp.1$density,comp.1$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.1a$density,comp.1a$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.2$density,comp.2$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.3$density,comp.3$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.4$density,comp.4$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.4a$density,comp.4a$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.5$density,comp.5$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.6$density,comp.6$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.7$density,comp.7$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.8$density,comp.8$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.9$density,comp.9$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.10$density,comp.10$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.11$density,comp.11$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.12$density,comp.12$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.13$density,comp.13$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.14$density,comp.14$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.15$density,comp.15$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.17$density,comp.17$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.18$density,comp.18$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.19$density,comp.19$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.20$density,comp.20$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.21$density,comp.21$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.22$density,comp.22$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.23$density,comp.23$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.24$density,comp.24$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.25$density,comp.25$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.26$density,comp.26$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.27$density,comp.27$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.28$density,comp.28$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.29$density,comp.29$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.31$density,comp.31$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.32$density,comp.32$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.33$density,comp.33$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.34$density,comp.34$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.35$density,comp.35$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.36$density,comp.36$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.37$density,comp.37$avg_wc_fst,method='spearman',use='complete.obs'))
cor.exon.fst <- append(cor.exon.fst,cor(comp.Z$density,comp.Z$avg_wc_fst,method='spearman',use='complete.obs'))

for (i in cor.exon.fst){print(i)}

par(mfrow=c(2,1))
plot(comp$density,type='l')
plot(comp$rate,type='l')

cor.rho.fst <- c()
cor.rho.fst <- append(cor.rho.fst,cor(comp.1$rate,comp.1$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.1a$rate,comp.1a$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.2$rate,comp.2$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.3$rate,comp.3$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.4$rate,comp.4$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.4a$rate,comp.4a$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.5$rate,comp.5$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.6$rate,comp.6$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.7$rate,comp.7$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.8$rate,comp.8$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.9$rate,comp.9$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.10$rate,comp.10$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.11$rate,comp.11$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.12$rate,comp.12$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.13$rate,comp.13$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.14$rate,comp.14$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.15$rate,comp.15$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.17$rate,comp.17$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.18$rate,comp.18$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.19$rate,comp.19$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.20$rate,comp.20$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.21$rate,comp.21$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.22$rate,comp.22$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.23$rate,comp.23$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.24$rate,comp.24$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.25$rate,comp.25$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.26$rate,comp.26$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.27$rate,comp.27$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.28$rate,comp.28$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.29$rate,comp.29$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.31$rate,comp.31$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.32$rate,comp.32$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.33$rate,comp.33$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.34$rate,comp.34$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.35$rate,comp.35$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.36$rate,comp.36$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.37$rate,comp.37$avg_wc_fst,method='spearman',use='complete.obs'))
cor.rho.fst <- append(cor.rho.fst,cor(comp.Z$rate,comp.Z$avg_wc_fst,method='spearman',use='complete.obs'))

for (i in cor.rho.fst){print(i)}

### Make plots of correlations by chromosome length

par(mfrow=c(3,1))
plot(chrom.length,cor.exon.rho,pch=20,col='slategray',xlab='Chromosome Length (bp)',ylab='r Exon Density ~ Rate')
abline(h=0,lty=2,col='grey')
plot(chrom.length,cor.exon.fst,pch=20,col='slategray',xlab='Chromosome Length (bp)',ylab='r Exon Density ~ Fst')
abline(h=0,lty=2,col='grey')
plot(chrom.length,cor.rho.fst,pch=20,col='slategray',xlab='Chromosome Length (bp)',ylab='r Rate ~ Fst')
abline(h=0,lty=2,col='grey')

## output at 8 x 3.5


