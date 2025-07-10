# Hirundo genomic divergence landscape analysis
# fd.R - compare fd to Fst genome scans

### Working directory and dependencies----
library(data.table)
library(tidyverse)
library(scales)

setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/fd')

### Read in data----
fd.aet.smi <- read.csv('./rustica_aethiopica_smithii_1Mb_fdfixed.csv',header=T)
fd.smi.dim <- read.csv('./rustica_smithii_dimidiata_1Mb_fdfixed.csv',header=T)

fst <- read.table('../pixy/results/pixy.all.1mb.fst.txt',header=T)
fst.aet.smi <- fst %>% filter(pop1 == 'HA' & pop2 == 'HS')
fst.smi.dim <- fst %>% filter(pop1 == 'HS' & pop2 == 'HD')

# Format and merge
colnames(fd.aet.smi)[1] <- "chromosome"
colnames(fd.aet.smi)[2] <- "window_pos_1"
colnames(fd.aet.smi)[3] <- "window_pos_2"

colnames(fd.smi.dim)[1] <- "chromosome"
colnames(fd.smi.dim)[2] <- "window_pos_1"
colnames(fd.smi.dim)[3] <- "window_pos_2"

fd.aet.smi <- fd.aet.smi %>% filter(D >= 0) # Removing windows where fd cannot be interpreted (i.e., negative D)
fd.smi.dim <- fd.smi.dim %>% filter(D >= 0) # Removing windows where fd cannot be interpreted (i.e., negative D)

comp.aet.smi <- merge(fst.aet.smi, fd.aet.smi, by=c("chromosome","window_pos_1"),sort=F)
comp.smi.dim <- merge(fst.smi.dim, fd.smi.dim, by=c("chromosome","window_pos_1"),sort=F)

plot(fst.aet.smi$avg_wc_fst,fst.smi.dim$avg_wc_fst,pch=20)

### Plot fd and examine relationship between Fst and fd ----
par(mfrow=c(1,2))
plot(comp.aet.smi$avg_wc_fst,comp.aet.smi$fd,pch=20)
cor.test(comp.aet.smi$avg_wc_fst,comp.aet.smi$fd,method='spearman')
plot(comp.smi.dim$avg_wc_fst,comp.smi.dim$fd,pch=20)
cor.test(comp.smi.dim$avg_wc_fst,comp.smi.dim$fd,method='spearman')

### Summarize distributions of fd ----

mean(comp.aet.smi$fd)
sd(comp.aet.smi$fd)

mean(comp.smi.dim$fd)
sd(comp.smi.dim$fd)

