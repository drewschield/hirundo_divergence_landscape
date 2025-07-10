#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Note: this script expects a command line argument in the form of 'chr1', for example.

## Load libraries
library(tidyverse)
library(data.table)
library(R.utils)
library(vcfR)
library(rehh)

## Convert VCF to haplo format
rus <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.rustica.',args[1],'.vcf.gz'), polarize_vcf = TRUE)
aet <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.aethiopica.',args[1],'.vcf.gz'), polarize_vcf = TRUE)
smi <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.smithii.',args[1],'.vcf.gz'), polarize_vcf = TRUE)
neo <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.neoxena.',args[1],'.vcf.gz'), polarize_vcf = TRUE)
jav <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.javanica.',args[1],'.vcf.gz'), polarize_vcf = TRUE)
dim <- data2haplohh(hap_file = paste0('./vcf/hirundo_genus.faywu.snps.polarized.dimidiata.',args[1],'.vcf.gz'), polarize_vcf = TRUE)

## Perform SFS scans
rus.sfs <- calc_sfs_tests(haplohh=rus,window_size=5000)
aet.sfs <- calc_sfs_tests(haplohh=aet,window_size=5000)
smi.sfs <- calc_sfs_tests(haplohh=smi,window_size=5000)
neo.sfs <- calc_sfs_tests(haplohh=neo,window_size=5000)
jav.sfs <- calc_sfs_tests(haplohh=jav,window_size=5000)
dim.sfs <- calc_sfs_tests(haplohh=dim,window_size=5000)

## Write results to files
write.table(rus.sfs,file=paste0('./results/sfs.rustica.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")
write.table(aet.sfs,file=paste0('./results/sfs.aethiopica.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")
write.table(smi.sfs,file=paste0('./results/sfs.smithii.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")
write.table(neo.sfs,file=paste0('./results/sfs.neoxena.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")
write.table(jav.sfs,file=paste0('./results/sfs.javanica.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")
write.table(dim.sfs,file=paste0('./results/sfs.dimidiata.',args[1],'.5kb.txt'),row.names=F,quote=F,sep="\t")

## End Analysis
quit(save="no")
