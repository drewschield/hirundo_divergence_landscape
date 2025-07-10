# Hirundo Genomic Divergence Landscape

![Genome scans](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps used to study correlated genomic landscapes across a continuum of evoluitonary divergence in swallows (genus _Hirundo_). This workflow is a companion to the methods described in Schield et al. (_in Press_).

__Note__: this repository assumes a specific file organization. Your environment may vary and will need to be adjusted accordingly.

------------------------------------------------------------------------------------------
## Contents

* [Part 0 - General information](#part-0---general-information)
* [Part 1 - Processing, mapping, and variant calling](#part-1---processing-mapping-and-variant-calling)
* [Part 2 - Phylogeny](#part-2---phylogeny)
* [Part 3 - Pixy: calculation of Fst, dxy, and pi](#part-3---pixy-calculation-of-fst-dxy-and-pi)
* [Part 4 - Estimation of historical demography using SMC++](#part-4---estimation-of-historical-demography-using-smc++)


------------------------------------------------------------------------------------------
## Software and dependencies

The analysis sections below use the following software and dependencies and assume they are on the user path unless otherwise specified:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bwa](http://bio-bwa.sourceforge.net/)
* [htslib](http://www.htslib.org/)
* [Samtools](http://www.htslib.org/)
* [bgzip](http://www.htslib.org/)
* [tabix](http://www.htslib.org/)
* [GATK](https://gatk.broadinstitute.org/hc/en-us)
* [Picard](https://broadinstitute.github.io/picard/)
* [bcftools](http://www.htslib.org/)
* [Plink](https://www.cog-genomics.org/plink/)
* [vcftools](https://vcftools.github.io/index.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [ADMIXTURE](https://dalexander.github.io/admixture/publications.html)
* [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
* [SMC++](https://github.com/popgenmethods/smcpp)
* [pyrho](https://github.com/popgenmethods/pyrho)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [VCF-kit](https://vcf-kit.readthedocs.io/en/latest/)
* [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/)
* [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html)
* [MashMap](https://github.com/marbl/MashMap)
* [Repeatmasker](https://www.repeatmasker.org/)
* [R](https://cran.r-project.org/)

[Back to top](#contents)

------------------------------------------------------------------------------------------
## Part 0 - General information

### Overview

This project aims to characterize the evolutionary processes shaping the genomic landscape of divergence among species of Hirundo.

Key analyses will involve the inference of historical demography within/among species, comparative analysis of recombination rate, and investigations of differentiation and sequence divergence among species.

### Sampling

We have whole genome data (generated in 2018 & 2021) for hundreds of barn swallows (_Hirundo rustica_) that were analyzed as part of our study on the role of sexual selection in speciation.

We also have whole genome data (generated in 2018) for Hirundo congeners from the different clades recovered by Dor et al. 2010:
* _Hirundo smithii_ (n = 10; barn swallow clade)
* _Hirundo aethiopica_ (n = 5; barn swallow clade)
* _Hirundo angolensis_ (n = 1; barn swallow clade)
* _Hirundo albigularis_ (n = 1; barn swallow clade)
* _Hirundo nigrita_ (n = 1; barn swallow clade)
* _Hirundo neoxena_ (n = 10; Pacific clade)
* _Hirundo tahitica_ (n = 5; Pacific clade)
* _Hirundo dimidiata_ (n = 8; pearl-breasted clade)
* _Hirundo atrocerulea_ (n = 1; blue clade)

Also included in this dataset is an outgroup sample:
* _Petrochelidon pyrrhonota_ (n = 1; cliff swallow)

### Data locations

The data to analyze for this project can be processed alongside existing directories for the barn swallow speciation genomics project.

Processing directories for fastq, filtered fastq, bam, gvcf, and vcf files are in:

```
/media/drewschield/VernalBucket/hirundo/fastq
/media/drewschield/VernalBucket/hirundo/fastq_filtered
/media/drewschield/VernalBucket/hirundo/bam
/media/drewschield/VernalBucket/hirundo/gvcf
/media/drewschield/VernalBucket/hirundo/vcf
```

These directories contain the files for all ingroup barn swallow samples AND one H. smithii samples (RS_5).

Raw fastq data for the congeners are in `/media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus`.

### General organization

We'll set up a main working directory within the `/data3/` drive.

```
cd /data3
mkdir hirundo_divergence_landscape
cd hirundo_divergence_landscape
```

[Back to top](#contents)

------------------------------------------------------------------------------------------
## Part 1 - Processing, mapping, and variant calling

### Overview

In this part of the workflow, we'll quality trim raw reads, perform mapping to the reference genome, and call variants for the _Hirundo_ congeners.

------------------------------------------------------------------------------------------
### Set up environment

We'll process the data using the existing directories on the `VernalBucket` drive:

```
/media/drewschield/VernalBucket/hirundo/fastq_filtered
/media/drewschield/VernalBucket/hirundo/bam
/media/drewschield/VernalBucket/hirundo/gvcf
/media/drewschield/VernalBucket/hirundo/vcf
```

The raw fastq data for the congeners are in `/media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus/`.

Get into working directory:
```
cd /data3/hirundo_divergence_landscape/
mkdir log
```

------------------------------------------------------------------------------------------
### 1. Quality trimming using Trimmomatic

#### 1. Format sample list

`listTrimmomatic.txt`

Note: RS_5 (H. smithii) has already been processed, so is not included in the sample list.
Note: In preliminary runs, samples RS_1 and RS_4 stopped prematurely due to corruption in the compression of the read 2 file (tested using `gunzip -t $file`). We will run these as single end analyses.

#### 2. Format script to run with `parallel`, calling the sample list

`runTrimmomatic_parallel.sh`:

```
indv=$1
trimmomatic PE -phred33 -threads 8 /media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus/fastq/${indv}_*_1.fq.gz /media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus/fastq/${indv}_*_2.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_1_P.trim.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_1_U.trim.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_2_P.trim.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 > ./log/runTrimmomatic.$indv.log
```

#### 3. Run script with `parallel`

```
parallel --progress --joblog ./log/logfile.trimmomatic -j 6 --workdir . ./runTrimmomatic_parallel.sh :::: listTrimmomatic.txt
```

#### 4. Run single end mode for samples with read 2 issues

```
trimmomatic SE -phred33 -threads 8 /media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus/fastq/RS_1_*_1.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/RS_1_1_P.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 > ./log/runTrimmomatic.RS_1.log
trimmomatic SE -phred33 -threads 8 /media/drewschield/DataStore1/data/hirundo_wgs/hirundo_genus/fastq/RS_4_*_1.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/RS_4_1_P.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 > ./log/runTrimmomatic.RS_4.log
```

#### 5. Remove unpaired reads from output directory to save disk space

```
rm /media/drewschield/VernalBucket/hirundo/fastq_filtered/*U.trim.fq.gz
```

------------------------------------------------------------------------------------------
### 2. Mapping to reference genome

We'll map the trimmed reads to the chromosome-assigned reference genome `/media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta`.

#### 1. Format script to run `bwa` to map and index bam files

`runBWA_parallel.sh`

```
indv=$1
bwa mem -t 8 -R "@RG\tID:$indv\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$indv" /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_1_P.trim.fq.gz /media/drewschield/VernalBucket/hirundo/fastq_filtered/${indv}_2_P.trim.fq.gz | samtools sort -@ 8 -O bam -T $indv.temp -o /media/drewschield/VernalBucket/hirundo/bam/$indv.bam -
samtools index -@ 8 /media/drewschield/VernalBucket/hirundo/bam/$indv.bam
```

#### 2. Run script with `parallel`

```
parallel --progress --joblog ./log/logfile.bwa -j 6 --workdir . ./runBWA_parallel.sh :::: listTrimmomatic.txt
```

#### 3. Run separate analyses for RS_1 and RS_4

```
bwa mem -t 8 -R "@RG\tID:RS_1\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:RS_1" /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta /media/drewschield/VernalBucket/hirundo/fastq_filtered/RS_1_1_P.trim.fq.gz | samtools sort -@ 8 -O bam -T RS_1.temp -o /media/drewschield/VernalBucket/hirundo/bam/RS_1.bam -
samtools index -@ 8 /media/drewschield/VernalBucket/hirundo/bam/RS_1.bam
bwa mem -t 8 -R "@RG\tID:RS_4\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:RS_4" /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta /media/drewschield/VernalBucket/hirundo/fastq_filtered/RS_4_1_P.trim.fq.gz | samtools sort -@ 8 -O bam -T RS_4.temp -o /media/drewschield/VernalBucket/hirundo/bam/RS_4.bam -
samtools index -@ 8 /media/drewschield/VernalBucket/hirundo/bam/RS_4.bam
```

------------------------------------------------------------------------------------------
### 3. Variant calling

We'll call genomic variants using `gatk`, starting with individual variant calls using `HaplotypeCaller`, followed by cohort variant calls using `GenotypeGVCFs`.

#### Format sample list

`listGATK.txt`

#### 1. Individual variant calls using HaplotypeCaller

##### 1. Retrieve `GATK` executables to working directory

```
cp -r /data3/hirundo/gatk* .
```

##### 2. Format script to run GATK HaplotypeCaller

`runGATKHaplotypeCaller_parallel.sh`

```
indv=$1
./gatk-4.0.8.1/gatk HaplotypeCaller -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta --ERC GVCF -I /media/drewschield/VernalBucket/hirundo/bam/$indv.bam -O /media/drewschield/VernalBucket/hirundo/gvcf/$indv.raw.snps.indels.g.vcf > ./log/runGATKHaplotypeCaller.$indv.log
bgzip /media/drewschield/VernalBucket/hirundo/gvcf/$indv.raw.snps.indels.g.vcf
```

##### 3. Run script with `parallel`

```
parallel --progress --joblog ./log/logfile.GATKHaplotypeCaller -j 12 --workdir . ./runGATKHaplotypeCaller_parallel.sh :::: listGATK.txt
```

##### 4. Tabix index output gVCFs

```
for i in `cat listGATK.txt`; do tabix -p vcf /media/drewschield/VernalBucket/hirundo/gvcf/$i.raw.snps.indels.g.vcf.gz; done
```

#### 2. Cohort variant calls using GenotypeGVCFs

#### Format sample list (with paths to gVCF files)

This includes barn swallow samples and the H. smithii sample used as an outgroup in previous analyses (RS_5), plus the newly generated gVCFs for the other congeners.

`listGVCF.list`

Note: this list does not include samples in the `sample.remove.list` from the previous Hirundo speciation genomics study. This list included samples with high missing data and/or that were mislabeled during the library preparation step.

#### Format scaffold interval lists

These will be used to run parallel cohort variant calling on subsets of the genome, to then be merged later.

```
./log/scaffold.list1.intervals
./log/scaffold.list2.intervals
./log/scaffold.list3.intervals
./log/scaffold.list4.intervals
./log/scaffold.list5.intervals
./log/scaffold.list6.intervals
./log/scaffold.list7.intervals
./log/scaffold.list8.intervals
./log/scaffold.list9.intervals
./log/scaffold.list10.intervals
./log/scaffold.list11.intervals
./log/scaffold.list12.intervals
./log/scaffold.list13.intervals
```

#### Run GenotypeGVCFs on sets of intervals

```
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list1.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list1.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list1.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list2.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list2.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list2.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list3.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list3.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list3.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list4.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list4.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list4.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list5.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list5.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list5.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list6.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list6.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list6.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list7.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list7.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list7.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list8.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list8.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list8.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list9.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list9.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list9.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list10.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list10.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list10.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list11.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list11.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list11.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list12.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list12.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list12.vcf.log &
nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list13.intervals -V ./listGVCF.list -allSites -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list13.vcf.gz > ./log/GenotypeGVCFs.hirundo_genus.allsites.raw.scaffold.list13.vcf.log &
```

------------------------------------------------------------------------------------------
### 4. Sex identification

We'll use relative read depths on the Z chromosome and autosomes to infer genetic sex of each individual.

We'll extract mapping data for Chromosome 3, taking the median as an 'autosomal median', then compare levels of coverage for individuals to this value across the Z chromosome scaffold.

We already have these data for the barn swallow samples, so here we'll focus on the congeners.

#### Set up environment

```
cd /data3/hirundo_divergence_landscape
mkdir analysis
cd analysis
mkdir sex_identification
cd sex_identification
mkdir mosdepth_results
```

#### 1. Format sample and scaffold input files

Make a list of the 48 congener/outgroup samples: `sample.list`.

Make 'genome' files for chromosome 3 and the Z chromosome, with scaffold ID and length:

```
chrom.chr3.genome
chrom.chrZ.genome
```

Converted genome files to BED sliding window files, in 10 kb windows:

```
bedtools makewindows -g chrom.chr3.genome -w 10000 > chrom.chr3.10kb.bed
bedtools makewindows -g chrom.chrZ.genome -w 10000 > chrom.chrZ.10kb.bed
```

#### 2. Write Python script to quantify autosome median and log2 Z/Autosome depth ratio

`identifySex.py`

```
import sys
from numpy import array
from numpy import median
from numpy import mean
from numpy import log2

auto = []
za_norm = []

for line in open(sys.argv[1], 'r'):
	acov = float(line.split()[3])
	auto.append(acov)

auto_med = median(array(auto))

for line in open(sys.argv[2], 'r'):
	zcov = line.split()[3]
	if float(zcov) > 0.00:
		za = log2(float(zcov)/float(auto_med))
		za_norm.append(za)

za_mean = mean(array(za_norm))

sample = str(sys.argv[1])
sample = sample.split('/')[1]
sample = sample.split('.')[0]

if za_mean < -0.2:
	message = 'is likely female'
	print za_mean, sample, message
else:
	message = 'is likely male'
	print za_mean, sample, message
```

#### 3. Write wrapper script to run mosdepth and identifySex.py

`runSexIdentification.sh`

```
list=$1
for i in `cat $list`; do
	mosdepth -t 4 --fast-mode -n -b chrom.chr3.10kb.bed -f /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta ./mosdepth_results/$i.chr3 /media/drewschield/VernalBucket/hirundo/bam/$i.bam
	mosdepth -t 4 --fast-mode -n -b chrom.chrZ.10kb.bed -f /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta ./mosdepth_results/$i.chrZ /media/drewschield/VernalBucket/hirundo/bam/$i.bam
	gunzip ./mosdepth_results/$i.chr3.regions.bed.gz 
	gunzip ./mosdepth_results/$i.chrZ.regions.bed.gz
	python identifySex.py mosdepth_results/$i.chr3.regions.bed mosdepth_results/$i.chrZ.regions.bed
done
```

#### 4. Run script

```
nohup sh runSexIdentification.sh sample.list > runSexIdentification.log &
```

#### 5. Extract general coverage information from mosdepth results

```
for indv in `cat sample.list`; do mean=`awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ./mosdepth_results/$indv.chr3.regions.bed`; echo $indv $mean; done
```

------------------------------------------------------------------------------------------
### 5. Variant filtering

In this section we'll impose various filtering steps to produce high-quality SNPs for downstream analysis.

#### 1. Hard filters in GATK

We'll first flag genotypes to be filtered based on the following specs:
* QD < 2.0
* FS > 60.0
* MQ < 40.0
* MQRankSum < -12.5
* ReadPosRankSum < -8.0

Then we'll set the annotated genotypes as missing in the all-sites VCFs.

Run GATK VariantFiltration to annotate genotypes that don't pass the hard filters.
```
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list1.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list1.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list1.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list2.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list2.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list3.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list3.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list3.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list4.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list4.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list4.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list5.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list5.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list5.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list6.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list6.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list6.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list7.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list7.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list7.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list8.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list8.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list8.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list9.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list9.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list9.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list10.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list10.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list10.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list11.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list11.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list11.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list12.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list12.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list12.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list13.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.scaffold.list13.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list13.log &
```

Format `./listVCFHardFilter.txt` with paths to interval HardFilter VCFs.

Run Picard to merge VCFs.
```
nohup java -jar picard.jar MergeVcfs -I ./listVCFHardFilter.txt -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.vcf.gz > ./log/Picard_MergeVCFs.HardFilter.log &
```

Remove raw interval VCFs taking up disk space.
```
rm /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.raw.scaffold.list*
```

Run bcftools to mask indels/hard filtered sites.
```
bcftools filter --threads 36 -e 'TYPE="indel" || FILTER="QD2" || FILTER="FS60" || FILTER="MQ40" || FILTER="MQRankSum-12.5" || FILTER="ReadPosRankSum-8"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.vcf.gz
bcftools index --threads 48 -t /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.vcf.gz
```

#### 2. Filtering female heterozygous sites on the sex chromosomes

Females are hemizygous ZW, so should not have heterozygous genotype calls on either the Z or W chromosome.

We'll identify any heterozygous calls in known females then conservatively mask these in all individuals.

##### 1. Format scaffold-assigned BED files for parsing chromosomes

```
./log/scaffold.assigned-autosome.bed
./log/scaffold.assigned-chrZ.bed
./log/scaffold.assigned-chrW.bed
```

##### 2. Format female individual list

This is based on the sex identification procedure above.

```
./log/sample.female.list
```

##### 3. Parse autosomes, Z chromosome, W chromosome

For the W chromosome, we'll pare down to females, only.

```
bcftools view --threads 16 -R ./log/scaffold.assigned-autosome.bed -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.auto.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.vcf.gz
bcftools view --threads 16 -R ./log/scaffold.assigned-chrZ.bed -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.vcf.gz
bcftools view --threads 16 -S ./log/sample.female.list -R ./log/scaffold.assigned-chrW.bed -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.vcf.gz
```

##### 4. Extract biallelic SNPs from the Z and W chromosome VCFs

```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.vcf.gz
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.vcf.gz
```

##### 5. Identify female heterozygous sites using Python script

Run `./sexChrFemaleHeterozygous.py` to extract female heterozygous genotype positions on the sex chromosomes.
```
gunzip /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.vcf.gz
python sexChrFemaleHeterozygous.py ./log/sample.female.list /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.vcf ./log/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.txt
gunzip /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.snps.vcf.gz
python sexChrFemaleHeterozygous.py ./log/sample.female.list /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.snps.vcf ./log/hirundo_genus.allsites.HardFilter.recode.chrW.snps.FemaleZhetSites.txt
```

##### 6. Convert output to BED format and index with GATK

```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.txt > ./log/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_genus.allsites.HardFilter.recode.chrW.snps.FemaleZhetSites.txt > ./log/hirundo_genus.allsites.HardFilter.recode.chrW.snps.FemaleZhetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_genus.allsites.HardFilter.recode.chrW.snps.FemaleZhetSites.bed
```

##### 7. Run GATK VariantFiltration to annotate female heterozygous sites and mask with bcftools

Z chromosome:
```
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.vcf.gz
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.vcf.gz --mask ./log/hirundo_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed --mask-name ZHET -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_genus.allsites.HardFilter.recode.chrZ.filter.log &
bcftools filter --threads 16 -e 'FILTER="ZHET"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrZ.filter.vcf.gz
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz
```

W chromosome:
```
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.vcf.gz
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.vcf.gz --mask ./log/hirundo_genus.allsites.HardFilter.recode.chrW.snps.FemaleZhetSites.bed --mask-name WHET -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_genus.allsites.HardFilter.recode.chrW.filter.log &
bcftools filter --threads 16 -e 'FILTER="WHET"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.chrW.filter.vcf.gz
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz
```

#### 3. Repeat masking

Repeat annotations are in `/data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed`.

##### 1. Annotate repeats with GATK VariantFiltration

```
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.auto.vcf.gz
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.HardFilter.recode.auto.vcf.gz --mask /data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_genus.auto.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz --mask /data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_genus.chrZ.log &
nohup ./gatk-4.0.8.1/gatk VariantFiltration -V /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz --mask /data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_genus.chrW.log &
```

##### 2. Recode repeat annotations as missing genotypes and index VCFs

```
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.tmp-rep.vcf.gz
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.vcf.gz
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.tmp-rep.vcf.gz
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.tmp-rep.vcf.gz
tabix -p vcf /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz 
```

#### 4. Extract chromosome-specific all-sites VCFs

##### 1. Set up environment

```
cd /media/drewschield/VernalBucket/hirundo/vcf/
mkdir chrom-specific-genus
```

##### 2. Format chromosome-scaffold lookup table

Note: this does not include assigned, but unordered scaffolds.

`./log/chromosome-scaffold.table.txt`

First couple of lines:
```
NC_053451.1	chr1
NC_053453.1	chr1A
NC_053450.1	chr2
NC_053452.1	chr3
NC_053454.1	chr4
```

Also format autosome-specific file (sex chromosomes are last two lines):
```
head -n -2 ./log/chromosome-scaffold.table.txt > ./log/chromosome-scaffold.table.auto.txt
```

##### 3. Parse chromosome-specific autosome VCFs

```
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; echo parsing $chrom; bcftools view --threads 24 -r $scaff -O z -o /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.${chrom}.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.vcf.gz; done < ./log/chromosome-scaffold.table.auto.txt
```

##### 4. Parse Z and W chromosome

```
bcftools view --threads 16 -r NC_053488.1 -O z -o /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.chrZ.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz
bcftools view --threads 16 -r NC_053487.1 -O z -o /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.chrW.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz
```

##### 5. Index VCFs

```
for i in /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/*.vcf.gz; do tabix -p vcf $i; done
```

#### 5. Extract biallelic SNPs

##### 1. Extract SNPs

```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.vcf.gz
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.vcf.gz
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrW.vcf.gz
```

##### 2. Concatenate autosome and Z chromosome SNP VCFs

```
bcftools concat -O z -o /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto+chrZ.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.snps.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.snps.vcf.gz
```

[Back to top](#contents)

------------------------------------------------------------------------------------------
## Part 2 - Phylogeny

### Overview

We'll use genome-wide SNPs to estimate the phylogeny of the genus. This will frame our understanding for downstream analyses.

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundo_divergence_landscape/analysis
mkdir phylogeny
cd phylogeny
mkdir vcf
mkdir raxml
mkdir svdq
```

------------------------------------------------------------------------------------------
### 1. Format input data

We'll extract SNPs with complete data (i.e., 0% missing genotypes) and reformat in phylip and nexus formats for analysis.

#### 1. Format focal sample list

We'll select a representative of each species (and subspecies in the case of H. rustica), along with an outgroup Petrochelidon pyrrhonota.

`./sample.list`

Also, format a table with a second column that includes taxon info:

`./sample.info`

#### 2. Extract SNPs using vcftools

We'll filter so that we retain SNPs with:
* no missing genotypes
* thinned by 10kb on autosomes

```
bcftools view --threads 16 -S sample.list -m2 -M2 -U -v snps -i 'F_MISSING=0.0' -O z -o ./hirundo_genus.auto-tmp.snps.focal.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.snps.vcf.gz
vcftools --gzvcf ./hirundo_genus.auto-tmp.snps.focal.vcf.gz --thin 10000 --recode --stdout | bgzip -c > ./hirundo_genus.auto.snps.focal.vcf.gz
bcftools view --threads 16 -S sample.list -m2 -M2 -U -v snps -i 'F_MISSING=0.0' -O z -o ./hirundo_genus.chrZ.snps.focal.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.snps.vcf.gz
bcftools concat -O z -o ./hirundo_genus.snps.focal.vcf.gz ./hirundo_genus.auto.snps.focal.vcf.gz ./hirundo_genus.chrZ.snps.focal.vcf.gz
```

Clean up:
```
rm hirundo_genus.auto-tmp.snps.focal.vcf.gz
rm hirundo_genus.auto.snps.focal.vcf.gz
rm hirundo_genus.chrZ.snps.focal.vcf.gz
```

#### 3. Convert VCF to alignment formats

We'll format a phylip input file for RAxML and a nexus for SVDquartets.

We'll use a couple of python scripts to do this. Retrieve the scripts:
```
curl -LO https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
curl -LO https://raw.githubusercontent.com/btmartin721/raxml_ascbias/master/ascbias.py
```

Format phylip and nexus input:
```
python vcf2phylip.py -n -i ./hirundo_genus.snps.focal.vcf.gz
```

Get rid of invariant sites:
```
python ascbias.py -p ./hirundo_genus.snps.focal.min4.phy
mv out.phy hirundo.snps.focal.phy
```

------------------------------------------------------------------------------------------
### 2. Perform phylogenetic analysis

We'll perform a concatenated analysis using RAxML and a coalescent-based analysis using SVDquartets.

#### 1. Concatenated analysis

Because we are analyzing only variant sites (including heterozygous sites), we'll specify a substitution model for genotypes that accounts for SNP ascertainment bias.

See this post for some useful info on the GTGTR4+G+ASC_LEWIS model: https://groups.google.com/g/raxml/c/WyelRC1X1Sg

```
cd ./raxml
raxml-ng --threads 12 --all --msa ../hirundo.snps.focal.phy --model GTGTR4+G+ASC_LEWIS --tree pars{10} --bs-trees 100
mv *.raxml.* ./raxml
cd ..
```

#### 2. Coalescent-based analysis

(View nexus input to verify which sample is the outgroup).

```
cd ./svdq
paup
exe ../hirundo_genus.snps.focal.min4.nexus;
outgroup 16;
set outroot=mono;
svdq;
svdq showScores=no seed=19856136 bootstrap nreps=100 treeFile=hirundo_genus.bootstrap.tre;
savetrees file=hirundo_genus.svdq.tre savebootp=nodelabels;
cd ..
```

[Back to top](#contents)

------------------------------------------------------------------------------------------
## Part 3 - Pixy: calculation of Fst, dxy, and pi

### Overview

We'll calculate Fst and dxy between populations, and pi within populations, using Pixy, with chromosome-specific all-sites VCFs as input data.

The input data are in `/media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/`.

------------------------------------------------------------------------------------------
### 1. Set up environment

#### Make analysis directory

```
cd /media/mother/extradrive3/hirundo_divergence_landscape/analysis
mkdir pixy
cd pixy
mkdir results
mkdir log
```

#### Format popmap for analysis

`popmap.pixy`

First few lines:
```
HR1601	HRS
HR1604	HRS
HR1607	HRS
HR1613	HRS
HR1615	HRS
```

The abbreviations correspond to:

```
HRS = Hirundo rustica savignii
HRT = Hirundo rustica tytleri
HRR = Hirundo rustica rustica
HA = Hirundo aethiopica
HS = Hirundo smithii
HN = Hirundo neoxena
HT = Hirundo tahitica
HD = Hirundo dimidiata
```

#### Format chromosome list

`./chrom.list`

------------------------------------------------------------------------------------------
### 2. Format sliding window files

#### 1. First, format 'genome' file with scaffold lengths

`/media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome.genome`

First few lines:

```
NC_053451.1	119023421
NC_053453.1	76187387
NC_053450.1	156035725
NC_053452.1	116801625
NC_053454.1	73257097
```

#### 2. Write sliding window BED files from genome file

```
cd /media/mother/extradrive3/hirundo_divergence_landscape/log/
bedtools makewindows -g chromosome.genome -w 1000000 -s 100000 > chromosome.window.1mb-100kb.bed
bedtools makewindows -g chromosome.genome -w 100000 -s 10000 > chromosome.window.100kb-10kb.bed
bedtools makewindows -g chromosome.genome -w 50000 -s 5000 > chromosome.window.50kb-5kb.bed
```

#### 3. Extract chromosome-specific sliding window files for analyses

```
cd /media/mother/extradrive3/hirundo_divergence_landscape/analysis/pixy/
mkdir windows
head -n -1 /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome-scaffold.table.txt > /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome-scaffold.table.auto+chrZ.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; grep $scaff /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome.window.1mb-100kb.bed > ./windows/window.1mb-100kb.$chrom.bed; grep $scaff /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome.window.100kb-10kb.bed > ./windows/window.100kb-10kb.$chrom.bed; grep $scaff /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome.window.50kb-5kb.bed > ./windows/window.50kb-5kb.$chrom.bed; done < /media/mother/extradrive3/hirundo_divergence_landscape/log/chromosome-scaffold.table.auto+chrZ.txt
```

------------------------------------------------------------------------------------------
### 3. Perform analysis in Pixy

#### 1. Format script to run analyses in sliding windows on each chromosome

runPixy.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.1mb-100kb.$chrom.bed --output_folder results --output_prefix pixy.1mb-100kb.$chrom
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.100kb-10kb.$chrom.bed --output_folder results --output_prefix pixy.100kb-10kb.$chrom
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.50kb-5kb.$chrom.bed --output_folder results --output_prefix pixy.50kb-5kb.$chrom
done
```

#### 2. Run Pixy

```
conda activate pixy
nohup sh runPixy.sh chrom.list > ./log/runPixy.log &
```

#### 3. Concatenate sliding window results across chromosomes

We just want the non-overlapping windows for this for 

awk 'NR % 10 == 0' input > output

```
head -n 1 ./results/pixy.50kb-5kb.chr1_fst.txt > pixy.all.50kb.fst.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.50kb-5kb.${chrom}_fst.txt | awk 'NR % 10 == 0' >> pixy.all.50kb.fst.txt; done
head -n 1 ./results/pixy.50kb-5kb.chr1_dxy.txt > pixy.all.50kb.dxy.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.50kb-5kb.${chrom}_dxy.txt | awk 'NR % 10 == 0' >> pixy.all.50kb.dxy.txt; done
head -n 1 ./results/pixy.50kb-5kb.chr1_pi.txt > pixy.all.50kb.pi.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.50kb-5kb.${chrom}_pi.txt >> pixy.all.50kb.pi.txt; done
```

Concatenate the 1Mb - 100kb windows genome-wide:

```
head -n 1 ./results/pixy.1mb-100kb.chr1_fst.txt > pixy.all.1mb-100kb.fst.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.1mb-100kb.${chrom}_fst.txt >> pixy.all.1mb-100kb.fst.txt; done
head -n 1 ./results/pixy.1mb-100kb.chr1_dxy.txt > pixy.all.1mb-100kb.dxy.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.1mb-100kb.${chrom}_dxy.txt >> pixy.all.1mb-100kb.dxy.txt; done
head -n 1 ./results/pixy.1mb-100kb.chr1_pi.txt > pixy.all.1mb-100kb.pi.txt; for chrom in `cat scaffold.order.list`; do tail -n+2 ./results/pixy.1mb-100kb.${chrom}_pi.txt >> pixy.all.1mb-100kb.pi.txt; done
```

------------------------------------------------------------------------------------------
### 4. Perform analysis in Pixy (inclusive of n = 1 congeners)

To broadly characterize between-species differentiation and genetic diversity, we'll also perform analysis on non-overlapping 1Mb windows for a popmap that includes H. atrocerulea, H. angolensis, H. nigrita, and H. albigularis.

#### Set up environment

```
mkdir results_all
```

#### 1. Format popmap for inclusive set for analysis

`popmap.all.pixy`

First few lines:
```
HR1601	HRS
HR1604	HRS
HR1607	HRS
HR1613	HRS
HR1615	HRS
```

The abbreviations correspond to:

```
HRS = Hirundo rustica savignii
HRT = Hirundo rustica tytleri
HRR = Hirundo rustica rustica
HA = Hirundo aethiopica
HAN = Hirundo angolensis
HNI = Hirundo nigrita
HS = Hirundo smithii
HAL = Hirundo albigularis
HN = Hirundo neoxena
HT = Hirundo tahitica
HD = Hirundo dimidiata
HAT = Hirundo atrocerulea
```

#### 2. Format script to run analyses in sliding windows on each chromosome

runPixy.all.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats pi dxy fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.all.pixy --window 1000000 --output_folder results_all --output_prefix pixy.all.1mb.$chrom
done
```

#### 3. Run Pixy

```
conda activate pixy
nohup sh runPixy.all.sh chrom.list > ./log/runPixy.all.log &
```

#### 4. Concatenated results for non-overlapping 1Mb windows

##### Format ordered chromosome scaffold list

`scaffold.order.list`

First few lines:
```
chr1
chr1A
chr2
chr3
chr4
chr4A
```

##### Write script to concatenate results

concatenatePixyOrder.sh

```
echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.1mb.pi.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.1mb.${chrom}_pi.txt | tail -n +2 >> pixy.all.1mb.pi.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy no_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.1mb.dxy.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.1mb.${chrom}_dxy.txt | tail -n +2 >> pixy.all.1mb.dxy.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.all.1mb.fst.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.1mb.${chrom}_fst.txt | tail -n +2 >> pixy.all.1mb.fst.txt
done
```

`sh concatenatePixyOrder.sh`

------------------------------------------------------------------------------------------
### 4(b). Perform analysis in Pixy (inclusive of n = 1 congeners) in 50kb windows

#### 1. Format script to run analyses in sliding windows on each chromosome

runPixy.all.50kb.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats pi dxy fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.all.pixy --window 50000 --output_folder results_all --output_prefix pixy.all.50kb.$chrom
done
```

#### 2. Run Pixy

```
conda activate pixy
nohup sh runPixy.all.50kb.sh chrom.list > ./log/runPixy.all.50kb.log &
```

#### 3. Concatenate results for non-overlapping 1Mb windows

concatenatePixyOrder.50kb.sh

```
echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.50kb.pi.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.50kb.${chrom}_pi.txt | tail -n +2 >> pixy.all.50kb.pi.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy no_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.50kb.dxy.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.50kb.${chrom}_dxy.txt | tail -n +2 >> pixy.all.50kb.dxy.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.all.50kb.fst.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.50kb.${chrom}_fst.txt | tail -n +2 >> pixy.all.50kb.fst.txt
done
```

`sh concatenatePixyOrder.50kb.sh`

------------------------------------------------------------------------------------------
### 4(c). Perform analysis in Pixy (inclusive of n = 1 congeners) in 5kb windows

#### 1. Format script to run analyses in sliding windows on each chromosome

Note: currently setting this up to only calculate Fst (for permutations to define outlier 50kb windows from random sets of 5kb windows).

runPixy.all.5kb.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.all.pixy --window 5000 --output_folder results_all --output_prefix pixy.all.5kb.$chrom
done
```

#### 2. Run Pixy

```
conda activate pixy
nohup sh runPixy.all.5kb.sh chrom.list > ./log/runPixy.all.5kb.log &
```

#### 3. Concatenate results for non-overlapping 5kb windows

concatenatePixyOrder.5kb.sh

```
echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.all.5kb.fst.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.5kb.${chrom}_fst.txt | tail -n +2 >> pixy.all.5kb.fst.txt
done
```

`sh concatenatePixyOrder.5kb.sh`


------------------------------------------------------------------------------------------
### 5. Perform analysis of barn swallow subspecies in Pixy

To broadly characterize differentiation and genetic diversity within barn swallow, we'll also perform analysis on non-overlapping 1Mb windows for a popmap that includes the six subspecies

#### Set up environment

```
mkdir results_subspecies
```

#### 1. Format popmap for inclusive set for analysis

`popmap.rustica.pixy`

First few lines:
```
HR2600_27318	HRE
HR2600_27714	HRE
HR2600_54936	HRE
HR2640_97014	HRE
HR2640_97129	HRE
```

The abbreviations correspond to:

```
HRS = Hirundo rustica savignii
HRT = Hirundo rustica tytleri
HRR = Hirundo rustica rustica
HRV = Hirundo rustica transitiva
HRG = Hirundo rustica gutturalis
HRE = Hirundo rustica erythrogaster
```

#### 2. Format script to run analyses in sliding windows on each chromosome

runPixy.rustica.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats pi dxy fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.rustica.pixy --window 1000000 --output_folder results_rustica --output_prefix pixy.rustica.1mb.$chrom
done
```

#### 3. Run Pixy

```
conda activate pixy
nohup sh runPixy.rustica.sh chrom.list > ./log/runPixy.rustica.log &
```

#### 4. Concatenated results for non-overlapping 1Mb windows

##### Write script to concatenate results

concatenatePixyOrderSubspecies.sh:

```
echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.rustica.1mb.pi.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_rustica/pixy.rustica.1mb.${chrom}_pi.txt | tail -n +2 >> pixy.rustica.1mb.pi.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy no_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.rustica.1mb.dxy.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_rustica/pixy.rustica.1mb.${chrom}_dxy.txt | tail -n +2 >> pixy.rustica.1mb.dxy.txt
done

echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.rustica.1mb.fst.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_rustica/pixy.rustica.1mb.${chrom}_fst.txt | tail -n +2 >> pixy.rustica.1mb.fst.txt
done
```

`sh concatenatePixyOrderSubspecies.sh`

[Back to top](#contents)

------------------------------------------------------------------------------------------
## Part 4 - Estimation of historical demography using SMC++

### Overview

We'll infer population histories using the sequential Markovian coalescent model implemented in SMC++.

We'll focus all analysis on autosomes to not deal with sex-linkage, except for Z chromosome-specific analyses to test whether inference of Z-linked demography influences downstream recombination rate estimates.

We'll use the following assumptions to scale inferences of effective population size:
* Generation time = 1 year
* Per generation generation rate = 4.6e-9 (from Smeds et al. 2016 germline mutation rate estimation paper)

To enable comparison to other Hirundo species, we'll sample H. rustica rustica from Karasuk, Russia as a representative population.

The input data are in `/media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.snps.vcf.gz`.

------------------------------------------------------------------------------------------
### Install SMC++

The SMC++ documentation can be found at:

https://github.com/popgenmethods/smcpp#installation-instructions

SMC++ relies on several dependencies that are somewhat specific. We'll install everything within a virtual environment.

#### Set up environment

```
cd ~/tmp
mkdir smc++-install
cd smc++-install
```

#### Install Python3.8 & libraries (important for some specific dependencies)

```
sudo apt install python3.8
sudo apt-get install python3.8-dev # Collects 3.8 libraries
```

#### Create and activate new virtual environment for SMC++:

```
virtualenv -p /usr/bin/python3.8 smc++-env
source smc++-env/bin/activate
```

(To get out of virtual environment, type and enter `deactivate`).

#### Install library requirements for SMC++:

```
sudo apt-get install -y python3-dev libgmp-dev libmpfr-dev libgsl0-dev
```

#### Install SMC++

```
pip install git+https://github.com/popgenmethods/smcpp
```

#### Test that build worked

```
smc++ vcf2smc -h
```

------------------------------------------------------------------------------------------
### 1. Set up environment

#### Make analysis directory

```
cd /data3/hirundo_divergence_landscape/analysis
mkdir smc++
cd smc++
mkdir out
mkdir out-chrZ
mkdir analysis
mkdir vcf
mkdir log
```

#### Format popmaps

Individual species:
```
./popmap.rustica
./popmap.aethiopica
./popmap.smithii
./popmap.neoxena
./popmap.tahitica
./popmap.dimidiata
```

All together:
```
cat popmap.* > sample.list
```

All together; males only: `sample.male.list`

#### Format distinguished individual lists

We'll sample 5 individuals per species for autosomes:

```
for i in popmap.*; do pop=`echo $i | cut -d'.' -f2`; shuf -n 5 $i > distinguished.$pop; done
```

#### Format scaffold-chromosome table

`./chromosome-scaffold.table.auto.txt`

------------------------------------------------------------------------------------------
### 2. Extract and format SNP data

#### 1. Extract SNPs for autosomes

```
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; bcftools view --threads 16 -r $scaff -S sample.list -O z -o ./vcf/hirundo_genus.snps.$chrom.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto.snps.vcf.gz; done < ./chromosome-scaffold.table.auto.txt
```

#### 2. Extract SNPs for the Z chromosome

```
bcftools view --threads 16 -r NC_053488.1 -S sample.male.list -O z -o ./vcf/hirundo_genus.snps.chrZ.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.chrZ.snps.vcf.gz
```

#### 3. Index VCFs

```
for vcf in ./vcf/*.vcf.gz; do tabix -C -p vcf $vcf; done
```

------------------------------------------------------------------------------------------
### 3. Convert VCF to SMC input

#### 1. Run conversions for single population inference

We'll generate SMC++ input files for each species and chromosome, iterating to specify each of the 5 distinguished individuals per species.

```
source ~/tmp/smc++-install/smc++-env/bin/activate
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.rustica`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/rustica.$chrom.$indv.smc.gz $scaff rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.aethiopica`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/aethiopica.$chrom.$indv.smc.gz $scaff aethiopica:RS_34,RS_35,RS_36,RS_37,RS_38 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.smithii`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/smithii.$chrom.$indv.smc.gz $scaff smithii:RS_5,RS_1,RS_2,RS_3,RS_4,RS_6,RS_7,RS_8,RS_9,RS_10 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.neoxena`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/neoxena.$chrom.$indv.smc.gz $scaff neoxena:RS_11,RS_12,RS_13,RS_14,RS_15,RS_16,RS_17,RS_18,RS_19,RS_20 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.tahitica`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/tahitica.$chrom.$indv.smc.gz $scaff tahitica:RS_29,RS_30,RS_31,RS_32,RS_33 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.dimidiata`; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/dimidiata.$chrom.$indv.smc.gz $scaff dimidiata:RS_21,RS_22,RS_23,RS_24,RS_25,RS_26,RS_27,RS_28 -d $indv $indv; done; done < ./chromosome-scaffold.table.auto.txt
deactivate
```

#### 2. Run conversions for Z chromosome

Important! Write these results to `out-chrZ` directory.

We'll choose 2 distinguished male individuals per species.

```
source ~/tmp/smc++-install/smc++-env/bin/activate
for indv in HRVN96108 HRVN96106; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/rustica.chrZ.$indv.smc.gz NC_053488.1 rustica:HRVN96108,HRVN96104,HRVN96106,HRVN96102 -d $indv $indv; done
for indv in RS_34 RS_36; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/aethiopica.chrZ.$indv.smc.gz NC_053488.1 aethiopica:RS_34,RS_35,RS_36,RS_37 -d $indv $indv; done
for indv in RS_2 RS_4; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/smithii.chrZ.$indv.smc.gz NC_053488.1 smithii:RS_1,RS_2,RS_3,RS_4,RS_7,RS_8 -d $indv $indv; done
for indv in RS_14 RS_16; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/neoxena.chrZ.$indv.smc.gz NC_053488.1 neoxena:RS_14,RS_15,RS_16,RS_18,RS_20 -d $indv $indv; done
for indv in RS_30 RS_32; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/tahitica.chrZ.$indv.smc.gz NC_053488.1 tahitica:RS_30,RS_31,RS_32,RS_33 -d $indv $indv; done
for indv in RS_22 RS_28; do smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.chrZ.vcf.gz ./out-chrZ/dimidiata.chrZ.$indv.smc.gz NC_053488.1 dimidiata:RS_22,RS_23,RS_25,RS_28 -d $indv $indv; done
deactivate
```

#### 3. Run conversions for two-population inference (i.e., to estimate split times); autosomes only

```
source ~/tmp/smc++-install/smc++-env/bin/activate
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-aethiopica12.$chrom.smc.gz $scaff rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 aethiopica:RS_34,RS_35,RS_36,RS_37,RS_38; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-aethiopica21.$chrom.smc.gz $scaff aethiopica:RS_34,RS_35,RS_36,RS_37,RS_38 rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-smithii12.$chrom.smc.gz $scaff rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 smithii:RS_5,RS_1,RS_2,RS_3,RS_4,RS_6,RS_7,RS_8,RS_9,RS_10; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-smithii21.$chrom.smc.gz $scaff smithii:RS_5,RS_1,RS_2,RS_3,RS_4,RS_6,RS_7,RS_8,RS_9,RS_10 rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-neoxena12.$chrom.smc.gz $scaff rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 neoxena:RS_11,RS_12,RS_13,RS_14,RS_15,RS_16,RS_17,RS_18,RS_19,RS_20; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-neoxena21.$chrom.smc.gz $scaff neoxena:RS_11,RS_12,RS_13,RS_14,RS_15,RS_16,RS_17,RS_18,RS_19,RS_20 rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-dimidiata12.$chrom.smc.gz $scaff rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 dimidiata:RS_21,RS_22,RS_23,RS_24,RS_25,RS_26,RS_27,RS_28; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.rustica-dimidiata21.$chrom.smc.gz $scaff dimidiata:RS_21,RS_22,RS_23,RS_24,RS_25,RS_26,RS_27,RS_28 rustica:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298; done < ./chromosome-scaffold.table.auto.txt
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.neoxena-tahitica12.$chrom.smc.gz $scaff neoxena:RS_11,RS_12,RS_13,RS_14,RS_15,RS_16,RS_17,RS_18,RS_19,RS_20 tahitica:RS_30,RS_31,RS_32,RS_33; smc++ vcf2smc -c 50000 ./vcf/hirundo_genus.snps.$chrom.vcf.gz ./out/two-pop.neoxena-tahitica21.$chrom.smc.gz $scaff tahitica:RS_30,RS_31,RS_32,RS_33 neoxena:RS_11,RS_12,RS_13,RS_14,RS_15,RS_16,RS_17,RS_18,RS_19,RS_20; done < ./chromosome-scaffold.table.auto.txt
deactivate
```

------------------------------------------------------------------------------------------
### 4. Run single population SMC++ analysis

Note: using '-c 50000' flag during estimate helps reduce warnings about long runs of homozygosity

Run SMC++ with cubic splines:
```
source ~/tmp/smc++-install/smc++-env/bin/activate
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base rustica.cubic 4.6e-9 ./out/rustica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base aethiopica.cubic 4.6e-9 ./out/aethiopica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base smithii.cubic 4.6e-9 ./out/smithii.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base neoxena.cubic 4.6e-9 ./out/neoxena.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base tahitica.cubic 4.6e-9 ./out/tahitica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 200000 -o ./analysis/ --base dimidiata.cubic 4.6e-9 ./out/dimidiata.*.smc.gz
```

Plot results and write .csv:
```
smc++ plot ./analysis/rustica.cubic-SMC.pdf ./analysis/rustica.cubic.final.json -g 1 -c
smc++ plot ./analysis/aethiopica.cubic-SMC.pdf ./analysis/aethiopica.cubic.final.json -g 1 -c
smc++ plot ./analysis/smithii.cubic-SMC.pdf ./analysis/smithii.cubic.final.json -g 1 -c
smc++ plot ./analysis/neoxena.cubic-SMC.pdf ./analysis/neoxena.cubic.final.json -g 1 -c
smc++ plot ./analysis/tahitica.cubic-SMC.pdf ./analysis/tahitica.cubic.final.json -g 1 -c
smc++ plot ./analysis/dimidiata.cubic-SMC.pdf ./analysis/dimidiata.cubic.final.json -g 1 -c
```

Note: not specifying '-c 50000' to mask long runs of homozygosity (unsurprisingly) systematically underestimates population size at more recent timescales.

------------------------------------------------------------------------------------------
### 5. Run single population SMC++ analysis on Z chromosome

Run SMC++ with cubic splines:
```
source ~/tmp/smc++-install/smc++-env/bin/activate
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base rustica.cubic-chrZ 4.6e-9 ./out-chrZ/rustica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base aethiopica.cubic-chrZ 4.6e-9 ./out-chrZ/aethiopica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base smithii.cubic-chrZ 4.6e-9 ./out-chrZ/smithii.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base neoxena.cubic-chrZ 4.6e-9 ./out-chrZ/neoxena.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base tahitica.cubic-chrZ 4.6e-9 ./out-chrZ/tahitica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 -o ./analysis/ --base dimidiata.cubic-chrZ 4.6e-9 ./out-chrZ/dimidiata.*.smc.gz
```

Plot results and write .csv:
```
smc++ plot ./analysis/rustica.cubic-chrZ-SMC.pdf ./analysis/rustica.cubic-chrZ.final.json -g 1 -c
smc++ plot ./analysis/aethiopica.cubic-chrZ-SMC.pdf ./analysis/aethiopica.cubic-chrZ.final.json -g 1 -c
smc++ plot ./analysis/smithii.cubic-chrZ-SMC.pdf ./analysis/smithii.cubic-chrZ.final.json -g 1 -c
smc++ plot ./analysis/neoxena.cubic-chrZ-SMC.pdf ./analysis/neoxena.cubic-chrZ.final.json -g 1 -c
smc++ plot ./analysis/tahitica.cubic-chrZ-SMC.pdf ./analysis/tahitica.cubic-chrZ.final.json -g 1 -c
smc++ plot ./analysis/dimidiata.cubic-chrZ-SMC.pdf ./analysis/dimidiata.cubic-chrZ.final.json -g 1 -c
```

------------------------------------------------------------------------------------------
### 5. Run single population SMC++ analysis into more distant past

These single-population marginal estimates will be used for two-population split time inference.

```
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base rustica.cubic-long 4.6e-9 ./out/rustica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base aethiopica.cubic-long 4.6e-9 ./out/aethiopica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base smithii.cubic-long 4.6e-9 ./out/smithii.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base neoxena.cubic-long 4.6e-9 ./out/neoxena.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base tahitica.cubic-long 4.6e-9 ./out/tahitica.*.smc.gz
smc++ estimate --cores 8 --spline cubic -c 50000 --timepoints 1000 5000000 -o ./analysis/ --base dimidiata.cubic-long 4.6e-9 ./out/dimidiata.*.smc.gz
```

------------------------------------------------------------------------------------------
### 5. Run two-population SMC++ analysis

Run SMC++ split to refine the single-population marginal estimates into an estimate of the joint demography:

```
smc++ split --cores 8 --em-iterations 20 -w 100 -o ./analysis/ --base two-pop.rustica-aethiopica.cubic ./analysis/rustica.cubic-long.final.json ./analysis/aethiopica.cubic-long.final.json ./out/two-pop.rustica-aethiopica*.smc.gz
smc++ split --cores 8 --em-iterations 20 -w 100 -o ./analysis/ --base two-pop.rustica-smithii.cubic ./analysis/rustica.cubic-long.final.json ./analysis/smithii.cubic-long.final.json ./out/two-pop.rustica-smithii*.smc.gz
smc++ split --cores 8 --em-iterations 20 -w 100 -o ./analysis/ --base two-pop.rustica-neoxena.cubic ./analysis/rustica.cubic-long.final.json ./analysis/neoxena.cubic-long.final.json ./out/two-pop.rustica-neoxena*.smc.gz
smc++ split --cores 8 --em-iterations 20 -w 100 -o ./analysis/ --base two-pop.rustica-dimidiata.cubic ./analysis/rustica.cubic-long.final.json ./analysis/dimidiata.cubic-long.final.json ./out/two-pop.rustica-dimidiata*.smc.gz
smc++ split --cores 8 --em-iterations 20 -w 100 -o ./analysis/ --base two-pop.neoxena-tahitica.cubic ./analysis/neoxena.cubic-long.final.json ./analysis/tahitica.cubic-long.final.json ./out/two-pop.neoxena-tahitica*.smc.gz
```

Note: toggling the '-w 100' option helped in a couple of cases with estimation at deeper timescales.

Plot results and write .csv:
```
smc++ plot ./analysis/two-pop.rustica-aethiopica.cubic-SMC.pdf ./analysis/two-pop.rustica-aethiopica.cubic.final.json -g 1 -c
smc++ plot ./analysis/two-pop.rustica-smithii.cubic-SMC.pdf ./analysis/two-pop.rustica-smithii.cubic.final.json -g 1 -c
smc++ plot ./analysis/two-pop.rustica-neoxena.cubic-SMC.pdf ./analysis/two-pop.rustica-neoxena.cubic.final.json -g 1 -c
smc++ plot ./analysis/two-pop.rustica-dimidiata.cubic-SMC.pdf ./analysis/two-pop.rustica-dimidiata.cubic.final.json -g 1 -c
smc++ plot ./analysis/two-pop.neoxena-tahitica.cubic-SMC.pdf ./analysis/two-pop.neoxena-tahitica.cubic.final.json -g 1 -c
```

[Back to top](#contents)