list=$1
for i in `cat $list`; do
	mosdepth -t 4 --fast-mode -n -b chrom.chr3.10kb.bed -f /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta ./mosdepth_results/$i.chr3 /media/drewschield/VernalBucket/hirundo/bam/$i.bam
	mosdepth -t 4 --fast-mode -n -b chrom.chrZ.10kb.bed -f /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.fasta ./mosdepth_results/$i.chrZ /media/drewschield/VernalBucket/hirundo/bam/$i.bam
	gunzip ./mosdepth_results/$i.chr3.regions.bed.gz 
	gunzip ./mosdepth_results/$i.chrZ.regions.bed.gz
	python identifySex.py mosdepth_results/$i.chr3.regions.bed mosdepth_results/$i.chrZ.regions.bed
done
