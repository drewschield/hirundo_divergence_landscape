list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats pi dxy fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.all.pixy --window 1000000 --output_folder results_all --output_prefix pixy.all.1mb.$chrom
done
