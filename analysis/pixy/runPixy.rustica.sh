list=$1
for chrom in `cat $list`; do
	pixy --n_cores 16 --stats pi dxy fst --vcf /media/mother/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.rustica.pixy --window 1000000 --output_folder results_rustica --output_prefix pixy.rustica.1mb.$chrom
done
