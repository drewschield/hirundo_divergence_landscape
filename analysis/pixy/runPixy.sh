list=$1
for chrom in `cat $list`; do
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.1mb-100kb.$chrom.bed --output_folder results --output_prefix pixy.1mb-100kb.$chrom
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.100kb-10kb.$chrom.bed --output_folder results --output_prefix pixy.100kb-10kb.$chrom
	pixy --n_cores 8 --stats pi dxy fst --vcf /media/drewschield/VernalBucket/hirundo/vcf/chrom-specific-genus/hirundo_genus.allsites.final.$chrom.vcf.gz --populations popmap.pixy --bed_file ./windows/window.50kb-5kb.$chrom.bed --output_folder results --output_prefix pixy.50kb-5kb.$chrom
done
