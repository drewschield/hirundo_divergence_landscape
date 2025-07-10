for pop in rustica aethiopica smithii neoxena tahitica dimidiata; do
	pyrho optimize --numthreads 8 --tablefile ./lookup/${pop}_lookuptable.hdf --vcffile ./vcf/$pop.snps.chr1A.vcf.gz --outfile ./test/$pop.chr1A.b10-w25.rmap --blockpenalty 10 --windowsize 25 --ploidy 2 --logfile .
	pyrho optimize --numthreads 8 --tablefile ./lookup/${pop}_lookuptable.hdf --vcffile ./vcf/$pop.snps.chr1A.vcf.gz --outfile ./test/$pop.chr1A.b10-w50.rmap --blockpenalty 10 --windowsize 50 --ploidy 2 --logfile .
	pyrho optimize --numthreads 8 --tablefile ./lookup/${pop}_lookuptable.hdf --vcffile ./vcf/$pop.snps.chr1A.vcf.gz --outfile ./test/$pop.chr1A.b20-w25.rmap --blockpenalty 20 --windowsize 25 --ploidy 2 --logfile .
	pyrho optimize --numthreads 8 --tablefile ./lookup/${pop}_lookuptable.hdf --vcffile ./vcf/$pop.snps.chr1A.vcf.gz --outfile ./test/$pop.chr1A.b20-w50.rmap --blockpenalty 20 --windowsize 50 --ploidy 2 --logfile .
done
