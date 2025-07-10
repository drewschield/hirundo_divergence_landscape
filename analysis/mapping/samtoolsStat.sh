for indv in `cat sample.list`; do
	echo calculating mapping statistics for $indv
	samtools stat -@ 16 /media/mother/VernalBucket/hirundo/bam/$indv.bam > stats/$indv.stat.txt
done
