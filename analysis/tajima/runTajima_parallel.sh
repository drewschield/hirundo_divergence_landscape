pop=$1
#bcftools view -S popmap.$pop -O z /media/mother/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto+chrZ.snps.vcf.gz | vk tajima 1000000 1000000 -  > ./results/tajima.$pop.1mb.txt
#bcftools view -S popmap.$pop -O z /media/mother/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto+chrZ.snps.vcf.gz | vk tajima 50000 50000 -  > ./results/tajima.$pop.50kb.txt
bcftools view -S popmap.$pop -O z /media/mother/VernalBucket/hirundo/vcf/hirundo_genus.allsites.final.auto+chrZ.snps.vcf.gz | vk tajima 5000 5000 -  > ./results/tajima.$pop.5kb.txt
