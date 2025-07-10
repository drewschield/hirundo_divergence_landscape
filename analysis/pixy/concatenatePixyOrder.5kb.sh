echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.all.5kb.fst.txt
for chrom in `cat scaffold.order.list`; do
	cat ./results_all/pixy.all.5kb.${chrom}_fst.txt | tail -n +2 >> pixy.all.5kb.fst.txt
done
