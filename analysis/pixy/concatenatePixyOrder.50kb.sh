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
