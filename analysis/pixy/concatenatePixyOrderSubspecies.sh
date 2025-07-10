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
