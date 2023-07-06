for window in 100kb 10kb 1kb; do
	echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.pi.$window.txt
	for chrom in `cat scaffold.all.list`; do
		cat ./results/pixy_${chrom}_${window}_pi.txt | tail -n +2 >> pixy.all.pi.$window.txt
	done
	echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy no_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy.all.dxy.$window.txt
	for chrom in `cat scaffold.all.list`; do
		cat ./results/pixy_${chrom}_${window}_dxy.txt | tail -n +2 >> pixy.all.dxy.$window.txt
	done
	echo "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps" > pixy.all.fst.$window.txt
	for chrom in `cat scaffold.all.list`; do
		cat ./results/pixy_${chrom}_${window}_fst.txt | tail -n +2 >> pixy.all.fst.$window.txt
	done
done

