for mod in SI AM IM SC AM2m IM2m SC2m; do
	for i in ${mod}_*; do
		model=`echo $i | cut -d'_' -f1`;
		iter=`echo $i | cut -d'_' -f2`; 
		ll=`grep 'Optimized log-likelihood' $i/$i.txt | tail -n1 | cut -d' ' -f3`;
		aic=`grep 'AIC' $i/$i.txt | tail -n1 | cut -d' ' -f2`;
		echo $model $iter $ll $aic; done | sort -g -k3 | tail -n1;
done

