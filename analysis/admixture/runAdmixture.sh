ped=$1
for K in 1 2 3 4 5 6 7 8 9 10; do
	admixture -j16 --cv $ped $K | tee log${K}.out
done