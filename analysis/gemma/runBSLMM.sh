in=$1
rel=$2
burn=$3
run=$4
field=$5
out=$6

for i in $(seq 1 10); do
	gemma -bfile $in -k $rel -bslmm 1 -w $burn -s $run -rpace 100 -wpace 25000 -maf 0.01 -hmin 0 -hmax 1 -rmin 0 -rmax 1 -n $field -o gwas_full_bslmm.$out.run${i}
done