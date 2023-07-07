for pop in rustica tytleri gutturalis rustica-tytleri rustica-gutturalis tytleri-gutturalis; do
	cat ./r2/r2.$pop.auto.hap.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}' | sort -k1,1n > ./summary/pop.$pop.auto.hap.ld.summary
	cat ./r2/r2.$pop.chrZ.hap.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}' | sort -k1,1n > ./summary/pop.$pop.chrZ.hap.ld.summary
done

