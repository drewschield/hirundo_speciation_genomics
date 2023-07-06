while read i; do
	scaff=`echo "$i" | cut -f 1`
	chr=`echo "$i" | cut -f 2`
	echo "HRVN96101	/media/drewschield/VernalBucket/hirundo/bam/HRVN96101.bam	$scaff" > bamlist.$chr
	echo "HRVN96107	/media/drewschield/VernalBucket/hirundo/bam/HRVN96107.bam	$scaff" >> bamlist.$chr
	echo "HRVN96108	/media/drewschield/VernalBucket/hirundo/bam/HRVN96108.bam	$scaff" >> bamlist.$chr
	echo "HRVN96103	/media/drewschield/VernalBucket/hirundo/bam/HRVN96103.bam	$scaff" >> bamlist.$chr
	echo "HRVN96104	/media/drewschield/VernalBucket/hirundo/bam/HRVN96104.bam	$scaff" >> bamlist.$chr
	echo "HRVN96105	/media/drewschield/VernalBucket/hirundo/bam/HRVN96105.bam	$scaff" >> bamlist.$chr
	echo "HRVN96106	/media/drewschield/VernalBucket/hirundo/bam/HRVN96106.bam	$scaff" >> bamlist.$chr
	echo "HRVN96300	/media/drewschield/VernalBucket/hirundo/bam/HRVN96300.bam	$scaff" >> bamlist.$chr
	echo "HRVN96102	/media/drewschield/VernalBucket/hirundo/bam/HRVN96102.bam	$scaff" >> bamlist.$chr
	echo "HRVN96298	/media/drewschield/VernalBucket/hirundo/bam/HRVN96298.bam	$scaff" >> bamlist.$chr
	echo "HRVN96110	/media/drewschield/VernalBucket/hirundo/bam/HRVN96110.bam	$scaff" >> bamlist.$chr
	echo "HRVN96113	/media/drewschield/VernalBucket/hirundo/bam/HRVN96113.bam	$scaff" >> bamlist.$chr
	echo "HRVN96116	/media/drewschield/VernalBucket/hirundo/bam/HRVN96116.bam	$scaff" >> bamlist.$chr
	echo "HRVN96117	/media/drewschield/VernalBucket/hirundo/bam/HRVN96117.bam	$scaff" >> bamlist.$chr
	echo "HRVN96118	/media/drewschield/VernalBucket/hirundo/bam/HRVN96118.bam	$scaff" >> bamlist.$chr
	echo "HRVN96128	/media/drewschield/VernalBucket/hirundo/bam/HRVN96128.bam	$scaff" >> bamlist.$chr
	echo "HRVN96126	/media/drewschield/VernalBucket/hirundo/bam/HRVN96126.bam	$scaff" >> bamlist.$chr
	echo "HRVN96124	/media/drewschield/VernalBucket/hirundo/bam/HRVN96124.bam	$scaff" >> bamlist.$chr
	echo "HRVN96123	/media/drewschield/VernalBucket/hirundo/bam/HRVN96123.bam	$scaff" >> bamlist.$chr
	echo "HRVN96122	/media/drewschield/VernalBucket/hirundo/bam/HRVN96122.bam	$scaff" >> bamlist.$chr
	echo "HR100	/media/drewschield/VernalBucket/hirundo/bam/HR100.bam	$scaff" >> bamlist.$chr
	echo "HR113	/media/drewschield/VernalBucket/hirundo/bam/HR113.bam	$scaff" >> bamlist.$chr
	echo "HR130	/media/drewschield/VernalBucket/hirundo/bam/HR130.bam	$scaff" >> bamlist.$chr
	echo "HR131	/media/drewschield/VernalBucket/hirundo/bam/HR131.bam	$scaff" >> bamlist.$chr
	echo "HR1646	/media/drewschield/VernalBucket/hirundo/bam/HR1646.bam	$scaff" >> bamlist.$chr
	echo "HR1648	/media/drewschield/VernalBucket/hirundo/bam/HR1648.bam	$scaff" >> bamlist.$chr
	echo "HR3	/media/drewschield/VernalBucket/hirundo/bam/HR3.bam	$scaff" >> bamlist.$chr
	echo "HR15	/media/drewschield/VernalBucket/hirundo/bam/HR15.bam	$scaff" >> bamlist.$chr
	echo "HR27	/media/drewschield/VernalBucket/hirundo/bam/HR27.bam	$scaff" >> bamlist.$chr
	echo "HR43	/media/drewschield/VernalBucket/hirundo/bam/HR43.bam	$scaff" >> bamlist.$chr
	echo "HR48	/media/drewschield/VernalBucket/hirundo/bam/HR48.bam	$scaff" >> bamlist.$chr
	echo "HRNB030	/media/drewschield/VernalBucket/hirundo/bam/HRNB030.bam	$scaff" >> bamlist.$chr
	echo "HRNB022	/media/drewschield/VernalBucket/hirundo/bam/HRNB022.bam	$scaff" >> bamlist.$chr
	echo "HRNB031	/media/drewschield/VernalBucket/hirundo/bam/HRNB031.bam	$scaff" >> bamlist.$chr
	echo "HRNB027	/media/drewschield/VernalBucket/hirundo/bam/HRNB027.bam	$scaff" >> bamlist.$chr
	echo "HRNB028	/media/drewschield/VernalBucket/hirundo/bam/HRNB028.bam	$scaff" >> bamlist.$chr
	echo "HRNB040	/media/drewschield/VernalBucket/hirundo/bam/HRNB040.bam	$scaff" >> bamlist.$chr
	echo "HRNB033	/media/drewschield/VernalBucket/hirundo/bam/HRNB033.bam	$scaff" >> bamlist.$chr
	echo "HRNB019	/media/drewschield/VernalBucket/hirundo/bam/HRNB019.bam	$scaff" >> bamlist.$chr
	echo "HRNB038	/media/drewschield/VernalBucket/hirundo/bam/HRNB038.bam	$scaff" >> bamlist.$chr
	echo "HRNB029	/media/drewschield/VernalBucket/hirundo/bam/HRNB029.bam	$scaff" >> bamlist.$chr
	echo "HRNB041	/media/drewschield/VernalBucket/hirundo/bam/HRNB041.bam	$scaff" >> bamlist.$chr
	echo "HR302	/media/drewschield/VernalBucket/hirundo/bam/HR302.bam	$scaff" >> bamlist.$chr
	echo "HR303	/media/drewschield/VernalBucket/hirundo/bam/HR303.bam	$scaff" >> bamlist.$chr
	echo "HR307	/media/drewschield/VernalBucket/hirundo/bam/HR307.bam	$scaff" >> bamlist.$chr
	echo "HR309	/media/drewschield/VernalBucket/hirundo/bam/HR309.bam	$scaff" >> bamlist.$chr
	echo "HR316	/media/drewschield/VernalBucket/hirundo/bam/HR316.bam	$scaff" >> bamlist.$chr
	echo "HR318	/media/drewschield/VernalBucket/hirundo/bam/HR318.bam	$scaff" >> bamlist.$chr
	echo "HR319	/media/drewschield/VernalBucket/hirundo/bam/HR319.bam	$scaff" >> bamlist.$chr
	echo "HR324	/media/drewschield/VernalBucket/hirundo/bam/HR324.bam	$scaff" >> bamlist.$chr
	echo "HR326	/media/drewschield/VernalBucket/hirundo/bam/HR326.bam	$scaff" >> bamlist.$chr
	echo "HR327	/media/drewschield/VernalBucket/hirundo/bam/HR327.bam	$scaff" >> bamlist.$chr
	echo "HR1086	/media/drewschield/VernalBucket/hirundo/bam/HR1086.bam	$scaff" >> bamlist.$chr
	echo "HR1097	/media/drewschield/VernalBucket/hirundo/bam/HR1097.bam	$scaff" >> bamlist.$chr
	echo "HR1102	/media/drewschield/VernalBucket/hirundo/bam/HR1102.bam	$scaff" >> bamlist.$chr
	echo "HR1094	/media/drewschield/VernalBucket/hirundo/bam/HR1094.bam	$scaff" >> bamlist.$chr
	echo "HR1099	/media/drewschield/VernalBucket/hirundo/bam/HR1099.bam	$scaff" >> bamlist.$chr
	echo "HR1093	/media/drewschield/VernalBucket/hirundo/bam/HR1093.bam	$scaff" >> bamlist.$chr
	echo "HR1103	/media/drewschield/VernalBucket/hirundo/bam/HR1103.bam	$scaff" >> bamlist.$chr
	echo "HR1105	/media/drewschield/VernalBucket/hirundo/bam/HR1105.bam	$scaff" >> bamlist.$chr
	echo "HR1098	/media/drewschield/VernalBucket/hirundo/bam/HR1098.bam	$scaff" >> bamlist.$chr
	echo "HR1100	/media/drewschield/VernalBucket/hirundo/bam/HR1100.bam	$scaff" >> bamlist.$chr
	echo "HR1089	/media/drewschield/VernalBucket/hirundo/bam/HR1089.bam	$scaff" >> bamlist.$chr
	echo "HR1104	/media/drewschield/VernalBucket/hirundo/bam/HR1104.bam	$scaff" >> bamlist.$chr
	echo "HR1091	/media/drewschield/VernalBucket/hirundo/bam/HR1091.bam	$scaff" >> bamlist.$chr
	echo "HR1095	/media/drewschield/VernalBucket/hirundo/bam/HR1095.bam	$scaff" >> bamlist.$chr
	echo "HR1003	/media/drewschield/VernalBucket/hirundo/bam/HR1003.bam	$scaff" >> bamlist.$chr
	echo "HR1015	/media/drewschield/VernalBucket/hirundo/bam/HR1015.bam	$scaff" >> bamlist.$chr
	echo "HR1019	/media/drewschield/VernalBucket/hirundo/bam/HR1019.bam	$scaff" >> bamlist.$chr
	echo "HR1021	/media/drewschield/VernalBucket/hirundo/bam/HR1021.bam	$scaff" >> bamlist.$chr
	echo "HR1036	/media/drewschield/VernalBucket/hirundo/bam/HR1036.bam	$scaff" >> bamlist.$chr
	echo "HR1040	/media/drewschield/VernalBucket/hirundo/bam/HR1040.bam	$scaff" >> bamlist.$chr
	echo "HR1041	/media/drewschield/VernalBucket/hirundo/bam/HR1041.bam	$scaff" >> bamlist.$chr
	echo "HR1051	/media/drewschield/VernalBucket/hirundo/bam/HR1051.bam	$scaff" >> bamlist.$chr
	echo "HR1058	/media/drewschield/VernalBucket/hirundo/bam/HR1058.bam	$scaff" >> bamlist.$chr
	echo "HR1075	/media/drewschield/VernalBucket/hirundo/bam/HR1075.bam	$scaff" >> bamlist.$chr
	echo "HR1107	/media/drewschield/VernalBucket/hirundo/bam/HR1107.bam	$scaff" >> bamlist.$chr
	echo "HR1158	/media/drewschield/VernalBucket/hirundo/bam/HR1158.bam	$scaff" >> bamlist.$chr
	echo "HR1162	/media/drewschield/VernalBucket/hirundo/bam/HR1162.bam	$scaff" >> bamlist.$chr
	echo "HR1171	/media/drewschield/VernalBucket/hirundo/bam/HR1171.bam	$scaff" >> bamlist.$chr
	echo "HR1177	/media/drewschield/VernalBucket/hirundo/bam/HR1177.bam	$scaff" >> bamlist.$chr
	echo "HR1198	/media/drewschield/VernalBucket/hirundo/bam/HR1198.bam	$scaff" >> bamlist.$chr
	echo "HR1223	/media/drewschield/VernalBucket/hirundo/bam/HR1223.bam	$scaff" >> bamlist.$chr
	echo "HR1224	/media/drewschield/VernalBucket/hirundo/bam/HR1224.bam	$scaff" >> bamlist.$chr
	echo "HR1777	/media/drewschield/VernalBucket/hirundo/bam/HR1777.bam	$scaff" >> bamlist.$chr
	echo "HR1791	/media/drewschield/VernalBucket/hirundo/bam/HR1791.bam	$scaff" >> bamlist.$chr
	echo "HR1798	/media/drewschield/VernalBucket/hirundo/bam/HR1798.bam	$scaff" >> bamlist.$chr
	echo "HR1808	/media/drewschield/VernalBucket/hirundo/bam/HR1808.bam	$scaff" >> bamlist.$chr
	echo "HR1822	/media/drewschield/VernalBucket/hirundo/bam/HR1822.bam	$scaff" >> bamlist.$chr
	echo "HR1838	/media/drewschield/VernalBucket/hirundo/bam/HR1838.bam	$scaff" >> bamlist.$chr
	echo "HR1848	/media/drewschield/VernalBucket/hirundo/bam/HR1848.bam	$scaff" >> bamlist.$chr
done < ../chromosome-scaffold.table.txt
