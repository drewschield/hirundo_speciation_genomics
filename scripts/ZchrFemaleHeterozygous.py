"""
identifyFemaleZhetSites.py:
Script to read in list of female sample IDs, then parse a VCF file by female sample columns,
and extract positions on the Z Chromosome with heterozygous genotype calls for ANY female.

python identify_female_Zhet_sites.py <female_list> <VCF> <output>
"""

import sys

out = open(sys.argv[3],'w')

females = []
for line in open(sys.argv[1],'r'):
	indv = line.split()[0]
	females.append(indv)

header_split = []
indexes = []
for line in open(sys.argv[2],'r'):
	if line.split()[0] == '#CHROM':
		for h in line.split():
			header_split.append(h)
			
		for f in females:
			indexes.append(header_split.index(f))
	
	if '##' not in str(line.split()[0]):
		if 'NC_' in str(line.split()[0]):
			genotypes = []
			chrom = line.split()[0]
			site = line.split()[1]
			for i in indexes:
				geno = line.split()[i]
				geno = geno.split(':')[0]
				genotypes.append(geno)
			if '0/1' in str(genotypes):
				out.write(str(chrom)+'\t'+str(site)+'\n')
		if 'NW_' in str(line.split()[0]):
			genotypes = []
			chrom = line.split()[0]
			site = line.split()[1]
			for i in indexes:
				geno = line.split()[i]
				geno = geno.split(':')[0]
				genotypes.append(geno)
			if '0/1' in str(genotypes):
				out.write(str(chrom)+'\t'+str(site)+'\n')