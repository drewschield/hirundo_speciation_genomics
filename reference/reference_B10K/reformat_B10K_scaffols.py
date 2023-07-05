import sys
from Bio import SeqIO

out = open(sys.argv[3],'w')
scaff_list = []
reformat = []
sort = ['chromosome 1','chromosome 1A','chromosome 2','chromosome 3','chromosome 4','chromosome 4A','chromosome 5','chromosome 6','chromosome 7','chromosome 8','chromosome 9','chromosome 10','chromosome 11','chromosome 12','chromosome 13','chromosome 14','chromosome 15','chromosome 17','chromosome 18','chromosome 19','chromosome 20','chromosome 21','chromosome 22','chromosome 23','chromosome 24','chromosome 25','chromosome 26','chromosome 27','chromosome 28','chromosome 29','chromosome 31','chromosome 32','chromosome 33','chromosome 34','chromosome 35','chromosome 36','chromosome 37','chromosome Z','chromosome W']

with open(sys.argv[1],'r') as input:
	next(input)
	for line in input:
		b10k_scaff = line.split('\t')[0]
		scaff_list.append(b10k_scaff)
		b10k_chrom = line.split('\t')[1]
		if "unlocalized" in b10k_chrom:
			b10k_chrom = b10k_chrom.split(',')[0]
		finch_chrom = line.split('\t')[3].strip()

		entry = b10k_scaff+'-'+b10k_chrom+'-'+finch_chrom
		reformat.append(entry)

for chrom in sort:
	for i in SeqIO.parse(sys.argv[2],'fasta'):
		header = i.description
		code = i.seq
		scaff = str(header).split()[0]
		prev = str(header).split('bHirRus1 ')[1]
		if "unlocalized" in str(header):
			prev = prev.split(' unlocalized')[0]
		else:
			prev = prev.split(',')[0]
		for r in reformat:
			if scaff == r.split('-')[0] and prev == r.split('-')[1]:
				corr = r.split('-')[2]
				header_fix = str(header).replace(prev,corr)
				if corr == chrom:
					i.description = header_fix
					print chrom, header_fix
					SeqIO.write(i,out,"fasta")

for i in SeqIO.parse(sys.argv[2],'fasta'):
	scaff = i.id
	if scaff not in scaff_list:
		SeqIO.write(i,out,"fasta")