import sys
import re
from os.path import basename
import csv


POS_FILE=sys.argv[1]
OUT_FILE=sys.argv[1][:-4]+'_summarize_diversity.txt'
#POS_FILE='CX179_n5_mapq30_sorted_485-585_mutreadinhotspot_posfile_sorted.txt'
#OUT_FILE='CX179_n5_mapq30_sorted_485-585_mutreadinhotspot_posfile_sorted.txt_summarize_diversity.txt'
def getmutation_read(POS_FILE):
	POS = open(POS_FILE)
	read_dict = {}
	for line in POS:
		line = line.rstrip().split('\t')
		read=line[1]
		snp=line[0]
		direction=line[3]
		allele=line[6]
		ref_allele=line[7]
		#put alleles same as ref alleles in lower case
		if allele==ref_allele:
			pass
		else:
			if direction not in read_dict.keys():
				read_dict[direction]={}
			if read not in read_dict[direction].keys():
				read_dict[direction][read]={}
			if snp not in read_dict[direction][read].keys():
				read_dict[direction][read][snp]={}
			read_dict[direction][read][snp]=allele
	return read_dict
	POS.close()

#
def getlist_combinations(read_dict):
	snp_alleles=[]
	#read_dict=getmutation_read(POS_FILE)
	for direction in read_dict:
		for readname in read_dict[direction]:
			i=[]
			for key, value in sorted(read_dict[direction][readname].items()):
				tup=(key,value)
				i.append(tup)
			i=tuple(i)
			snp_alleles.append(i)
	return snp_alleles

def get_summary_combination(snp_alleles):
	mut_dict={}
	#snp_alleles=getlist_combinations(read_dict)
	for combination in snp_alleles:
		if combination not in mut_dict:
			mut_dict[combination]=1
		else:
			mut_dict[combination]+=1
	return mut_dict



def format_results(mut_dict):
	OUT=open(OUT_FILE,'w')
	for key, value in mut_dict.items():
		
		toto=[]
		for el in key:
			toto.append('-'.join(el))
			#titi=('-').join(toto)
		OUT.write(str(value)+'	'+'	'.join(toto)+'\n')
	OUT.close()
	

##Get list of mutations per read
print 'Get Mutations of each read'
sys.stdout.flush()
read_dict = getmutation_read(POS_FILE)


##Get list of possible combinations of mutations
print 'Get list of possible combinations of mutations'
sys.stdout.flush()
snp_alleles = getlist_combinations(read_dict)

##Get summary of combinations of mutations
print 'Get summary of combinations of mutations'
sys.stdout.flush()
mut_dict = get_summary_combination(snp_alleles)

##Format results and write output
print 'Format results and write output'
sys.stdout.flush()
format_results(mut_dict)
