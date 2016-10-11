import numpy
import sys

#SMP_FILE='CX179_n5_mapq30_sorted_485-585_mutreadinhotspot_posfile_sorted.txt_summarize_diversity.txt'
SMP_FILE=sys.argv[1]
#PARENT_FILE='CX188_n5_mapq30_sorted_485-585_mutreadinhotspot_posfile_sorted.txt_summarize_diversity.txt'
PARENT_FILE=sys.argv[2]

reads_sample=int(sys.argv[3])
#reads_sample=533884
reads_parent=int(sys.argv[4])
#reads_parent=453828

OUT_FILE1=sys.argv[1][:-4]+'_stats_diversity.txt'
OUT_FILE=sys.argv[1][:-4]+'_combinations_diversity_filtered.txt'

def make_dictionary(myfile):
	SMP=open(myfile)
	sample={}
	for line in SMP:
		line = line.rstrip().split('\t')
		total=line[0]
		combi=tuple(line[1:])
		if combi not in sample.keys():
			sample[combi]=total
	return sample



#def allkeys_summary(sample,parent):
#	all_keys=set(sample.keys()) | set(parent.keys())
#	summary={}
#	for key in all_keys:
#		if key not in summary:
#			if key not in sample:
#				summary[key]=-int(parent[key])
#			elif key not in parent:
#				 summary[key]=int(sample[key])
#			else:
#				summary[key]=int(sample[key])-int(parent[key])
#	return summary


def get_stats_values(dic,reads):
	val=[]
	OUT_stats=open(OUT_FILE1,'w')
	for value in dic.values():
		val.append(float(value)/reads)
	OUT_stats.write('Mean_frequency_of_combinations_in_parent '+str(numpy.mean(val))+'\n'+'Standard_deviation_frequency_of_combinations_in_parent '+str(numpy.std(val))+'\n')
	OUT_stats.close()
	return numpy.mean(val), numpy.std(val) 


def get_values_over_threshold_sample(dic,average,SD,reads_sample):
	val={}
	for key, value in dic.items():
		if (float(value)/reads_sample)-average >=SD:
			val[key]=int(value)
	
	return val 

def filter_values(val,parent_dict,reads_sample, reads_parent, SD):
	OUT=open(OUT_FILE,'w')
	OUT.write('combination\tsample_frequency\tparent_frequency\n')
	OUT_stats=open(OUT_FILE1,'a+')
	valok=0
	for key, value in val.items():
		if key in parent_dict:
			par_val=parent_dict[key]
		else:
			par_val=0
		if float(value)/reads_sample-float(par_val)/reads_parent >= SD:
			valok+=1
			OUT.write(str(key)+'	'+str(float(value)/reads_sample)+'	'+str(float(par_val)/reads_parent)+'\n')
	OUT_stats.write('Diversity_after_filtering_is '+ str(valok)+'\n')
	OUT_stats.close()
	OUT.close()





#make_list_values_sample(sample)

def format_results(summary):
	OUT=open(OUT_FILE,'w')
	
	for key, value in summary.items():
		
		toto=[]
		for el in key:
			toto.append(el)
			#titi=('-').join(toto)
		print(str(value)+'	'+'	'.join(toto)+'\n')
	OUT.close()
	

##Make dictionnary of combination files
print 'Make dictionnary of combination files'

sample=make_dictionary(SMP_FILE)
parent=make_dictionary(PARENT_FILE)


##Get statistics on parent combination frequency
print 'Get statistics on parent combination frequency'
average, SD=get_stats_values(parent, reads_parent)


##Filter sample combinations for which frequency - mean parent frequency is over SD of frequency in parent
print 'Filter sample combinations for which frequency - mean parent frequency is over SD of frequency in parent'

first_filter=get_values_over_threshold_sample(sample,average,SD,reads_sample)

##Filter sample combinations for which frequency - parent frequency is over SD of frequency in parent
print 'Filter sample combinations for which frequency - parent frequency is over SD of frequency in parent'

filter_values(first_filter,parent,reads_sample, reads_parent, SD)
#make_list_values(parent, reads_parent)

#shared_keys=set(sample.keys()) & set(parent.keys())
#
#summary={}
#for key in shared_keys:
#	if key not in summary:
#		summary[key]=int(sample[key])-int(parent[key])

