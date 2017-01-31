import sys
import re
from os.path import basename

###READ SAM DATA FOR SELECTED READS
###GET MD FIELD
###SPLIT MD FIELD BY ":"
###GET THE THIRD FIELD
###COUNT NUMBER OF LETTERS
###OUTPUT READ NAME; MD FIELD; NUMBER OF LETTERS

SAM=sys.argv[1]
h1=int(sys.argv[2])
h2=int(sys.argv[3])
OUT=basename(sys.argv[1])+'_mutations_per_reads.txt'

def is16(n):
	if n & 16:
		return True
	else:
		return False

def snp_position(readstart,info, h1,h2):
	#snp_position=0
	#if is16(flag):
	#	snp_position=readstart
#get indexes of letters
	spt=filter(None, re.split(r'(\d+)', info))
	positions=[]
	for i in range(len(spt)):
		if spt[i].isalpha():
			#indexes=[j for j, k in enumerate(spt) if k == i]
			#print spt[i]
			#print i
			pos=int(readstart)+1

			for j in range(i):
				
				#print spt[j]
				if spt[j].isalpha():
					pos+=1
				elif spt[j].startswith('^'):
					pos+=1
				else:
					pos+=int(spt[j])
			positions.append(pos)

		elif spt[i].startswith('^'):
			pos=int(readstart)+1

			for j in range(i):
				
				#print spt[j]
				if spt[j].isalpha():
					pos+=1
				elif spt[j].startswith('^'):
					pos+=1
				else:
					pos+=int(spt[j])
			positions.append(pos)				

	return 	len([s for s in positions if (int(s)>=h1 and int(s)<=h2) ])#, '\t'.join(str(e) for e in positions)

def number_mut(SAM,OUT):
	SAM_FILE=open(SAM)
	OUT_FILE=open(OUT,'w')
	#OUT_FILE.write('readname\tchrom\tread_start\tcigar\tMD\tmutations\thotspot_mut\tmut_pos\n')
	OUT_FILE.write('readname\tchrom\tread_start\tcigar\tMD\tmutations\thotspot_mut\n')

	linenum=0
	for line in SAM_FILE:
		newline=line.rstrip().split('\t')
		readname=newline[0]
		chrom=newline[2]
		readstart=newline[3]
		cigar=newline[5]
		MD=newline[18]
		MD_spl=MD.split(':')
		info=MD_spl[2]

		#hotspot_number,snp_pos=snp_position(readstart,info, h1,h2)
		hotspot_number=snp_position(readstart,info, h1,h2)
		
		letters=0
		other=0

		for i in info:
			if i.isalpha():
				letters+=1
			else:
				other+=1
		#OUT_FILE.write(str(readname) + '\t' + str(chrom) + '\t' + str(readstart) + '\t' + str(cigar) + '\t'+ str(MD) + '\t' + str(letters) +'\t'+ str(hotspot_number)+'\t'+ str(snp_pos)+'\n')
		OUT_FILE.write(str(readname) + '\t' + str(chrom) + '\t' + str(readstart) + '\t' + str(cigar) + '\t'+ str(MD) + '\t' + str(letters) +'\t'+ str(hotspot_number)+'\n')

		linenum+=1
	OUT_FILE.close()


print 'Get number of mutation per read in the hotspot'
sys.stdout.flush()
number_mut(SAM,OUT)
