#!/bin/R
#CRIPSR-X
#LF
#DURGA
#Mai 2016

#process results from getnumber_mutations_read.py
args <- commandArgs(trailingOnly=T)


#SET WORKING DIRECTORY
wd=as.character(args[1])
setwd(wd)


#GET LIST OF INTEREST FILES
H1=args[2]
H2=args[3]
h1=args[4]
h2=args[5]

par=paste(as.character(args[6]),'_n',as.character(args[7]),'_mapq30_sorted_',as.character(H1),'-',as.character(h2),'.sam_multiguide_mutations_per_reads.txt', sep='')
files<-list.files(pattern=paste('sorted_',as.character(H1),'-',as.character(h2),'.sam_multiguide_mutations_per_reads.txt', sep=''))

#PROCESS EACH FILE AND GET RESULTS

parent<-as.data.frame(read.table(par, header=T, sep="\t"))
results=array(dim=c(length(files),20))
colnames(results)=c('Sample','H1','H2', 'Mut_readshotspot1','total_reads_hotspot1','proportion1', 'enrichment1', 'Ave_mut_read_hotspot1','SD_mut_read_hotspot1','H1','H2', 'Mut_readshotspot2','total_reads_hotspot2','proportion2', 'enrichment2','Ave_mut_read_hotspot2','SD_mut_read_hotspot2','Ave_mut_read','SD_mut_read','mut_both')

for (i in 1: length(files)){
	temp=as.data.frame(read.table(files[i],header=T, sep="\t"))
	results[i,1]=unlist(strsplit(files[i],"_"))[1]
	results[i,2]=H1
	results[i,3]=H2
	results[i,4]=length(which(temp$hotspot1_mut != 0))
	results[i,5]=length(which(temp$read_start>=H1 & temp$read_start<=H2))
	results[i,6]=length(which(temp$hotspot1_mut != 0))/length(which(temp$read_start>=H1 & temp$read_start<=H2))
	results[i,7]=(length(which(temp$hotspot1_mut != 0))/length(which(temp$read_start>=H1 & temp$read_start<=H2)))/(length(which(parent$hotspot1_mut != 0))/length(which(parent$read_start>=H1 & parent$read_start<=H2)))
	results[i,8]=mean(temp$hotspot1_mut[which(temp$hotspot1_mut != 0)])
	results[i,9]=sd(temp$hotspot1_mut[which(temp$hotspot1_mut != 0)])
	results[i,10]=h1
	results[i,11]=h2
	results[i,12]=length(which(temp$hotspot2_mut != 0))
	results[i,13]=length(which(temp$read_start>=h1 & temp$read_start<=h2))
	results[i,14]=length(which(temp$hotspot2_mut != 0))/length(which(temp$read_start>=h1 & temp$read_start<=h2))
	results[i,15]=(length(which(temp$hotspot2_mut != 0))/length(which(temp$read_start>=h1 & temp$read_start<=h2)))/(length(which(parent$hotspot2_mut != 0))/length(which(parent$read_start>=h1 & parent$read_start<=h2)))
	results[i,16]=mean(temp$hotspot2_mut[which(temp$hotspot2_mut != 0)])
	results[i,17]=sd(temp$hotspot2_mut[which(temp$hotspot2_mut != 0)])
	results[i,18]=mean(temp$mutations[which(temp$mutations != 0)])
	results[i,19]=sd(temp$mutations[which(temp$mutations != 0)])
	results[i,20]=length(which(temp$hotspot1_mut != 0 & temp$hotspot2_mut != 0 ))
}

write.table(results,paste('Readstatistics_hotspot_',as.character(H1),'-',as.character(H2),'_',as.character(h1),'-',as.character(h2),'.txt', sep=''), quote=F,append=F,col.names=T,row.names=F,sep="\t")

#Sample
#H1
#H2
#Number of mutated reads in hotspot
#Total reads in hotspot
#ave number of mutations in read in hotspot among mutated reads
#SD
#ave number of mutations in read among mutated reads
#SD
#number of 



