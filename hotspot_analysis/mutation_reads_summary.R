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

par=paste(as.character(args[4]),'_n',as.character(args[5]),'_mapq30_sorted_',as.character(H1),'-',as.character(H2),'.sam_mutations_per_reads.txt', sep='')
files<-list.files(pattern=paste('sorted_',as.character(H1),'-',as.character(H2),'.sam_mutations_per_reads.txt', sep=''))

#PROCESS EACH FILE AND GET RESULTS

parent<-as.data.frame(read.table(par, header=T, sep="\t"))
results=array(dim=c(length(files),11))
colnames(results)=c('Sample','H1','H2', 'Mut_readshotspot','total_reads_hotspot','proportion', 'enrichment', 'Ave_mut_read_hotspot','SD_mut_read_hotspot','Ave_mut_read','SD_mut_read')

for (i in 1: length(files)){
	temp=as.data.frame(read.table(files[i],header=T, sep="\t"))
	results[i,1]=unlist(strsplit(files[i],"_"))[1]
	results[i,2]=H1
	results[i,3]=H2
	results[i,4]=length(which(temp$hotspot_mut != 0))
	results[i,5]=dim(temp)[1]
	results[i,6]=length(which(temp$hotspot_mut != 0))/dim(temp)[1]
	results[i,7]=(length(which(temp$hotspot_mut != 0))/dim(temp)[1])/(length(which(parent$hotspot_mut != 0))/dim(parent)[1])
	results[i,8]=mean(temp$hotspot_mut[which(temp$hotspot_mut != 0)])
	results[i,9]=sd(temp$hotspot_mut[which(temp$hotspot_mut != 0)])
	results[i,10]=mean(temp$mutations[which(temp$mutations != 0)])
	results[i,11]=sd(temp$mutations[which(temp$mutations != 0)])
}

write.table(results,paste('Readstatistics_hotspot_',as.character(H1),'-',as.character(H2),'.txt', sep=''), quote=F,append=F,col.names=T,row.names=F,sep="\t")

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



