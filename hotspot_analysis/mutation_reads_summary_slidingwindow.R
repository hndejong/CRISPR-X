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

par=as.character(args[2])
sample=as.character(args[3])

#GET LIST OF INTEREST FILES
#H1=args[2]
#H2=args[3]

# Get list of files to be analyzed
par_files=list.files(pattern=paste0('CX',par))
sample_files=list.files(pattern=paste0('CX',sample))


# remove log files from list of files
par_files=par_files[!grepl("log", par_files)]
sample_files=sample_files[!grepl("log", sample_files)]

# PROCESS EACH FILE AND GET RESULTS

#parent<-as.data.frame(read.table(par, header=T, sep="\t"))
results=array(dim=c(length(sample_files),11))
colnames(results)=c('Sample','H1','H2', 'Mut_readshotspot','total_reads_hotspot','proportion', 'enrichment', 'Ave_mut_read_hotspot','SD_mut_read_hotspot','Ave_mut_read','SD_mut_read')

for (i in 1: length(sample_files)){
	temp=as.data.frame(read.table(par_files[i],header=T, sep="\t"))
	parent=as.data.frame(read.table(sample_files[i],header=T, sep="\t"))
	results[i,1]=unlist(strsplit(sample_files[i],"_"))[1]
	results[i,2]=unlist(strsplit(unlist(strsplit(unlist(strsplit(par_files[i], '_'))[5],'[.]'))[1],'-'))[1]
	results[i,3]=unlist(strsplit(unlist(strsplit(unlist(strsplit(par_files[i], '_'))[5],'[.]'))[1],'-'))[2]
	results[i,4]=length(which(temp$hotspot_mut != 0))
	results[i,5]=dim(temp)[1]
	results[i,6]=length(which(temp$hotspot_mut != 0))/dim(temp)[1]
	results[i,7]=(length(which(temp$hotspot_mut != 0))/dim(temp)[1])/(length(which(parent$hotspot_mut != 0))/dim(parent)[1])
	results[i,8]=mean(temp$hotspot_mut[which(temp$hotspot_mut != 0)])
	results[i,9]=sd(temp$hotspot_mut[which(temp$hotspot_mut != 0)])
	results[i,10]=mean(temp$mutations[which(temp$mutations != 0)])
	results[i,11]=sd(temp$mutations[which(temp$mutations != 0)])
}

write.table(results,paste('CX',as.character(sample),'_readstatistics_hotspot_slidingwindow','.txt', sep=''), quote=F,append=F,col.names=T,row.names=F,sep="\t")

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



