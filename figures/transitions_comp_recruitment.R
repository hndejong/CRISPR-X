#!/bin/R
#LF
#November 2016
#durga

# This scripts performs the analysis of transitions obtained with CRIPR-X system with different recruitment methods

#------------- LIBRARIES
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)

#------------- COLORS
palette=brewer.pal(9,"Set1")
mycol_DNA=c(palette[2],palette[3],palette[1],brewer.pal(8,"Set2")[6],palette[4],palette[9])


#------------ FUNCTIONS

plot_transitions_points=function(df, pattern){
	df_short=df[df$recruit==as.character(pattern),]
	#df_short.m=melt(df_short,measure.vars=3:8)
	plot=ggplot(df_short, aes(x=REF,y=value, colour=variable))
	plot=plot+geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=I(1/3))
	plot=plot+facet_grid(.~REF, scales="free")
	plot=plot+ scale_colour_manual(values=mycol_DNA,name="Transition/Transversion",labels=c("A", "G","C","T","deletion","insertion"))
	plot=plot+labs(title=as.character(pattern),x='Reference', y='Frequency of Alternative Alleles')+ theme_bw()
	plot=plot+ theme(axis.title.x=element_text(size=15)) 
	plot=plot+theme(axis.title.y=element_text(size=15))+ theme(legend.text=element_text(size=12))
	plot=plot+ theme(legend.title=element_text(size=15))
	plot=plot+ theme(axis.text.x= element_text(size=10))
	plot=plot+ theme(axis.text.y= element_text(size=10))
	plot=plot+theme(axis.text.x=element_blank())
	return(plot)
}

#------------- WORKING DIRECTORY
setwd('/users/lfresard/CRISPR-X/analysis/transition_recruitment/CX338-351')

#------------- MAIN

# list in count files

counts_file<-list.files(pattern="count")

big_enrichment2=c()
parents=c()
for (i in 1:length(counts_file)){
	temp.count<-read.table(counts_file[i], header=T)
	n<-as.numeric(gsub("CX","", unlist(strsplit(counts_file[i],"_"))[1]))
	if (n >=338 && n<=343){
		chr="GFP684"
		parent=read.table("./parents/CX344_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if( n>=345 && n<=350){
		chr="HBG2"
		parent=read.table("./parents/CX351_n3_mapq30_sorted_qual30.count", header=T)
	}
	parent=parent[parent$CHR==chr,]
	parent$A_freq=parent$A/parent$Depth
	parent$G_freq=parent$G/parent$Depth
	parent$C_freq=parent$C/parent$Depth
	parent$T_freq=parent$T/parent$Depth
	parent$del_freq=parent$del/parent$Depth
	parent$ins_freq=parent$ins/parent$Depth

	temp.count=temp.count[temp.count$CHR==chr,]
	temp.count$A_freq=temp.count$A/temp.count$Depth
	temp.count$G_freq=temp.count$G/temp.count$Depth
	temp.count$C_freq=temp.count$C/temp.count$Depth
	temp.count$T_freq=temp.count$T/temp.count$Depth
	temp.count$del_freq=temp.count$del/temp.count$Depth
	temp.count$ins_freq=temp.count$ins/temp.count$Depth
	if (n >=338 && n<=343){
		temp.count=temp.count[temp.count$bp >=130 & temp.count$bp <=230,]
		parent=parent[parent$bp >=130 & parent$bp <=230,]

	}
	else if (n >=345 && n<=350){
		temp.count=temp.count[temp.count$bp >=275 & temp.count$bp <=375,]
		parent=parent[parent$bp >=275 & parent$bp <=375,]
	}
	temp.count=temp.count[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
	for (j in 1:(dim(temp.count)[1])){
		for (k in 3:8){
			if (abs(temp.count[j,k]-temp.count[j,9])<10e-12){temp.count[j,k]=0}
		}
	}
	parent=parent[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
	for (j in 1:(dim(parent)[1])){
		for (k in 3:8){
			if (abs(parent[j,k]-parent[j,9])<10e-10){parent[j,k]=0}
		}
	}
	parent=parent[which(parent$X._REF>=0.9),]
	bp=intersect(temp.count$bp, parent$bp)
	temp.count=temp.count[which(temp.count$bp %in% bp),]
	parent=parent[which(parent$bp %in% bp),]


	enrichment2=cbind(temp.count[,1:2],temp.count[,3:8]-parent[,3:8])
	enrichment2=as.data.frame(enrichment2)
	enrichment2$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(enrichment2)[1])

	if (n==338 || n==339 || n==345 ||n==346){
		enrichment2$recruit=rep('MS2', dim(enrichment2)[1])
	}else if ( n==340 ||n==341|| n==347 || n==348){
		enrichment2$recruit=rep('C-Fusion', dim(enrichment2)[1])
	}
	else if (n==342 ||n==343 || n==349 || n==350){
		enrichment2$recruit=rep('N-Fusion', dim(enrichment2)[1])
	}
	parents=rbind(parents,parent)
	big_enrichment2=rbind(big_enrichment2,enrichment2)
}

big_enrichment2.m=melt(big_enrichment2,measure.vars=3:8)
big_enrichment2.m$value[which(big_enrichment2.m$value <=0)]=0
parents.m=melt(parents,measure.vars=3:8)

MS2_plot=plot_transitions_points(big_enrichment2.m,'MS2')
C_fusion_plot=plot_transitions_points(big_enrichment2.m,'C-Fusion')
N_fusion_plot=plot_transitions_points(big_enrichment2.m,'N-Fusion')

pdf('transitions_recruitment.pdf')
MS2_plot
C_fusion_plot
N_fusion_plot
dev.off()


CX338=big_enrichment2[big_enrichment2$sample=='CX338',]
CX338.m=melt(CX338, measure.vars=3:8)
CX338.m$value[which(CX338.m$value <=0)]=0

g=ggplot(CX338.m, aes(x=bp,y=value/sum(value),fill=variable))+geom_bar(stat="identity")
p = ggplot(CX338.m, aes(x = bp)) + 
    geom_bar(aes(CX338.m,y = (..value..)/sum(..value..))) + 
    scale_y_continuous(labels = percent)



######

counts_file<-list.files(pattern="count")

percent_change=c()
parents=c()
for (i in 1:length(counts_file)){
	temp.count<-read.table(counts_file[i], header=T)
	n<-as.numeric(gsub("CX","", unlist(strsplit(counts_file[i],"_"))[1]))
	if (n >=338 && n<=343){#next
		chr="GFP684"
		parent=read.table("./parents/CX344_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if( n>=345 && n<=350){
		chr="HBG2"
		parent=read.table("./parents/CX351_n3_mapq30_sorted_qual30.count", header=T)
	}
	parent=parent[parent$CHR==chr,]
	parent$A_freq=parent$A/parent$Depth
	parent$G_freq=parent$G/parent$Depth
	parent$C_freq=parent$C/parent$Depth
	parent$T_freq=parent$T/parent$Depth
	parent$del_freq=parent$del/parent$Depth
	parent$ins_freq=parent$ins/parent$Depth

	temp.count=temp.count[temp.count$CHR==chr,]
	temp.count$A_freq=temp.count$A/temp.count$Depth
	temp.count$G_freq=temp.count$G/temp.count$Depth
	temp.count$C_freq=temp.count$C/temp.count$Depth
	temp.count$T_freq=temp.count$T/temp.count$Depth
	temp.count$del_freq=temp.count$del/temp.count$Depth
	temp.count$ins_freq=temp.count$ins/temp.count$Depth
	if (n >=338 && n<=343){
		temp.count=temp.count[temp.count$bp >=130 & temp.count$bp <=230,]
		parent=parent[parent$bp >=130 & parent$bp <=230,]

	}
	else if (n >=345 && n<=350){
		temp.count=temp.count[temp.count$bp >=275 & temp.count$bp <=375,]
		parent=parent[parent$bp >=275 & parent$bp <=375,]
	}
	for (j in 1:(dim(temp.count)[1])){
		for (k in 15:20){
			if (abs(temp.count[j,k]-temp.count[j,14])<10e-12){temp.count[j,k]=0}
		}
	}
	#parent=parent[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
	for (j in 1:(dim(parent)[1])){
		for (k in 15:20){
			if (abs(parent[j,k]-parent[j,14])<10e-10){parent[j,k]=0}
		}
	}
	
	parent=parent[which(parent$X._REF>=0.9),]
	bp=intersect(temp.count$bp, parent$bp)
	temp.count=temp.count[which(temp.count$bp %in% bp),]
	parent=parent[which(parent$bp %in% bp),]
	
	
	# For each calculated frequency, substract parent frequenc
	temp.count[,15:20]=temp.count[,15:20]-parent[,15:20]
	
	# put negative values to 0
	temp.count[,15:20][temp.count[, 15:20] < 0] = 0
	
	# recalculate counts
	temp.count[,4:9]=temp.count[,15:20]*temp.count$Depth
	
	# simplify table
	temp.count=temp.count[,c(3,4,5,6,7,8,9)]
	
	sum_counts=ddply(temp.count, .(REF), colwise(sum))
	sum_counts.m=melt(sum_counts, measure.vars=2:7)
	
	# remove unchanged bases
	sum_counts.m=sum_counts.m[!as.character(sum_counts.m$REF)==as.character(sum_counts.m$variable),]
	sum_counts.m$percent=sum_counts.m$value/sum(sum_counts.m$value)*100
	
	sum_counts.m$change = apply( sum_counts.m[ , c(1,2)] , 1 , paste , collapse = "->" )

	sum_counts.m$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(sum_counts.m)[1])

	if (n==338 || n==339 || n==345 ||n==346){
		sum_counts.m$recruit=rep('MS2', dim(sum_counts.m)[1])
	}else if ( n==340 ||n==341|| n==347 || n==348){
		sum_counts.m$recruit=rep('C-Fusion', dim(sum_counts.m)[1])
	}
	else if (n==342 ||n==343 || n==349 || n==350){
		sum_counts.m$recruit=rep('N-Fusion', dim(sum_counts.m)[1])
	}
	sum_counts.m$locus=rep(chr,dim(sum_counts.m)[1])
	parents=rbind(parents,parent)
	percent_change=rbind(percent_change,sum_counts.m)
}

write.table( percent_change, 'Percent_change_CX338-351_normalized.txt',col.names=T, row.names=F, quote=F, append=F, sep='\t')

g1=ggplot(percent_change, aes(x=change,y=percent, fill=recruit, group=sample)) +geom_bar(stat="identity",position=position_dodge())+
	labs(x="Base Change", y="Percentage", title="GFP")


g2=ggplot(percent_change, aes(x=change,y=percent, fill=recruit, group=sample)) +geom_bar(stat="identity",position=position_dodge())+
	labs(x="Base Change", y="Percentage", title="HBG2")


