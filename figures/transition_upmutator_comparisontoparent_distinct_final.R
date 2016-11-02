#!/bin/R
#LF
#durga

library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(plyr)

palette=brewer.pal(9,"Set1")
mycol_DNA=c(palette[2],palette[3],palette[1],brewer.pal(8,"Set2")[6],palette[4],palette[9])
#[1] "A_freq"   "G_freq"   "C_freq"   "T_freq"   "del_freq" "ins_freq"

setwd('/users/lfresard/CRISPR-X/data/mpileup/upmutator')

counts_file<-list.files(pattern="count")

####Filter FOR HOTSPOT 100bp
#big_df=c()
big_enrichment=c()
big_enrichment2=c()
parents=c()

for (i in 1:length(counts_file)){
#for (i in 1:1){
	temp.count<-read.table(counts_file[i], header=T)
	#filter for different possible chromosomes
		#get sample id to filter on it
	n<-as.numeric(gsub("CX","", unlist(strsplit(counts_file[i],"_"))[1]))
	if (n >=179 && n<=184){
		chr="GFP684"
		parent=read.table("./parents/CX188_n5_mapq30_sorted_qual30.count", header=T)
	}
	else if( n>=211 && n<=220){
		chr="HBG2"
		parent=read.table("./parents/CX221_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n>=222 && n <=231){
		chr="GSTP1"
		parent=read.table("./parents/CX232_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=233 && n<=242){
		chr="FTL"
		parent=read.table("./parents/CX243_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=244 && n<=249){
		chr="CD45_Prom1"
		parent=read.table("./parents/CX250_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=251 && n<=258){
		chr="CD45_Prom2"
		parent=read.table("./parents/CX259_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=260 && n<=263){
		chr="CD274_Prom_SNP12"
		parent=read.table("./parents/CX264_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=270 && n<=275){
		chr="CD274_Prom"
		parent=read.table("./parents/CX276_n3_mapq30_sorted_qual30.count", header=T)
	}
	else if (n >=277 && n<=286){
		chr="CD14_Prom"
		parent=read.table("./parents/CX289_n3_mapq30_sorted_qual30.count", header=T)
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

	if (n >=179 && n<=181){
		temp.count=temp.count[temp.count$bp >=485 & temp.count$bp <=585,]
		parent=parent[parent$bp >=485 & parent$bp <=585,]

	}
	else if (n >=182 && n<=184){
		temp.count=temp.count[temp.count$bp >=130 & temp.count$bp <=230,]
		parent=parent[parent$bp >=130 & parent$bp <=230,]

	}
	else if (n >=211 && n<=212){
		temp.count=temp.count[temp.count$bp >=275 & temp.count$bp <=375,]
		parent=parent[parent$bp >=275 & parent$bp <=375,]
	}
	else if (n >=213 && n<=214){
		temp.count=temp.count[temp.count$bp >=278 & temp.count$bp <=378,]
		parent=parent[parent$bp >=278 & parent$bp <=378,]
	}
	else if (n >=215 && n<=216){
		temp.count=temp.count[temp.count$bp >=318 & temp.count$bp <=418,]
		parent=parent[parent$bp >=318 & parent$bp <=418,]
	}
	else if (n >=217 && n<=218){
		temp.count=temp.count[temp.count$bp >=147 & temp.count$bp <=247,]
		parent=parent[parent$bp >=147 & parent$bp <=247,]
	}
	else if (n >=219 && n<=220){
		temp.count=temp.count[temp.count$bp >=443 & temp.count$bp <=543,]
		parent=parent[parent$bp >=443 & parent$bp <=543,]
	}
	else if (n >=222 && n<=223){
		temp.count=temp.count[temp.count$bp >=269 & temp.count$bp <=369,]
		parent=parent[parent$bp >=269 & parent$bp <=369,]
	}
	else if (n >=224 && n<=225){
		temp.count=temp.count[temp.count$bp >=263 & temp.count$bp <=363,]
		parent=parent[parent$bp >=263 & parent$bp <=363,]
	}
	else if (n >=226 && n<=227){
		temp.count=temp.count[temp.count$bp >=374 & temp.count$bp <=474,]
		parent=parent[parent$bp >=374 & parent$bp <=474,]
	}
	else if (n >=228 && n<=229){
		temp.count=temp.count[temp.count$bp >=405 & temp.count$bp <=505,]
		parent=parent[parent$bp >=405 & parent$bp <=505,]
	}
	else if (n >=230 && n<=231){
		temp.count=temp.count[temp.count$bp >=179 & temp.count$bp <=279,]
		parent=parent[parent$bp >=179 & parent$bp <=279,]
	}
	else if (n >=233 && n<=234){
		temp.count=temp.count[temp.count$bp >=398 & temp.count$bp <=498,]
		parent=parent[parent$bp >=398 & parent$bp <=498,]
	}
	else if (n >=235 && n<=236){
		temp.count=temp.count[temp.count$bp >=512 & temp.count$bp <=612,]
		parent=parent[parent$bp >=512 & parent$bp <=612,]
	}
	else if (n >=237 && n<=238){
		temp.count=temp.count[temp.count$bp >=582 & temp.count$bp <=682,]
		parent=parent[parent$bp >=582 & parent$bp <=682,]
	}
	else if (n >=239 && n<=240){
		temp.count=temp.count[temp.count$bp >=265 & temp.count$bp <=365,]
		parent=parent[parent$bp >=265 & parent$bp <=365,]
	}
	else if (n >=241 && n<=242){
		temp.count=temp.count[temp.count$bp >=242 & temp.count$bp <=342,]
		parent=parent[parent$bp >=242 & parent$bp <=342,]
	}
	else if (n >=244 && n<=245){
		temp.count=temp.count[temp.count$bp >=820 & temp.count$bp <=920,]
		parent=parent[parent$bp >=820 & parent$bp <=920,]
	}
	else if (n >=246 && n<=247){
		temp.count=temp.count[temp.count$bp >=222 & temp.count$bp <=322,]
		parent=parent[parent$bp >=222 & parent$bp <=322,]
	}
	else if (n >=248 && n<=249){
		temp.count=temp.count[temp.count$bp >=127 & temp.count$bp <=227,]
		parent=parent[parent$bp >=127 & parent$bp <=227,]
	}
	else if (n >=251 && n<=252){
		temp.count=temp.count[temp.count$bp >=768 & temp.count$bp <=868,]
		parent=parent[parent$bp >=768 & parent$bp <=868,]
	}	
	else if (n >=253 && n<=254){
		temp.count=temp.count[temp.count$bp >=704 & temp.count$bp <=804,]
		parent=parent[parent$bp >=704 & parent$bp <=804,]
	}
	else if (n >=255 && n<=256){
		temp.count=temp.count[temp.count$bp >=242 & temp.count$bp <=342,]
		parent=parent[parent$bp >=242 & parent$bp <=342,]
	}
	else if (n >=257 && n<=258){
		temp.count=temp.count[temp.count$bp >=224 & temp.count$bp <=324,]
		parent=parent[parent$bp >=224 & parent$bp <=324,]
	}
	else if (n >=260 && n<=261){
		temp.count=temp.count[temp.count$bp >=154 & temp.count$bp <=254,]
		parent=parent[parent$bp >=154 & parent$bp <=254,]
	}
	else if (n >=262 && n<=263){
		temp.count=temp.count[temp.count$bp >=1504 & temp.count$bp <=1604,]
		parent=parent[parent$bp >=1504 & parent$bp <=1604,]
	}
	else if (n >=270 && n<=271){
		temp.count=temp.count[temp.count$bp >=556 & temp.count$bp <=656,]
		parent=parent[parent$bp >=556 & parent$bp <=656,]
	}
	else if (n >=272 && n<=273){
		temp.count=temp.count[temp.count$bp >=279 & temp.count$bp <=379,]
		parent=parent[parent$bp >=279 & parent$bp <=379,]
	}
	else if (n >=274 && n<=275){
		temp.count=temp.count[temp.count$bp >=1043 & temp.count$bp <1143,]
		parent=parent[parent$bp >=1043 & parent$bp <=1143,]
	}
	else if (n >=277 && n<=278){
		temp.count=temp.count[temp.count$bp >=221 & temp.count$bp <=321,]
		parent=parent[parent$bp >=221 & parent$bp <=321,]
	}
	else if (n >=279 && n<=280){
		temp.count=temp.count[temp.count$bp >=225 & temp.count$bp <=325,]
		parent=parent[parent$bp >=225 & parent$bp <=325,]
	}
	else if (n >=281 && n<=282){
		temp.count=temp.count[temp.count$bp >=483 & temp.count$bp <=583,]
		parent=parent[parent$bp >=483 & parent$bp <=583,]
	}
	else if (n >=283 && n<=284){
		temp.count=temp.count[temp.count$bp >=572 & temp.count$bp <=672,]
		parent=parent[parent$bp >=572 & parent$bp <=672,]
	}
	else if (n >=285 && n<=286){
		temp.count=temp.count[temp.count$bp >=1089 & temp.count$bp <=1189,]
		parent=parent[parent$bp >=1089 & parent$bp <=1189,]
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
	#enrichment=cbind(temp.count[,1:2],(temp.count[,3:8]+10e-12)/(parent[,3:8]+10e-12))

	#enrichment=as.data.frame(enrichment)
	enrichment2=as.data.frame(enrichment2)
	#enrichment$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(enrichment)[1])
	enrichment2$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(enrichment2)[1])
	parents=rbind(parents,parent)

	#big_enrichment=rbind(big_enrichment,enrichment)
	big_enrichment2=rbind(big_enrichment2,enrichment2)

	#temp.count=as.data.frame(temp.count)
	#temp.count$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(temp.count)[1])

	#big_df=rbind(big_df,temp.count)
	#names=c(names,unlist(strsplit(counts_file[i]), "_"))[1])
	#res[,i]=
}
parents.m=melt(parents,measure.vars=3:8)
#big_enrichment.m=melt(big_enrichment,measure.vars=3:8)
big_enrichment2.m=melt(big_enrichment2,measure.vars=3:8)











big_enrichment2.m$value[which(big_enrichment2.m$value <=0)]=0

g<-ggplot(big_enrichment2.m, aes(x=REF,y=value, colour=variable))
g<-g+geom_boxplot(fill="grey",outlier.size=1.5)
g<-g+facet_grid(.~REF, scales="free")
g<-g+ scale_colour_manual(values=mycol_DNA,name="Transition/Transversion",labels=c("A", "G","C","T","deletion","insertion"))
g<-g+labs(x='Reference', y='Frequency of Alternative Alleles')+ theme_bw()
g<-g+ theme(axis.title.x=element_text(size=15)) 
g<-g+theme(axis.title.y=element_text(size=15))+ theme(legend.text=element_text(size=12))
g<-g+ theme(legend.title=element_text(size=15))
g<-g+ theme(axis.text.x= element_text(size=10))
g<-g+ theme(axis.text.y= element_text(size=10))
g<-g+theme(axis.text.x=element_blank())



h<-ggplot(big_enrichment2.m, aes(x=REF,y=value, colour=variable))
h<-h+geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=I(1/3))
h<-h+facet_grid(.~REF, scales="free")
h<-h+ scale_colour_manual(values=mycol_DNA,name="Transition/Transversion",labels=c("A", "G","C","T","deletion","insertion"))
h<-h+labs(x='Reference', y='Frequency of Alternative Alleles')+ theme_bw()
h<-h+ theme(axis.title.x=element_text(size=15)) 
h<-h+theme(axis.title.y=element_text(size=15))+ theme(legend.text=element_text(size=12))
h<-h+ theme(legend.title=element_text(size=15))
h<-h+ theme(axis.text.x= element_text(size=10))
h<-h+ theme(axis.text.y= element_text(size=10))
h<-h+theme(axis.text.x=element_blank())
h


pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_scale1_ok.pdf')
g+scale_y_continuous(limits=c(0,.25))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_points_scale1_ok.pdf')
h+scale_y_continuous(limits=c(0,.25))
dev.off()


pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_scale2_ok.pdf')
g+scale_y_continuous(limits=c(0,.025))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_points_scale2_ok.pdf')
h+scale_y_continuous(limits=c(0,.025))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_scale3_ok.pdf')
g+scale_y_continuous(limits=c(0,.0025))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled_points_scale3_ok.pdf')
h+scale_y_continuous(limits=c(0,.0025))
dev.off()


write.table(big_enrichment2.m, "TransitionTableupmutator_data_100bphotspot_minusparent-filtered_ok.txt",col.names=T, sep="\t", quote=F, append =F, row.names=F)

cdata<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,median=median(value))

write.table(cdata,"MedianTransitions_upmutator100bp_minusparentdata-filtered_ok.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)

cdata_max<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,max=max(value))
write.table(cdata_max,"MaxTransitions_upmutator_100bp_minusparentdata-filtered_ok.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)


#write.table(big_enrichment2.m, "TransitionTableupmutator_data_1000bphotspot_minusparent-filtered.txt",col.names=T, sep="\t", quote=F, append =F, row.names=F)
#cdata<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,median=median(value))
#write.table(cdata,"MedianTransitions_upmutator_100bp_minusparentdata-filtered.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)
#cdata_max<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,max=max(value))
#write.table(cdata_max,"MaxTransitions_upmutator_100bp_minusparentdata-filtered.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)

mean(1-parents$X._REF)
#[1] 0.000388935
sd(1-parents$X._REF)
#[1] 0.0003736355




g<-ggplot(big_enrichment2.m, aes(x=REF,y=value, colour=variable))
g<-g+geom_boxplot(fill="grey",outlier.size=1.5)
g<-g+facet_grid(.~REF, scales="free")
g<-g+ scale_colour_manual(values=mycol_DNA,name="Transition/Transversion",labels=c("A", "G","C","T","deletion","insertion"))
g<-g+labs(x='Reference', y='Frequency of Alternative Alleles')+ theme_bw()
g<-g+ theme(axis.title.x=element_text(size=15)) 
g<-g+theme(axis.title.y=element_text(size=15))+ theme(legend.text=element_text(size=12))
g<-g+ theme(legend.title=element_text(size=15))
g<-g+ theme(axis.text.x= element_text(size=10))
g<-g+ theme(axis.text.y= element_text(size=10))
g<-g+theme(axis.text.x=element_blank())

pdf('upmutator_bpenrichment_100bphotspot_minusparent_filtered.pdf')
g
dev.off()

big_enrichment2.m$value[which(big_enrichment2.m$value <=0)]=0

g<-ggplot(big_enrichment2.m, aes(x=REF,y=value, colour=variable))
g<-g+geom_boxplot(fill="grey",outlier.size=1.5)
g<-g+facet_grid(.~REF, scales="free")
g<-g+ scale_colour_manual(values=mycol_DNA,name="Transition/Transversion",labels=c("A", "G","C","T","deletion","insertion"))
g<-g+labs(x='Reference', y='Frequency of Alternative Alleles')+ theme_bw()
g<-g+ theme(axis.title.x=element_text(size=15)) 
g<-g+theme(axis.title.y=element_text(size=15))+ theme(legend.text=element_text(size=12))
g<-g+ theme(legend.title=element_text(size=15))
g<-g+ theme(axis.text.x= element_text(size=10))
g<-g+ theme(axis.text.y= element_text(size=10))
g<-g+theme(axis.text.x=element_blank())

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero.pdf')
g+scale_y_continuous(limits=c(0,.22))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled.pdf')
g+scale_y_continuous(limits=c(0,.05))
dev.off()

pdf('upmutator_bpenrichment_100bphotspot_minusparent_negativetozero_rescaled2.pdf')
g+scale_y_continuous(limits=c(0,.0025))
dev.off()


#big_df.m=melt(big_df,measure.vars=3:8)
#write.table(big_df.m, "TransitionTableupmutator_100bp.txt",col.names=T, sep="\t", quote=F, append =F, row.names=F)

#get median values
cdata<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,median=median(value))
write.table(cdata,"MedianTransitions_upmutator_100bp_minusparentdata.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)

cdata_max<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,max=max(value))
write.table(cdata_max,"MaxTransitions_upmutator_100bp_minusparentdata.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)


res=array(dim=c(length(samples),2))
samples=unique(big_enrichment2.m$sample)
for (i in 1: length(samples)){
	temp=big_enrichment2.m[which(big_enrichment2.m$sample==samples[i]),]
	toto=length(which(temp$value>=sd(1-parents$X._REF)))
	res[i,1]=samples[i]
	res[i,2]=toto/303*100
}
colnames(res)=c('sample','percent_comb_over_noise')

write.table(res,'diversitycombination_persample.txt', col.names=T, row.names=F, quote=F,append=F)


