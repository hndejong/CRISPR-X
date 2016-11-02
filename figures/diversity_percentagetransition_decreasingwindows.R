#!/bin/R
#LF
#durga

library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(plyr)
##For each sample have a decreasing window and calculate the frequency of transition that pass noise

palette=brewer.pal(9,"Set1")
mycol_DNA=c(palette[2],palette[3],palette[1],brewer.pal(8,"Set2")[6],palette[4],palette[9])
#[1] "A_freq"   "G_freq"   "C_freq"   "T_freq"   "del_freq" "ins_freq"

setwd('/users/lfresard/CRISPR-X/data/mpileup/upmutator')

counts_file<-list.files(pattern="count")

####Filter FOR HOTSPOT 100bp
big_enrichment2=list()
parents=list()

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
		S=485
		E=585
	}
	else if (n >=182 && n<=184){
		S=130
		E=230

	}
	else if (n >=211 && n<=212){
		S=275
		E=375
	}
	else if (n >=213 && n<=214){
		S=278
		E=378
	}
	else if (n >=215 && n<=216){
		S=318
		E=418
	}
	else if (n >=217 && n<=218){
		S=147
		E=247

	}
	else if (n >=219 && n<=220){
		S=443
		E=543
	}
	else if (n >=222 && n<=223){
		S=269
		E=369
	}
	else if (n >=224 && n<=225){
		S=263
		E=363
	}
	else if (n >=226 && n<=227){
		S=374
		E=474
	}
	else if (n >=228 && n<=229){
		S=405
		E=505
	}
	else if (n >=230 && n<=231){
		S=179
		E=279
	}
	else if (n >=233 && n<=234){
		S=499
		E=599
	}
	else if (n >=235 && n<=236){
		S=512
		E=612
	}
	else if (n >=237 && n<=238){
		S=582
		E=682
	}
	else if (n >=239 && n<=240){
		S=265
		E=365
	}
	else if (n >=241 && n<=242){
		S=242
		E=342
	}
	else if (n >=244 && n<=245){
		S=820
		E=920
	}
	else if (n >=246 && n<=247){
		S=222
		E=322
	}
	else if (n >=248 && n<=249){
		S=127
		E=227
	}
	else if (n >=251 && n<=252){
		S=768
		E=868
	}	
	else if (n >=253 && n<=254){
		S=704
		E=804
	}
	else if (n >=255 && n<=256){
		S=242
		E=342
	}
	else if (n >=257 && n<=258){
		S=224
		E=324
	}
	else if (n >=260 && n<=261){
		S=154
		E=254
	}
	else if (n >=262 && n<=263){
		S=1504
		E=1604
	}
	else if (n >=270 && n<=271){
		S=556
		E=656
	}
	else if (n >=272 && n<=273){
		S=279
		E=379
	}
	else if (n >=274 && n<=275){
		S=1043
		E=1143
	}
	else if (n >=277 && n<=278){
		S=221
		E=321
	}
	else if (n >=279 && n<=280){
		S=225
		E=325
	}
	else if (n >=281 && n<=282){
		S=483
		E=583
	}
	else if (n >=283 && n<=284){
		S=572
		E=672
	}
	else if (n >=285 && n<=286){
		S=1089
		E=1189
	}
###### For each window, calculate the frequency of possible observations#######	
	windows_dim=c()
	for (z in 1:length(seq(S,E,10)[1:5])){
		seq1=seq(S,E,10)[1:5]
		seq2=seq(S,E,10)[11:7]
		windows_dim[z]=length(seq(seq1[z],seq2[z]))
		tmp=temp.count[temp.count$bp >=seq1[z] & temp.count$bp <=seq2[z],]
		par=parent[parent$bp >=seq1[z] & parent$bp <=seq2[z],]


		tmp=temp.count[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
		for (j in 1:(dim(temp.count)[1])){
			for (k in 3:8){
				if (abs(tmp[j,k]-tmp[j,9])<10e-12){tmp[j,k]=0}
			}
		}
		par=par[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
		for (j in 1:(dim(par)[1])){
			for (k in 3:8){
				if (abs(par[j,k]-par[j,9])<10e-12){par[j,k]=0}
			}
		}

		par=par[which(par$X._REF>=0.9),]
		bp=intersect(tmp$bp, par$bp)

		tmp=tmp[which(tmp$bp %in% bp),]
		par=par[which(par$bp %in% bp),]


		enrichment2=cbind(tmp[,1:2],tmp[,3:8]-par[,3:8])

		enrichment2=as.data.frame(enrichment2)
		enrichment2$sample=rep(unlist(strsplit(counts_file[i],"_"))[1],dim(enrichment2)[1])
		parents[[as.character(z)]]=rbind(parents[[as.character(z)]],par)

		big_enrichment2[[as.character(z)]]=rbind(big_enrichment2[[as.character(z)]],enrichment2)

	}


}

#parents.m=melt(parents,measure.vars=3:8)
#big_enrichment.m=melt(big_enrichment,measure.vars=3:8)
res=list()
for (i in 1:length(big_enrichment2)){
	temp.m=melt(big_enrichment2[[i]],measure.vars=3:8)
	samples=unique(temp.m$sample)

	res[[as.character(i)]]=array(dim=c(length(samples),2))
	for (j in 1: length(samples)){
	temp=temp.m[which(temp.m$sample==samples[j]),]
	toto=length(which(temp$value>=sd(1-parents[[i]]$X._REF)))
	res[[as.character(i)]][j,1]=samples[j]
	res[[as.character(i)]][j,2]=toto/(windows_dim[i]*3)*100
	}
}


medians=c()

for (i in 1:length(res)){
	medians[i]=median(as.numeric(res[[i]][,2]))
	
}
res_med=data.frame(WINDOW=windows_dim, MEDIAN_PCT=medians)

png('ImpactWS_Percentage_combinations_upmutator.png')
barplot(res_med$MEDIAN_PCT, names.arg=as.character(res_med$WINDOW), border=F, ylim=c(0,20), xlab='Window Size (bp)', ylab='Percentage combinations over noise')
dev.off()

write.table(res_med,'ImpactWS_Percentage_combinations_upmutator.txt', col.names=T, row.names=F, quote=F,append=F)



