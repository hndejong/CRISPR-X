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

setwd('/users/lfresard/CRISPR-X/data/mpileup/bp_transition_GFP')

counts_file<-list.files(pattern="count")
temp.count<-read.table(counts_file[1], header=T)

chr="GFP684"


####Filter FOR HOTSPOT
big_enrichment=c()
big_enrichment2=c()
parents=c()

for (i in 1:length(counts_file)){
	temp.count<-read.table(counts_file[i], header=T)
	temp.count=temp.count[temp.count$CHR==chr,]

	temp.count$A_freq=temp.count$A/temp.count$Depth
	temp.count$G_freq=temp.count$G/temp.count$Depth
	temp.count$C_freq=temp.count$C/temp.count$Depth
	temp.count$T_freq=temp.count$T/temp.count$Depth
	temp.count$del_freq=temp.count$del/temp.count$Depth
	temp.count$ins_freq=temp.count$ins/temp.count$Depth

	if (unlist(strsplit(counts_file[i],"_"))[1]=="CX6"||unlist(strsplit(counts_file[i],"_"))[1]=="CX26"||unlist(strsplit(counts_file[i],"_"))[1]=="CX27"||unlist(strsplit(counts_file[i],"_"))[1]=="CX28"||unlist(strsplit(counts_file[i],"_"))[1]=="CX29"||unlist(strsplit(counts_file[i],"_"))[1]=="CX30"||unlist(strsplit(counts_file[i],"_"))[1]=="CX31"||unlist(strsplit(counts_file[i],"_"))[1]=="CX32"||unlist(strsplit(counts_file[i],"_"))[1]=="CX33"||unlist(strsplit(counts_file[i],"_"))[1]=="CX34"){
		parent<-read.table("./parents/CX9_n5_mapq30_sorted_qual30.count", header=T)
		parent=parent[parent$CHR==chr,]
		parent$A_freq=parent$A/parent$Depth
		parent$G_freq=parent$G/parent$Depth
		parent$C_freq=parent$C/parent$Depth
		parent$T_freq=parent$T/parent$Depth
		parent$del_freq=parent$del/parent$Depth
		parent$ins_freq=parent$ins/parent$Depth

		if (unlist(strsplit(counts_file[i],"_"))[1]=="CX6"){
			temp.count=temp.count[temp.count$bp >=268 & temp.count$bp <=288,]
			parent=parent[parent$bp>=268 & parent$bp<=288,]

		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX26"||unlist(strsplit(counts_file[i],"_"))[1]=="CX27" || unlist(strsplit(counts_file[i],"_"))[1]=="CX28"){
			temp.count=temp.count[temp.count$bp >=295 & temp.count$bp <=315,]
			parent=parent[parent$bp>=295 & parent$bp<=315,]
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX29"||unlist(strsplit(counts_file[i],"_"))[1]=="CX30" || unlist(strsplit(counts_file[i],"_"))[1]=="CX31"){
			temp.count=temp.count[temp.count$bp >=546 & temp.count$bp <=566,]
			parent=parent[parent$bp>=546 & parent$bp<=566,]
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX32"||unlist(strsplit(counts_file[i],"_"))[1]=="CX33" || unlist(strsplit(counts_file[i],"_"))[1]=="CX34"){
			temp.count=temp.count[temp.count$bp >=512 & temp.count$bp <=532,]
			parent=parent[parent$bp>=512 & parent$bp<=532,]
		}

	}
	else if (unlist(strsplit(counts_file[i],"_"))[1]=="CX10" || unlist(strsplit(counts_file[i],"_"))[1]=="CX11"){
		temp.count=temp.count[temp.count$bp >=268 & temp.count$bp <=288,]
		parent<-read.table("./parents/CX17_n3_mapq30_sorted_qual30.count", header=T)
		parent=parent[parent$CHR==chr,]
		parent$A_freq=parent$A/parent$Depth
		parent$G_freq=parent$G/parent$Depth
		parent$C_freq=parent$C/parent$Depth
		parent$T_freq=parent$T/parent$Depth
		parent$del_freq=parent$del/parent$Depth
		parent$ins_freq=parent$ins/parent$Depth
		parent=parent[parent$bp>=268 & parent$bp<=288,]
	}
	else if ( unlist(strsplit(counts_file[i],"_"))[1] %in% c("CX42","CX43", "CX44","CX45","CX46","CX47","CX48", "CX49", "CX50","CX51","CX52","CX53","CX54","CX55","CX56","CX57","CX58","CX59","CX60","CX61","CX62","CX63","CX64","CX65")){
		parent<-read.table("./parents/CX66_n5_mapq30_sorted_qual30.count", header=T)
		parent=parent[parent$CHR==chr,]
		parent$A_freq=parent$A/parent$Depth
		parent$G_freq=parent$G/parent$Depth
		parent$C_freq=parent$C/parent$Depth
		parent$T_freq=parent$T/parent$Depth
		parent$del_freq=parent$del/parent$Depth
		parent$ins_freq=parent$ins/parent$Depth
		if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX42"||unlist(strsplit(counts_file[i],"_"))[1]=="CX43" || unlist(strsplit(counts_file[i],"_"))[1]=="CX44"){
			temp.count=temp.count[temp.count$bp >=349 & temp.count$bp <=369,]
			parent=parent[parent$bp>=349 & parent$bp<=369,]
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX45"||unlist(strsplit(counts_file[i],"_"))[1]=="CX46" || unlist(strsplit(counts_file[i],"_"))[1]=="CX47"){
			temp.count=temp.count[temp.count$bp >=246 & temp.count$bp <=266,]
			parent=parent[parent$bp>=246 & parent$bp<=266,]
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX48"||unlist(strsplit(counts_file[i],"_"))[1]=="CX49" || unlist(strsplit(counts_file[i],"_"))[1]=="CX50"){
			temp.count=temp.count[temp.count$bp >=483 & temp.count$bp <=503,]
			parent=parent[parent$bp>=483 & parent$bp<=503,]
	
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX51"||unlist(strsplit(counts_file[i],"_"))[1]=="CX52" || unlist(strsplit(counts_file[i],"_"))[1]=="CX53"){
			temp.count=temp.count[temp.count$bp >=528 & temp.count$bp <=548,]
			parent=parent[parent$bp>=528 & parent$bp<=548,]
		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX54"||unlist(strsplit(counts_file[i],"_"))[1]=="CX55" || unlist(strsplit(counts_file[i],"_"))[1]=="CX56"){
			temp.count=temp.count[temp.count$bp >=590 & temp.count$bp <=610,]
			parent=parent[parent$bp>=590 & parent$bp<=610,]
		}
	
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX57"||unlist(strsplit(counts_file[i],"_"))[1]=="CX58" || unlist(strsplit(counts_file[i],"_"))[1]=="CX59"){
			temp.count=temp.count[temp.count$bp >=191 & temp.count$bp <=211,]
			parent=parent[parent$bp>=191 & parent$bp<=211,]
		}	

		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX60"||unlist(strsplit(counts_file[i],"_"))[1]=="CX61" || unlist(strsplit(counts_file[i],"_"))[1]=="CX62"){
			temp.count=temp.count[temp.count$bp >=338 & temp.count$bp <=358,]
			parent=parent[parent$bp>=338 & parent$bp<=358,]

		}
		else if ( unlist(strsplit(counts_file[i],"_"))[1]=="CX63"||unlist(strsplit(counts_file[i],"_"))[1]=="CX64" || unlist(strsplit(counts_file[i],"_"))[1]=="CX65"){
			temp.count=temp.count[temp.count$bp >=503 & temp.count$bp <=523,]
			parent=parent[parent$bp>=503 & parent$bp<=523,]

		}
	
	}
									
	temp.count=temp.count[,c('bp', 'REF', 'A_freq', 'G_freq', 'C_freq' ,'T_freq','del_freq' ,'ins_freq','X._REF')]
	for (j in 1:(dim(temp.count)[1])){
		for (k in 3:8){
			if (abs(temp.count[j,k]-temp.count[j,9])<10e-10){temp.count[j,k]=0}
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

	#big_enrichment=rbind(big_enrichment,enrichment)
	big_enrichment2=rbind(big_enrichment2,enrichment2)

	parents=rbind(parents,parent)

	#names=c(names,unlist(strsplit(counts_file[i]), "_"))[1])
	#res[,i]=
}

#big_enrichment.m=melt(big_enrichment,measure.vars=3:8)
big_enrichment2.m=melt(big_enrichment2,measure.vars=3:8)




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



pdf('GFP_bpenrichment_20bphotspot_minusparent_filtered.pdf')
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

#pdf('GFP_bpenrichment_20bphotspot_minusparent_negativetozero.pdf')
#g
#dev.off()




pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_scale1.pdf')
g+scale_y_continuous(limits=c(0,.25))
dev.off()

pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_points_scale1.pdf')
h+scale_y_continuous(limits=c(0,.25))
dev.off()


pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_scale2.pdf')
g+scale_y_continuous(limits=c(0,.025))
dev.off()

pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_points_scale2.pdf')
h+scale_y_continuous(limits=c(0,.025))
dev.off()

pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_scale3.pdf')
g+scale_y_continuous(limits=c(0,.025))
dev.off()

pdf('GFP_bpenrichment_20bphotspotminusparent_negativetozero_rescaled_points_scale3.pdf')
h+scale_y_continuous(limits=c(0,.0025))
dev.off()



#big_df.m=melt(big_df,measure.vars=3:8)
#big_df.m$value=big_df.m$value+10e-20
write.table(big_enrichment2.m, "TransitionTableGFP_data_20bphotspot_minusparent-filtered.txt",col.names=T, sep="\t", quote=F, append =F, row.names=F)

cdata<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,median=median(value))

write.table(cdata,"MedianTransitions_GFP_20bp_minusparentdata-filtered.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)

cdata_max<-ddply(big_enrichment2.m, c("REF", "variable"),summarise,max=max(value))
write.table(cdata_max,"MaxTransitions_GFP_20bp_minusparentdata-filtered.txt", col.names=T, sep="\t", quote=F, append=F, row.names=F)

big_enrichment2.m$value[which(big_enrichment2.m$value <=0)]=min(big_enrichment2.m$value[which(big_enrichment2.m$value >0)])

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
g<-g+scale_y_log10()

pdf('GFP_bpenrichment_20bphotspot_minusparent_logscale-filetered.pdf')
g
dev.off()
scale_y_log10(limits=c(5e-9, 0.25))

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
pdf('GFP_bpenrichment_20bphotspot_minusparent_filtered.pdf')
g
dev.off()

