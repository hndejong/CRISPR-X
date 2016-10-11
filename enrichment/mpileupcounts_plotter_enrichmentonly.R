library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)



###Handle arguments
args <- commandArgs(trailingOnly=T) 
mut=as.character(args[1])
cont=as.character(args[2])
chr=as.character(args[3])
dir=as.character(args[4])
###set working directory
setwd(dir)


###Load data
mutated<-read.table(mut, sep="\t", header =T)
mutated=as.data.frame(mutated)


mutated$A_freq=mutated$A/mutated$Depth
mutated$G_freq=mutated$G/mutated$Depth
mutated$C_freq=mutated$C/mutated$Depth
mutated$T_freq=mutated$T/mutated$Depth
mutated$del_freq=mutated$del/mutated$Depth
mutated$ins_freq=mutated$ins/mutated$Depth

control=read.table(cont, sep="\t", header=T)
control=as.data.frame(control)


control$A_freq=control$A/control$Depth
control$G_freq=control$G/control$Depth
control$C_freq=control$C/control$Depth
control$T_freq=control$T/control$Depth
control$del_freq=control$del/control$Depth
control$ins_freq=control$ins/control$Depth



#Filter on the right chromosome
mutated=mutated[which(mutated$CHR==chr),]
control=control[which(control$CHR==chr),]


#Get common bp between control and mutated
bp_com=intersect(mutated$bp, control$bp)

mutated=mutated[mutated$bp %in% bp_com ,]
control=control[control$bp %in% bp_com,]

#Create summary table with alternative allele frequency

comp=data.frame(BP=mutated$bp,mutated=(1-mutated$X._REF), control=(1-control$X._REF))


pdf(paste(substr(mut,0,5),'_', substr(cont,0,5),substr(cont,23,30),'_',chr,'_ratio.pdf', sep=''))
plot((comp$mutated/comp$control)~comp$BP, xlab="Position", ylab="Frequency of alternative alleles: Mutated/Control", main=paste(substr(mut,0,5), " vs ", substr(cont,0,5), sep=''), pch =20,cex=.5)
dev.off()






