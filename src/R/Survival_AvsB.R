#################################################################
#
# date: Ocotober 22, 2016
# platform: Ubuntu 16.04
# R.version : 3.2.2
# author: Villemin Jean-Philippe
# team: Epigenetic Component of Alternative Splicing - IGH
#
# annotFileWithSymbol
# Usage : 
# 
# 
#
#################################################################

library(optparse)
####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-s", "--survival"), type="character", default=NULL, help="Absolute File Input Path for Survival dataset", metavar="character"),
  make_option(c("-e", "--endPoint"), type="character", default=NULL, help="OS,DSS,DSI,tpm", metavar="character"),
  make_option(c("-a", "--listpatient1"), type="character", default=NULL, help="IDs", metavar="character"),
  make_option(c("-b", "--listpatient2"), type="character", default=NULL, help="IDs", metavar="character"), 
  make_option(c("-o", "--outputName"), type="character", default=NULL, help="IDs", metavar="character")


)



parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

######################################################
#########         MAIN                     ###########
######################################################

####################################################################################
#####################################  Package Loading  ############################
####################################################################################

library(data.table)
library(reshape2)
library(plyr)
library(tidyr)
library(stringr)
require(gridExtra)
library(survminer)
library(survival)

warnings()
time =""

if (opt$endPoint =="OS") {
	
	time = "OS.time"
}
if (opt$endPoint =="DSS") {
	
	time = "DSS.time"
}
if (opt$endPoint =="DFI") {
	
	time = "DFI.time"
}
if (opt$endPoint =="PFI") {
	
	time = "PFI.time"
}

print(time)
print(opt$endPoint)

#print(geneList$V1[2])
patient1 <- read.table(opt$listpatient1, header=FALSE,stringsAsFactors=FALSE)
head(patient1)

patient2 <- read.table(opt$listpatient2, header=FALSE,stringsAsFactors=FALSE)
head(patient2)

tobepasted = c(dirname(opt$listpatient1),"/",opt$outputName,".png")
tobepastedFinal =paste(tobepasted,sep="/",collapse="")
print(tobepastedFinal)
####################################################################################
#####################################  Package Loading  ############################
####################################################################################
###Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV2TPM.R -e OS -m /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/outputBasalTPM.csv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv
#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV2TPM.R -e OS -m /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/FINAL.SF.csv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv

#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV2TPM.R -e OS -m /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/outputBasalTPM.csv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv
#Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/SurvivalV2TPM.R -e $i -m /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/FINAL.SF.csv -s /home/jean-philippe.villemin/data/data/PROJECT/SURVIVAL/TCGA_CDR.csv

survival <- read.table(opt$survival,sep=";", header=TRUE,stringsAsFactors=FALSE, comment.char = "@", na.strings = "#N/A", quote="" )
survival_subset = subset(survival, select = c("bcr_patient_barcode",time,opt$endPoint))
survival_subset_unique_id <- survival_subset[!duplicated(survival_subset$bcr_patient_barcode), ]

df <- data.frame(NAME=character(),PVAL=double(),LOW_CUTOFF=double(),HIGH_CUTOFF=double(),LOWCOUNT=integer(),HIGHCOUNT=integer(),HR=double(),PVALHR=double(),CILOW=double(),CIHIGH=double())

#final <- merge(x = tpm_matrice_clean_unique_id, y = survival_subset_unique_id, by = "bcr_patient_barcode", all.x = TRUE)
final <- survival_subset_unique_id

# Remove NA
final<- final[complete.cases(final[,opt$endPoint]),]

print("wholematrix size")
print(dim(final))
print(head(final))


final$group[final$bcr_patient_barcode %in% patient2$V1]<-"group2"
final$group[final$bcr_patient_barcode %in% patient1$V1]<-"group1"

final <- rbind(final[final$bcr_patient_barcode %in% patient2$V1,], final[final$bcr_patient_barcode %in% patient1$V1,])

final_group_annotated <- final[!is.na(final$group),]

print("group2")
print(dim(final_group_annotated[final_group_annotated$group=="group2",]))
low_count <- nrow(final_group_annotated[final_group_annotated$group=="group2",])

print(dirname(opt$listpatient2))
print(opt$outputName)

finalname = paste(dirname(opt$listpatient2),opt$outputName,sep="/")
finalname2 = paste0(finalname,"_groupEPI.csv")

write.csv(final_group_annotated[final_group_annotated$group=="group2",],row.names=FALSE,file=finalname2)

print("group1")
print(dim(final_group_annotated[final_group_annotated$group=="group1",]))
high_count <- nrow(final_group_annotated[final_group_annotated$group=="group1",])

finalname = paste(dirname(opt$listpatient1),opt$outputName,sep="/")
finalname1 = paste0(finalname,"_groupMES.csv")
write.csv(final_group_annotated[final_group_annotated$group=="group1",],row.names=FALSE,file=finalname1)


print("final_group_annotated")
print(final_group_annotated)

print("total")
dim(final_group_annotated)

#head(final_group_annotated,5)

print("endpoint")
print(final_group_annotated[,opt$endPoint])
print("tpm")
print(final_group_annotated$res_tpm)
print("length")
print(length(final_group_annotated[,time]))

surv_object <- Surv(time = as.numeric(final_group_annotated[,time]), event = final_group_annotated[,opt$endPoint] )
#surv_object
fit1 <- survfit(surv_object ~ group, data = final_group_annotated)

test=surv_pvalue(fit1, final_group_annotated)
print("survpvalue")
print(test)
print(test$pval)
print(surv_pvalue(fit1)[,2])

fit.coxph <- coxph(surv_object ~ group, data = final_group_annotated, ties = 'exact')#
confidenceInterval <- summary(fit.coxph)$conf.int
print("summ Coxph")
print(summary(fit.coxph))

pvc <- summary(fit.coxph)$sctest[3]

#pvc <- coef(summary(fit.coxph))[,5]
hr <- coef(summary(fit.coxph))[,2]
low_confidence <-confidenceInterval[,3]
high_confidence<-confidenceInterval[,4]
#tobepasted = c("HR = ",round(hr,digits=2),"[",round(low_confidence,digits=2),"-",round(high_confidence,digits=2),"]")
tobepasted = c("HR = ",round(1/hr,digits=2),"[",round(1/high_confidence,digits=2),"-",round(1/low_confidence,digits=2),"]")

title <- paste(tobepasted, collapse="",sep="")

print(low_confidence)	
print(high_confidence)
print(hr)
print(pvc)

print(strsplit(opt$outputName, "_")[[1]][3])
df <- data.frame(ID=character(),NAME=character(),PVAL=double(),HR=double(),PVALHR=double(),CILOW=double(),CIHIGH=double())

df <- rbind(df, data.frame(ID=strsplit(opt$outputName, "_")[[1]][1],NAME = strsplit(opt$outputName, "_")[[1]][2], PVAL = test$pval[1],HR=hr,PVALHR=pvc,CILOW=low_confidence,CIHIGH=high_confidence))

tobepasted = c(dirname(opt$listpatient1),"/","PVAL","_",strsplit(opt$outputName, "_")[[1]][3],".csv")
towrite =paste(tobepasted,sep="/",collapse="")
out = file(towrite, 'a')
write.table(df,col.names=FALSE,file=out,quote=FALSE,sep = ",")

mypalette = c("A-LIKE" = "#00BFC4", "B-LIKE"="#F8766D")
#scale_fill_manual("legend", values = c("BASALB-LIKE" = "#F8766D", "BASALA-LIKE" = "#00BFC4")
final_group_annotated$group[final_group_annotated$group=="group1"]<- "B-LIKE"
final_group_annotated$group[final_group_annotated$group=="group2"]<- "A-LIKE"

ggsurv <- ggsurvplot(fit1, data = final_group_annotated, legend.title  = "Group :", pval = TRUE,xlab="Years",title=title,break.x.by =365.25,xscale=365.25,xlim=c(0,2000), legend.labs = c("B-LIKE" , "A-LIKE" ),palette=mypalette)
png(tobepastedFinal,width=10,height=10,units = 'cm', res = 300)
print({	ggpubr::ggpar(ggsurv,
							legend = "top",
							font.legend = list(size = 12, color = "black", face = "bold"),
							font.x = c(12 , "bold", "black"),           # font for x axises in the plot, the table and censor part
							font.y = c(12, "bold", "black"),       # font for y axises in the plot, the table and censor part
							font.tickslab = c(14, "plain", "black") 
					) })

dev.off()

