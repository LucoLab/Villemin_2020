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

library(data.table)
library(reshape2)

library(plyr)
library(reshape)
library(ggplot2)

library(tidyr)
library(stringr)
require(gridExtra)
library(optparse)
library(ggpubr)
library(dplyr)
require(gdata)

library(ggsignif)
####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-m", "--matricetpm"), type="character", default=NULL, help="Absolute File Input Path for PSI", metavar="character"),
  make_option(c("-a", "--listpatient1"), type="character", default=NULL, help="IDs", metavar="character"),
  make_option(c("-b", "--listpatient2"), type="character", default=NULL, help="IDs", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, help="Features ", metavar="character") ,
  make_option(c("-o", "--outputName"), type="character", default=NULL, help="IDs", metavar="character")

)

#Rscript /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/tpm_per_group.R  -f /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_01_17-0850_0_1000_0_BASALB_BASALA__BASAL_4/symbol_features_final.txt -m /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/Analysis/BasalB_vs_BasalA_CCLE_GSE73526_PRJNA210428_TCGA_BASAL/MatriceTPM_basal.csv -a /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_01_17-0850_0_1000_0_BASALB_BASALA__BASAL_4/10_Patients_MES.txt -b /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_01_17-0850_0_1000_0_BASALB_BASALA__BASAL_4/10_Patients_EPI.txt -o ExpresionSplicing


parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

patient1 <- read.table(opt$listpatient1, header=FALSE,stringsAsFactors=FALSE)
patient2 <- read.table(opt$listpatient2, header=FALSE,stringsAsFactors=FALSE)


# Create Matrice TPM using MolSubtype
type_cancer <- scan(opt$matricetpm, skip = 0 ,sep=",", nlines = 1, what = character())
id_patient <- scan(opt$matricetpm, skip = 1 ,sep=",", nlines = 1, what = character())


id_patient <- id_patient[-1] #delete column 1
id_patient <- id_patient[-1] #delete column 2
id_patient <- id_patient[-1] #delete column 3
id_patient <- id_patient[-1] #delete column 4
id_patient <- id_patient[-1] #delete column 5
id_patient <- id_patient[-1] #delete column 6
id_patient<-gsub("\\.", "-", id_patient)

id_patient <- as.data.frame(as.list(id_patient))
id_patient_transposed <-  transpose(id_patient)

type_cancer <- type_cancer[-1] #delete column 1
type_cancer <- type_cancer[-1] #delete column 2
type_cancer <- type_cancer[-1] #delete column 3
type_cancer <- type_cancer[-1] #delete column 4
type_cancer <- type_cancer[-1] #delete column 5
type_cancer <- type_cancer[-1] #delete column 6
type_cancer <- as.data.frame(as.list(type_cancer))
type_cancer_transposed <- transpose(type_cancer)


features <- read.table(opt$features,sep="\t", header=FALSE,stringsAsFactors=FALSE )
features <- unique(features)

tpm <- read.table(opt$matricetpm,skip = 2 ,sep=",", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )


head(tpm)
# Set Colnames ID PATIENTS
colnames(tpm)[1] = "Ensembl"
#tpm_filtered <- psi[psi$Event %in% features$V1,] 
colnames(tpm)[2] = "Symbol"
tpm_filtered <- tpm[tpm$Symbol %in% features$V1,] 

max = nrow(id_patient_transposed) + 6
print(max)
for( i in 7:max ){
	#print(i)
	#print(colnames(tpm_filtered)[i])
	colnames(tpm_filtered)[i] = id_patient_transposed[i-6,1]
	#print(colnames(tpm_filtered)[i])
	#print(id_patient_transposed[i-6,1])
	
}

# Filter only features of interest 

# and Patient of interest list 1
vector <- unlist(patient1, use.names=FALSE)

#vector <- c(vector, "Event")
tpm_filteredSubset1 = subset(tpm_filtered, select = vector)


# and Patient of interest list2 
vector <- unlist(patient2, use.names=FALSE)
#vector <- c(vector, "Event")
print(vector)
tpm_filteredSubset2 = subset(tpm_filtered, select = vector)

# Reformat (loosing PatientID but you know they are from patient1 list)
#BE CAREFULL transpose in order of the index !!!!!!!!
vector1 = c()
tpm_filtered_transposed1 <-  transpose(tpm_filteredSubset1)
colnames(tpm_filtered_transposed1) <- rownames(tpm_filteredSubset1)

for( i in rownames(tpm_filteredSubset1)){
	
	test<-subset(tpm_filtered, rownames(tpm_filtered) == i)
	vector1 <- c(vector1, test$Symbol)
}
colnames(tpm_filtered_transposed1) <- vector1

vector2 = c()
tpm_filtered_transposed2 <-  transpose(tpm_filteredSubset2)
colnames(tpm_filtered_transposed2) <- rownames(tpm_filteredSubset2)
#print(tpm_filtered_transposed2)
for( i in rownames(tpm_filteredSubset2) ){
	test<-subset(tpm_filtered, rownames(tpm_filtered) == i)
	vector2<- c(vector2, test$Symbol)

}
colnames(tpm_filtered_transposed2) <- vector2

print("tpm_filtered_transposed2")
#print(tpm_filtered_transposed2)

tpm_filtered_transposed_stack1 <- stack(as.data.frame(tpm_filtered_transposed1))
tpm_filtered_transposed_stack2 <- stack(as.data.frame(tpm_filtered_transposed2))

tpm_filtered_transposed_stack1$group <- "group1"
tpm_filtered_transposed_stack2$group <- "group2"
print("+++")

#  values                           ind
#1  0.996 AAMP_chr2:218267493-218267613

final  <- rbind(tpm_filtered_transposed_stack1, tpm_filtered_transposed_stack2)
dim(final)
#print(final)
toPaste = paste(opt$out,basename(opt$features),sep="_")
write.csv(final,row.names=FALSE,file=paste(toPaste,"_tpm_plot.csv",sep="") ,quote=FALSE)


# Reorganise by the diff
final2= final %>% group_by(ind, group) %>% summarise(avg = mean(values, na.rm=T)) %>% group_by(ind) %>% summarise(diff = diff(avg, na.rm=T)) %>% arrange(desc(diff))

# Reorganise by the mean
#final2= final %>% group_by(ind) %>% summarise(avg = mean(values, na.rm=T))  %>% arrange(avg)
head(final2)
final$ind <- factor(final$ind, levels = as.vector(final2$ind))

write.csv(final2,row.names=FALSE,file=paste(toPaste,"_tpm_diff_.csv",sep="") ,quote=FALSE)


# grouped boxplot
#ggplot(data, aes(x=variety, y=note, fill=treatment)) +  geom_boxplot()
png(paste(opt$out,"_plot.png",sep=""),width=2000,height=1000) 
#ggplot(final, aes(x=ind, y=values, fill=group)) +geom_boxplot()     +  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20))
ggplot(final, aes(x=ind, y=values, fill=group)) +geom_boxplot()  + ylim(0, 250)  +  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20))
 #+ geom_jitter(width=0.1,alpha=0.2) 

dev.off()

  #geom_boxplot() ++  geom_jitter() +  stat_compare_means(hide.ns = TRUE,label = "p.signif",aes(group = group),label.y = 1.1,size =10)

