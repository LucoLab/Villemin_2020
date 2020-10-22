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
library(ggsignif)

####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-m", "--matricetpm"), type="character", default=NULL, help="Absolute File Input Path for PSI", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, help="Features ", metavar="character") ,
  make_option(c("-o", "--outputName"), type="character", default=NULL, help="IDs", metavar="character")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options


features <- read.table(opt$features,sep="\t", header=FALSE,stringsAsFactors=FALSE )
features <- unique(features)

# Matrice TPM has two rows of annotation :
# First row : Subtype
# Second row : TCGA ID
# Matrice TPM has 6 columns annotation :
# ensembl_gene_id,hgnc_symbol,gene_biotype,Chr,Start,End

# Read Once to get the list of ID for patients 
id_patient <- scan(opt$matricetpm, skip = 1 ,sep=",", nlines = 1, what = character())
id_patient <- id_patient[-1:-6] #delete column 1
id_patient<-gsub("\\.", "-", id_patient)
id_patient <- as.data.frame(as.list(id_patient))
id_patient_transposed <-  transpose(id_patient)

type <- scan(opt$matricetpm, skip = 0 ,sep=",", nlines = 1, what = character())
type <- type[-1:-6]
type <- as.data.frame(as.list(type))
type_transposed <-  transpose(type)


# Read second times to get the list of values
tpm <- read.table(opt$matricetpm,skip = 2 ,sep=",", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )

# Set Colnames ID PATIENTS
colnames(tpm)[1] = "Ensembl"
colnames(tpm)[2] = "Symbol"
tpm_filtered <- tpm[tpm$Symbol %in% features$V1,] 

index_basalA <- which(type_transposed ==  "BASALA")
patient1 <- id_patient_transposed[index_basalA,]

index_basalB <- which(type_transposed ==  "BASALB")
patient2 <- id_patient_transposed[index_basalB,]


max = nrow(id_patient_transposed) + 6
print(max)
for( i in 7:max ){
	colnames(tpm_filtered)[i] = id_patient_transposed[i-6,1]
}


#  Patient of interest list 1
vector <- unlist(patient1, use.names=FALSE)
tpm_filteredSubset1 = subset(tpm_filtered, select = vector)


#  Patient of interest list2 
vector <- unlist(patient2, use.names=FALSE)
tpm_filteredSubset2 = subset(tpm_filtered, select = vector)

# Reformat (loosing PatientID but you know they are from patient1 list)
# BE CAREFULL transpose in order of the index !!!!!!!!
vector1 = c()
tpm_filtered_transposed1 <-  transpose(tpm_filteredSubset1)
colnames(tpm_filtered_transposed1) <- rownames(tpm_filteredSubset1)

for( i in rownames(tpm_filteredSubset1)){
	test    <-subset(tpm_filtered, rownames(tpm_filtered) == i)
	vector1 <- c(vector1, test$Symbol)
}
colnames(tpm_filtered_transposed1) <- vector1

vector2 = c()
tpm_filtered_transposed2 <-  transpose(tpm_filteredSubset2)
colnames(tpm_filtered_transposed2) <- rownames(tpm_filteredSubset2)
#print(tpm_filtered_transposed2)
for( i in rownames(tpm_filteredSubset2) ){
	test    <-subset(tpm_filtered, rownames(tpm_filtered) == i)
	vector2 <- c(vector2, test$Symbol)

}
colnames(tpm_filtered_transposed2) <- vector2

tpm_filtered_transposed_stack1 <- stack(as.data.frame(tpm_filtered_transposed1))
tpm_filtered_transposed_stack2 <- stack(as.data.frame(tpm_filtered_transposed2))

tpm_filtered_transposed_stack1$group <- "BASALA"
tpm_filtered_transposed_stack2$group <- "BASALB"

final  <- rbind(tpm_filtered_transposed_stack1, tpm_filtered_transposed_stack2)


# Reorganise by the diff
final2= final %>% group_by(ind, group) %>% summarise(avg = mean(values, na.rm=T)) %>% group_by(ind) %>% summarise(diff = diff(avg, na.rm=T)) %>% arrange(desc(diff))
final$ind <- factor(final$ind, levels = as.vector(final2$ind))

# Reorganise by the mean
#final2= final %>% group_by(ind) %>% summarise(avg = mean(values, na.rm=T))  %>% arrange(avg)
#final$ind <- factor(final$ind, levels = as.vector(final2$ind))

toPaste = paste(opt$out,basename(opt$features),sep="_")
write.csv(final2,row.names=FALSE,file=paste(toPaste,"_tpm_diff_cellLines.csv",sep="") ,quote=FALSE)
write.csv(final,row.names=FALSE,file=paste(toPaste,"_tpm_plot_cellLines.csv",sep="") ,quote=FALSE)

# Control Wilcoxon is ok for one example.
#print("ZEB1")
#zeb1 <- final[final$ind=="ZEB1", ]
#test = wilcox.test(zeb1$values ~ zeb1$group ,paired = FALSE)
#print(test)


# TRANSFORMATION EN LOG2(TPM+1)
final$values <- log2(final$values + 1)

max = max(final$values,na.rm=TRUE)
mypalette = c("BASALA" = "#00BFC4", "BASALB"="#F8766D")
print(max)
png(paste(opt$out,"_EXPRESSION_plot.png",sep=""), width=9,height=8,units = 'cm', res = 300)
ggplot(final, aes(x=ind, y=values, fill=group)) + geom_boxplot(outlier.shape = NA)  + ylim(0, max+1.5)  +  
theme(legend.title = element_blank() ,axis.text.x = element_text(angle = 90,size=12)) +    scale_fill_manual(values = mypalette)+    labs(x= "",y="LOG2(TPM+1)")  +
stat_compare_means(hide.ns = TRUE,label = "p.signif",aes(group = group),label.y = max+1 ,size =4)#,size =12
dev.off()
