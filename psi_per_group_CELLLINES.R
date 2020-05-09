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
  make_option(c("-m", "--matricePsi"), type="character", default=NULL, help="Absolute File Input Path for PSI", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, help="Features ", metavar="character") ,
  make_option(c("-o", "--outputName"), type="character", default=NULL, help="IDs", metavar="character")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options


id_patient <- scan(opt$matricePsi, skip = 0 ,sep="\t", nlines = 1, what = character())
id_patient <- id_patient[-1:-3] #delete column 1
id_patient <-gsub("Patient: ","",id_patient)
id_patient <- as.data.frame(as.list(id_patient))
id_patient_transposed <-  transpose(id_patient)


type <- scan(opt$matricePsi, skip = 3 ,sep="\t", nlines = 1, what = character())
type <- type[-1:-3]
type <-gsub("Group2: ","",type)
type <- as.data.frame(as.list(type))
type_transposed <-  transpose(type)

features <- read.table(opt$features,sep="\t", header=FALSE,stringsAsFactors=FALSE )

psi <- read.table(opt$matricePsi,skip = 4 ,sep="\t", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )


# Set Colnames ID PATIENTS
colnames(psi)[1] = "Event"
psi_filtered <- psi[psi$Event %in% features$V1,] 


index_basalA <- which(type_transposed ==  "BASALA")
patient1 <- id_patient_transposed[index_basalA,]

index_basalB <- which(type_transposed ==  "BASALB")
patient2 <- id_patient_transposed[index_basalB,]


max = nrow(id_patient_transposed) + 3
for( i in 4:max ){
	colnames(psi_filtered)[i] = id_patient_transposed[i-3,1]
}

#  Patient of interest list 1
vector <- unlist(patient1, use.names=FALSE)
psi_filteredSubset1 = subset(psi_filtered, select = vector)

#  Patient of interest list2 
vector <- unlist(patient2, use.names=FALSE)
psi_filteredSubset2 = subset(psi_filtered, select = vector)

# Reformat (loosing PatientID but you know they are from patient1 list)
#BE CAREFULL transpose in order of the index !!!!!!!!
vector1 = c()
psi_filtered_transposed1 <-  transpose(psi_filteredSubset1)
colnames(psi_filtered_transposed1) <- rownames(psi_filteredSubset1)

for( i in rownames(psi_filteredSubset1)){
	test    <- subset(psi_filtered, rownames(psi_filtered) == i)
	vector1 <- c(vector1, test$Event)
}
colnames(psi_filtered_transposed1) <- vector1

vector2 = c()
psi_filtered_transposed2 <-  transpose(psi_filteredSubset2)
colnames(psi_filtered_transposed2) <- rownames(psi_filteredSubset2)

for( i in rownames(psi_filteredSubset2) ){
	test    <- subset(psi_filtered, rownames(psi_filtered) == i)
	vector2 <- c(vector2, test$Event)

}
colnames(psi_filtered_transposed2) <- vector2


psi_filtered_transposed_stack1 <- stack(as.data.frame(psi_filtered_transposed1))
psi_filtered_transposed_stack2 <- stack(as.data.frame(psi_filtered_transposed2))

psi_filtered_transposed_stack1$group <- "BASALA"
psi_filtered_transposed_stack2$group <- "BASALB"


final  <- rbind(psi_filtered_transposed_stack1, psi_filtered_transposed_stack2)


write.csv(final,row.names=FALSE,file=paste(opt$out,"psi_plot.csv",sep="") ,quote=FALSE)

# To replace ID you need to change as character
final$ind <- as.character(final$ind)

final$ind[final$ind=="CTNND1_chr11:57789036-57789155"] <- "CTNND1_a"
final$ind[final$ind=="CTNND1_chr11:57791493-57791673"] <- "CTNND1_b"
final$ind[final$ind=="DNM2_chr19:10796060-10796199"] <- "DNM2_a"
final$ind[final$ind=="DNM2_chr19:10808568-10808580"] <- "DNM2_b"
final$ind  <-gsub("_chr.*","",final$ind)
# Re-factorize with the as.factor function or simple factor(fixed$Type)
final$ind <- as.factor(final$ind)

# Reorganise by the diff
#final2= final %>% group_by(ind, group) %>% summarise(avg = mean(values, na.rm=T)) %>% group_by(ind) %>% summarise(diff = diff(avg, na.rm=T)) %>% arrange(desc(diff))
# Reorganise by the mean
final2= final %>% group_by(ind) %>% summarise(avg = mean(values, na.rm=T))  %>% arrange(avg)
final$ind <- factor(final$ind, levels = as.vector(final2$ind))

#final$ind_renamed  <-gsub("_.*","",final$ind)

mypalette = c("BASALA" = "#00BFC4", "BASALB"="#F8766D")

png(paste(opt$out,"_PSI_plot.png",sep=""), width=16,height=8,units = 'cm', res = 300)
ggplot(final, aes(x=ind, y=values, fill=group)) + geom_boxplot(outlier.shape = NA)  + ylim(0, 1.1)  +  
theme(legend.title = element_blank() ,legend.position="top",legend.text= element_text(size=12),axis.text.y = element_text(size=12),axis.text.x = element_text(face="bold", color="black", angle=90,size=12)) +    scale_fill_manual(values = mypalette)+    labs(x= "",y="PSI")  +
stat_compare_means(hide.ns = TRUE,label = "p.signif",aes(group = group),label.y = 1 + 0.05 ,size=3)
dev.off()#,vjust=-0.1,hjust=-0.1
#