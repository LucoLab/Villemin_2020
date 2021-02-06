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
  make_option(c("-m", "--matricePsi"), type="character", default=NULL, help="Absolute File Input Path for PSI", metavar="character"),
  make_option(c("-a", "--listpatient1"), type="character", default=NULL, help="IDs", metavar="character"),
  make_option(c("-b", "--listpatient2"), type="character", default=NULL, help="IDs", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, help="Features ", metavar="character") ,
  make_option(c("-o", "--outputName"), type="character", default=NULL, help="IDs", metavar="character")

)

#Rscript /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/psi_per_group.R -f /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/features.tsv -m /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/Matrice/2dataset.166.tsv -a /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/test_Patients_MES.txt -b /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/test_Patients_EPI.txt 
#Rscript /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/psi_per_group.R -f /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/features.tsv -m /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/Matrice/output166.intersect.indexed.sorted.bed.tsv -a /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/test_Patients_MES.txt -b /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/test_Patients_EPI.txt 


parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

patient1 <- read.table(opt$listpatient1, header=FALSE,stringsAsFactors=FALSE)
patient2 <- read.table(opt$listpatient2, header=FALSE,stringsAsFactors=FALSE)

# Create Matrice PSI using MolSubtype
type_cancer <- scan(opt$matricePsi, skip = 3 ,sep="\t", nlines = 1, what = character())
id_patient <- scan(opt$matricePsi, skip = 0 ,sep="\t", nlines = 1, what = character())

#head(type_cancer)
id_patient <- id_patient[-1] #delete column 1
#id_patient <- id_patient[-1] #delete column 2 # BE CAREFULL IT CAN HAVE ONE OR TWO COL
id_patient<-gsub("Patient: ","",id_patient)
id_patient<-gsub("_.*","",id_patient)
id_patient <- as.data.frame(as.list(id_patient))

id_patient_transposed <-  transpose(id_patient)


features <- read.table(opt$features,sep="\t", header=FALSE,stringsAsFactors=FALSE )

psi <- read.table(opt$matricePsi,skip = 32 ,sep="\t", header=FALSE,stringsAsFactors=FALSE, comment.char = "@",  quote="" )

head(features)
# Set Colnames ID PATIENTS
colnames(psi)[1] = "Event"
psi_filtered <- psi[psi$Event %in% features$V1,] 


max = nrow(id_patient_transposed) + 1
for( i in 2:max ){
	#print(colnames(psi_filtered)[i])
	colnames(psi_filtered)[i] = id_patient_transposed[i-1,1]
	#print(colnames(psi_filtered)[i])
	#print(id_patient_transposed[i-1,1])
	
}
print("psi_filtered")
print(psi_filtered[psi_filtered$Event=="PLOD2_chr3:146077861-146077924",])

# Filter only features of interest 

# and Patient of interest list 1
vector <- unlist(patient1, use.names=FALSE)

#vector <- c(vector, "Event")
psi_filteredSubset1 = subset(psi_filtered, select = vector)

print(subset(psi_filteredSubset1, rownames(psi_filteredSubset1) == 123))

# and Patient of interest list2 
vector <- unlist(patient2, use.names=FALSE)
#vector <- c(vector, "Event")
print(vector)
psi_filteredSubset2 = subset(psi_filtered, select = vector)

print("psi_filteredSubset2")
print(psi_filteredSubset2["123",])

# Reformat (loosing PatientID but you know they are from patient1 list)
#BE CAREFULL transpose in order of the index !!!!!!!!
vector1 = c()
psi_filtered_transposed1 <-  transpose(psi_filteredSubset1)
colnames(psi_filtered_transposed1) <- rownames(psi_filteredSubset1)

for( i in rownames(psi_filteredSubset1)){
	
	test<-subset(psi_filtered, rownames(psi_filtered) == i)
	vector1 <- c(vector1, test$Event)
}
colnames(psi_filtered_transposed1) <- vector1

vector2 = c()
psi_filtered_transposed2 <-  transpose(psi_filteredSubset2)
colnames(psi_filtered_transposed2) <- rownames(psi_filteredSubset2)
print(psi_filtered_transposed2)
for( i in rownames(psi_filteredSubset2) ){
	test<-subset(psi_filtered, rownames(psi_filtered) == i)
	vector2<- c(vector2, test$Event)

}
colnames(psi_filtered_transposed2) <- vector2

print("psi_filtered_transposed2")
print(psi_filtered_transposed2)

psi_filtered_transposed_stack1 <- stack(as.data.frame(psi_filtered_transposed1))
psi_filtered_transposed_stack2 <- stack(as.data.frame(psi_filtered_transposed2))

psi_filtered_transposed_stack1$group <- "group1"
psi_filtered_transposed_stack2$group <- "group2"
#head(psi_filtered_transposed_stack1[psi_filtered_transposed_stack1$ind=="PLOD2_chr3:146077861-146077924",])
print("+++")
#head(psi_filtered_transposed_stack2[psi_filtered_transposed_stack2$ind=="PLOD2_chr3:146077861-146077924",])

#  values                           ind
#1  0.996 AAMP_chr2:218267493-218267613

final  <- rbind(psi_filtered_transposed_stack1, psi_filtered_transposed_stack2)
dim(final)
#print(final)

write.csv(final,row.names=FALSE,file=paste(opt$out,"psi_plot.csv",sep="") ,quote=FALSE)


# Reorganise by the diff
final2= final %>% group_by(ind, group) %>% summarise(avg = mean(values, na.rm=T)) %>% group_by(ind) %>% summarise(diff = diff(avg, na.rm=T)) %>% arrange(desc(diff))

# Reorganise by the mean
#final2= final %>% group_by(ind) %>% summarise(avg = mean(values, na.rm=T))  %>% arrange(avg)
head(final2)
final$ind <- factor(final$ind, levels = as.vector(final2$ind))



# grouped boxplot
#ggplot(data, aes(x=variety, y=note, fill=treatment)) +  geom_boxplot()
png(paste(opt$out,"_plot.png",sep=""),width=1000,height=1000) 
#ggplot(final, aes(x=ind, y=values, fill=group)) +geom_boxplot()     +  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20))
ggplot(final, aes(x=ind, y=values, fill=group)) +geom_boxplot()    +  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20))
 #+ geom_jitter(width=0.1,alpha=0.2) 

dev.off()

  #geom_boxplot() ++  geom_jitter() +  stat_compare_means(hide.ns = TRUE,label = "p.signif",aes(group = group),label.y = 1.1,size =10)

