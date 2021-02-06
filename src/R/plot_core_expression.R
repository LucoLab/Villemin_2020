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
#Rscript /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/plot.R  
#-x /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//sum_features_final_selectedatLeastOnceInaRun.txt 
#-g /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//features_final_absolute_diff.txt -z 173 
#-n /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//features_final_normalised.txt
 #-e /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//occurrencies.txt 
 #-d /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//probasMean.txt 
 #-c /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//probas.txt 
 #-a /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//features_final.txt 
 #-f /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//sum_features_final.txt 
# -b /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0//patients_added.txt 
 #-o /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_17-2337_0_1000_0/final_plot

library(optparse)
####################################################################################
######################### Parameters  ##############################################
####################################################################################

option_list = list(
  make_option(c("-a", "--input1"), type="character", default=NULL, help="input", metavar="character"),
  make_option(c("-b", "--input2"), type="character", default=NULL, help="Path to file with Patients associated to group", metavar="character"),
  make_option(c("-c", "--input3"), type="character", default=NULL, help="input", metavar="character"),
  make_option(c("-d", "--input4"), type="character", default=NULL, help="input", metavar="character"),
  make_option(c("-e", "--input5"), type="character", default=NULL, help="input", metavar="character"),
  make_option(c("-n", "--input7"), type="character", default=NULL, help="input", metavar="character"),
  make_option(c("-g", "--input8"), type="character", default=NULL, help="input", metavar="character"),

  make_option(c("-x", "--input9"), type="character", default=NULL, help="Subset of features at leat selected one time", metavar="character"),
  make_option(c("-f", "--input6"), type="character", default=NULL, help="All features", metavar="character"),

  make_option(c("-y", "--input10"), type="character", default=NULL, help="input", metavar="character"),

  make_option(c("-z", "--zfeature"), type="integer", default=NULL, help="Number of Best features you want to plot", metavar="character"),

  make_option(c("-o", "--outputName"), type="character", default=NULL, help="output", metavar="character")
  #make_option(c("-u", "--input11"), type="character", default=NULL, help="Path to ID best future after boruta", metavar="character")


)

# e occurrencies input 5
# d proba_means input 4
# c proba input 3
# a features_final  input 1
# b patients_added input 2
# f sum_features_final input 6
# n features_final  input 7
# z number of feature you want to plot
# x sum_features_final input9 maybe to replace with dataframe6
# y proba test 10

#Rscript /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/plot.R 
#-u /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//outputBorutaPy.txt 
#-y /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//probas_test.txt 
#-x /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//sum_features_final_selectedatLeastOnceInaRun.txt 
#-f /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//sum_features_final.txt (217 features)
#-g /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//features_final_absolute_diff.txt 
#-z 217 
#-n /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//features_final_normalised.txt 
#-e /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//occurrencies.txt 
#-d /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//probasMean.txt 
#-c /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//probas.txt 
#-a /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//features_final.txt 
#-b /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//patients_added.txt 
#-o /home/jean-philippe.villemin/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0/final_plot


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
library(ggrepel)

input2 <- read.table(opt$input2, header=FALSE,stringsAsFactors=FALSE)
head(input2)
dataframe2 <- data.frame(run=input2$V1,type=input2$V2,nombrePatients=input2$V3)

## 00 c'est le vert
png(paste(opt$out,"number_patients_added.png",sep=""),width=1000,height=1000  )
ggplot(dataframe2, aes(x=factor(run),y=nombrePatients,fill=type))  + scale_fill_manual("legend", values = c("BASALB-LIKE" = "#F8766D", "BASALA-LIKE" = "#00BFC4")) + geom_bar(stat="identity",colour="black",position=position_dodge())+  labs(y = "Number of Patients Added",x= "Batch Number") + theme(text = element_text(size=20))
dev.off()
#
#+ ylim(0, 125) palette =c("#00AFBB", "#E7B800")
########################################

input3 <- read.table(opt$input3, header=FALSE,stringsAsFactors=FALSE)
head(input3)
dataframe3 <- data.frame(run=input3$V1,probas=input3$V2)

dataframe3 <- dataframe3 %>% 
  mutate(mycolor = ifelse(probas > 0.6, "#F8766D", ifelse(probas < 0.4 , "#00BFC4", "black")))

png(paste(opt$out,"probas.png",sep="") , width=10,height=10,units = 'cm', res = 300) 
#ggplot(dataframe3, aes(x=run, y=probas,color=run) ) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) 
ggplot(dataframe3, aes(x=factor(run), y=probas),group=1 ) + geom_point(colour=dataframe3$mycolor,position=position_jitter(width=0.2), alpha=0.5)  +  

#geom_boxplot(fill=NA,outlier.shape=NA) +
 geom_hline(yintercept = 0.4 , color='#F0E442',size=1,linetype=1) +
 geom_hline(yintercept = 0.6 , color='#F0E442',size=1,linetype=1) +
 #stat_summary(fun.y=median, geom="line", colour="red",aes(group=1))  +
 labs(x= "",y="") +  
 theme(legend.position="none",text=element_text(size=12))
 #labs(x= "Round Number",y="Probability of being classified as BasalB-Like") 
dev.off()



tobepasted = c(dirname(opt$input4),"/","PVAL_NA",".csv")
towrite =paste(tobepasted,sep="/",collapse="")

input6<- read.table(towrite, header=FALSE,stringsAsFactors=FALSE,sep=',')
head(input6)
dataframe6 <- data.frame(run=input6$V2,endPoint=input6$V3,log10pvalue=-log10(input6$V4) )

png(paste(opt$out,"pvalKaplanMeieir_NA.png",sep=""),width=1000,height=1000  )
ggplot(dataframe6, aes(x = run , y = log10pvalue,color=endPoint)) + 
  geom_point(aes(color=endPoint)) + 
  geom_line(aes(color=endPoint))  +
 geom_hline(yintercept = -log10(0.05), color='black',size=1) +  theme(text = element_text(size=20))

dev.off()

tobepasted = c(dirname(opt$input4),"/","PVAL_ALL",".csv")
towrite =paste(tobepasted,sep="/",collapse="")

input6<- read.table(towrite, header=FALSE,stringsAsFactors=FALSE,sep=',')
head(input6)
dataframe6 <- data.frame(run=input6$V2,endPoint=input6$V3,log10pvalue=-log10(input6$V4) )

png(paste(opt$out,"pvalKaplanMeieir_ALL.png",sep=""),width=1000,height=1000  )
ggplot(dataframe6, aes(x = run , y = log10pvalue,color=endPoint)) + 
  geom_point(aes(color=endPoint)) + 
  geom_line(aes(color=endPoint))  +
 geom_hline(yintercept = -log10(0.05), color='black',size=1) +  theme(text = element_text(size=20))

dev.off()



