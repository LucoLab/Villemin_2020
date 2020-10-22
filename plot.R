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

  make_option(c("-o", "--outputName"), type="character", default=NULL, help="output", metavar="character"),
  make_option(c("-u", "--input11"), type="character", default=NULL, help="Path to ID best future after boruta", metavar="character")


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
#-a /home/jean-philippe.villemx= "Batch Number") in/data/data/PROJECT/zeroApriori/2020_04_09-1130_False_1000_0//features_final.txt 
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

input11 <- read.table(opt$input11, header=FALSE,stringsAsFactors=FALSE)
Boruta <- data.frame(features=input11$V1)

########################################

# Before I was using the whole set of features Dont use it anymore that was the wole stuff, 
#now use only the set selected at least once
input9 <- read.table(opt$input9, header=FALSE,stringsAsFactors=FALSE)
dataframe6 <- data.frame(features=input9$V1,sum_features=input9$V2,boruta=input9$V3)
print(opt$input9)
input6 <- read.table(opt$input6, header=FALSE,stringsAsFactors=FALSE)
dataframe6_allFeatures <- data.frame(features=input6$V1,sum_features=input6$V2,boruta=input6$V3)

mypalette = c("YES"="red",  "NO" = "gray")
label = c("YES"="Best Features", "NO" = "") 
png(paste(opt$out,"_sum_features.png",sep=""),width=22,height=5,units = 'cm', res = 300)
print("sum_features")
print(dataframe6)
ggplot(dataframe6,aes(x=features,y=sum_features,fill=boruta))   + scale_fill_manual(name="", values = mypalette,labels=label)  +
#ggplot(head(dataframe6,opt$zfeature),aes(x=features,y=sum_features,fill=features))   +
geom_hline(yintercept = median(dataframe6$sum_features), color='yellow',size=0.1) + labs( y = "Score" , x = "Feature") +
 geom_bar(stat="identity",colour="black")  + # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x=element_blank(),text = element_text(size=12))
dev.off()#


dataframe_subset_boruta <- subset(dataframe6, features %in% Boruta$features)
png(paste(opt$out,"_sum_features_boruta.png",sep=""),width=1500,height=1000) 
print("sum_features_boruta")
print(dataframe6)
ggplot(dataframe_subset_boruta,aes(x=features,y=sum_features,fill=features))   +
 geom_bar(stat="identity",colour="black")  + # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
 geom_hline(yintercept = median(dataframe6_allFeatures$sum_features), color='black',size=1) +
theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20))
dev.off()

#print("tail")
#tail(dataframe6,opt$zfeature)
#+ ylim(0, 1.5)
#png(paste(opt$out,"_sum_features_worst.png",sep=""),width=1000,height=1000) 
#ggplot(tail(dataframe6,opt$zfeature),aes(x=features,y=sum_features,fill=features))   + geom_bar(stat="identity",colour="black") + theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
#dev.off()
#+ ylim(0, 1.5)
tobepasted = c(dirname(opt$input4),"/",opt$zfeature,"_features_ordered_sum",".csv")
towrite =paste(tobepasted,sep="/",collapse="")#
write.table(head(dataframe6$features,opt$zfeature),row.names=FALSE,col.names=FALSE,file=towrite ,quote=FALSE)

# the worst
tobepasted = c(dirname(opt$input4),"/",opt$zfeature,"_worst_features_ordered_sum",".csv")
towrite =paste(tobepasted,sep="/",collapse="")
write.table(tail(dataframe6$features,opt$zfeature),row.names=FALSE,col.names=FALSE,file=towrite ,quote=FALSE)

########################################

input1 <- read.table(opt$input1, header=FALSE,stringsAsFactors=FALSE)

# All features at least selected as best feature
dataframe <- data.frame(run=input1$V1,gene=input1$V2,score=input1$V3,boruta=input1$V4 )

# add this part to filter the plot
filter <- as.data.frame(head(dataframe6$features,opt$zfeature))
colnames(filter)<-"id"
filter2 <- as.data.frame(tail(dataframe6$features,opt$zfeature))
colnames(filter2)<-"id_worst"
print("filter2")
head(filter2)

dataframe_subset <- subset(dataframe, gene %in% filter$id)

dataframe_subset2 <- subset(dataframe, gene %in% filter2$id_worst)
print("dataframe_subset2")
print(dataframe_subset2)

png(paste(opt$out,"features.png",sep=""),width=11,height=9,units = 'cm', res = 300)
ggplot(dataframe_subset[dataframe_subset$gene %in% c("PLOD2_chr3:146077861-146077924" , "ENAH_chr1:225504990-225505053"),], aes(x = run , y = score,group=gene,label=gene)) + 
   geom_point(color="red") + geom_line(color="red") + #+   ylim(0, 0.15) + geom_point(aes(color=gene)) 
    labs(y = "Score",x= "Round") +
  #geom_text(aes(label=ifelse(run==max(run),as.character(gene),'')),hjust=0,vjust=0) +
  geom_label_repel(aes(label = ifelse(run==max(run),as.character(gene),'')),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
 theme(legend.position="none",text = element_text(size=12)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
dev.off()

#png(paste(opt$out,"features_worst.png",sep=""),width=2000,height=1000  )
#ggplot(dataframe_subset2, aes(x = run , y = score,label=gene)) + 
#  geom_point(aes(color=gene)) + ylim(0, 0.15) + 
#  geom_line(aes(color=gene)) + 
#geom_label_repel(aes(label = ifelse(run==max(run),as.character(gene),'')),
#                  box.padding   = 0.35, 
#                  point.padding = 0.5,
#                  segment.color = 'grey50') +
# theme(legend.position="bottom",text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
#dev.off()

########################################

input7 <- read.table(opt$input7, header=FALSE,stringsAsFactors=FALSE)

dataframe <- data.frame(run=input7$V1,gene=input7$V2,score=input7$V3,boruta=input7$V4 )
print(opt$input7)
# add this part to filter the plot
filter <- as.data.frame(head(dataframe6$features,opt$zfeature))
colnames(filter)<-"id"
dataframe_subset <- subset(dataframe, gene %in% filter$id)

filter2 <- as.data.frame(tail(dataframe6$features,opt$zfeature))
colnames(filter2)<-"id_worst"
dataframe_subset2 <- subset(dataframe, gene %in% filter2$id_worst)

png(paste(opt$out,"features_normalised.png",sep=""),width=11,height=9,units = 'cm', res = 300)
ggplot(dataframe_subset[dataframe_subset$gene %in% c("PLOD2_chr3:146077861-146077924" , "ENAH_chr1:225504990-225505053"),], aes(x = run , y = score,label=gene)) + 
  geom_point(aes(color=gene)) + 
  geom_line(aes(color=gene)) +  geom_hline(yintercept = 0, color='black',size=1) + labs(y = "Score",x= "Round Number") +
geom_label_repel(aes(label = ifelse(run==max(run),as.character(gene),'')),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
 theme(legend.position="bottom",text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
dev.off()

#png(paste(opt$out,"features_normalised_worst.png",sep=""),width=2000,height=1000  )
#ggplot(dataframe_subset2, aes(x = run , y = score,label=gene)) + 
#  geom_point(aes(color=gene)) +   
#  geom_line(aes(color=gene)) +  geom_hline(yintercept = 0, color='black',size=1) +
#geom_label_repel(aes(label = ifelse(run==max(run),as.character(gene),'')),
#                  box.padding   = 0.35, 
#                  point.padding = 0.5,
#                  segment.color = 'grey50') +
# theme(legend.position="bottom",text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
#dev.off()
#ylim(-3, 3) +
########################################

input2 <- read.table(opt$input2, header=FALSE,stringsAsFactors=FALSE)
head(input2)
dataframe2 <- data.frame(run=input2$V1,type=input2$V2,nombrePatients=input2$V3)

## 00 c'est le vert
png(paste(opt$out,"number_patients_added.png",sep=""), width=11,height=9,units = 'cm', res = 300) 
ggplot(dataframe2, aes(x=factor(run),y=nombrePatients,fill=type))  + scale_fill_manual("legend", values = c("BASALB-LIKE" = "#F8766D", "BASALA-LIKE" = "#00BFC4")) + geom_bar(stat="identity",colour="black",position=position_dodge()) +
labs(y = "Number of Patients Added",x= "Round Number") +
 theme(text = element_text(size=12))
dev.off()
#
#+ ylim(0, 125) palette =c("#00AFBB", "#E7B800")
########################################

input3 <- read.table(opt$input3, header=FALSE,stringsAsFactors=FALSE)
head(input3)
dataframe3 <- data.frame(run=input3$V1,probas=input3$V2)

png(paste(opt$out,"probas.png",sep="") , width=10,height=10,units = 'cm', res = 300) 
#ggplot(dataframe3, aes(x=run, y=probas,color=run) ) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) 
ggplot(dataframe3, aes(x=factor(run), y=probas),group=1 ) + geom_point(position=position_jitter(width=0.2), alpha=0.5)  +  

geom_boxplot(fill=NA,outlier.shape=NA) +
 geom_hline(yintercept = 0.6 , color='yellow',size=1,linetype=2) +
 geom_hline(yintercept = 0.4 , color='yellow',size=1,linetype=2) +
 stat_summary(fun.y=median, geom="line", colour="red",aes(group=1))  +
 labs(x= "Round Number",y="Probability of being classified as BasalB-Like") +  
 theme(legend.position="none",text=element_text(size=12))

dev.off()



#input10 <- read.table(opt$input10, header=FALSE,stringsAsFactors=FALSE)
#head(input10)
#dataframe3 <- data.frame(run=input10$V1,probas=input10$V2)

#png(paste(opt$out,"probas_test.png",sep="") ,width=1000,height=1000  )
#ggplot(dataframe3, aes(x=run, y=probas,color=run) ) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) 
#ggplot(dataframe3, aes(x=factor(run), y=probas,color=run),group=1 ) + geom_point(position=position_jitter(width=0.2), alpha=0.5)  +  
#geom_boxplot(fill=NA) + stat_summary(fun.y=median, geom="line", colour="red",aes(group=1,size=2))  +  theme(text = element_text(size=20))

#dev.off()
########################################

#input4 <- read.table(opt$input4, header=FALSE,stringsAsFactors=FALSE)
#head(input4)
#dataframe4 <- data.frame(run=input4$V1,type=input4$V2,score=input4$V3 )

#png(paste(opt$out,"featuresMean.png",sep=""),width=1000,height=1000  )
#ggplot(dataframe4, aes(x = run , y = score,color=type)) + 
##  geom_point(aes(color=type)) + 
#  geom_line(aes(color=type))  
#dev.off()

########################################

input5 <- read.table(opt$input5, header=FALSE,stringsAsFactors=FALSE)
head(input5)

dataframe5 <- data.frame(feature=input5$V1,run=input5$V2 )
freq.dataframe <- data.frame(table(dataframe5$feature))
png(paste(opt$out,"_occurencies.png",sep=""),width=2000,height=1000  )
ggplot(dataframe5, aes(x=feature, y=run, color=feature)) +  geom_point(size=4) +  theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
dev.off()

png(paste(opt$out,"_num_occurencies.png",sep=""),width=2000,height=1000  )

##freq.dataframe$Freq <- factor(freq.dataframe$Freq, levels = unique(freq.dataframe$Freq))
#head(freq.dataframe[order(freq.dataframe$Freq),])
#ggplot(freq.dataframe, aes(x=Var1, y=Freq, fill=Var1)) + geom_bar(stat="identity",colour="black") + theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
ggplot(freq.dataframe, aes(x=reorder(Var1, -Freq),y=Freq))+ geom_bar(stat="identity",colour="black",fill="skyblue2") + theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +


#ggplot(dataframe5, aes(x = feature)) + geom_histogram(stat = "count")

dev.off()


#tobepasted = c(dirname(opt$input4),"/","RandomState_testing",".txt")
#towrite =paste(tobepasted,sep="/",collapse="")
#inputonsenfout<- read.table(towrite, header=FALSE,stringsAsFactors=FALSE,sep='\t')
#dataframeonsenfout <- data.frame(mesPredicted=inputonsenfout$V2)

#png(paste(opt$out,"RandomState_testing.png",sep="") ,width=1000,height=1000  )
#ggplot(dataframeonsenfout, aes(x=1, y=mesPredicted)) + geom_point(position=position_jitter(width=0.2), alpha=0.5)  +  
#geom_boxplot(fill=NA) +  theme(text = element_text(size=20))

#ggplot(data=dataframeonsenfout, aes(x=mesPredicted)) + geom_bar(aes(y = ..prop.., group = 1))
#dev.off()

#png(paste(opt$out,"_occurencies.png",sep=""),width=1000,height=1000  )


########################################

tobepasted = c(dirname(opt$input4),"/","PVAL_NA",".csv")
towrite =paste(tobepasted,sep="/",collapse="")

input6<- read.table(towrite, header=FALSE,stringsAsFactors=FALSE,sep=',')
head(input6)
dataframe6 <- data.frame(run=input6$V2,endPoint=input6$V3,log10pvalue=-log10(input6$V4) )

png(paste(opt$out,"pvalKaplanMeieir_NA.png",sep=""), width=9,height=8,units = 'cm', res = 300) 
ggplot(dataframe6, aes(x = factor( run), y = log10pvalue,color=endPoint ,group=endPoint)) + 
geom_point(aes(color=endPoint),size=2) +   geom_line(aes(color=endPoint),size=1 )+
 geom_hline(yintercept = -log10(0.05), color='black',size=1) + labs(x= "Round Number")+ theme(text = element_text(size=12))
#color="endpoint"
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
 geom_hline(yintercept = -log10(0.05), color='black',size=1) +  labs(x= "Batch Number") + theme(text = element_text(size=20))

dev.off()
########################################


input8 <- read.table(opt$input8, header=FALSE,stringsAsFactors=FALSE)
head(input8)
dataframe8 <- data.frame(feature=input8$V2,diff_score=input8$V3 )
dataframe_subset <- subset(dataframe8, feature %in% filter$id)
png(paste(opt$out,"_diff_run_stop.png",sep=""),width=1500,height=1000  )
ggplot(dataframe_subset, aes(x=feature,y=diff_score,fill=feature)) + geom_bar(stat="identity",colour="black",position=position_dodge())+  labs(y = "Score difference ") + theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20))
dev.off()


# I COMMENT FOR THE EXPRESSION THIS PART.

########################################
#tobepasted = c(dirname(opt$input4),"/","accuracies_selected",".txt")
#towrite =paste(tobepasted,sep="/",collapse="")

#input_accurary <- read.table(towrite, header=FALSE,stringsAsFactors=FALSE)

#dataframe_accurary <- data.frame(numComparison=input_accurary$V1,accuracy=input_accurary$V2 )
#png(paste(opt$out,"_accuracies.png",sep=""),width=1000,height=1000  )
#ggplot(dataframe_accurary, aes(x=feature, y=run, color=feature)) +  geom_point(size=4) +  theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
#ggplot(dataframe_accurary, aes(x = numComparison , y = accuracy)) + geom_point() + geom_line() + labs(y = "Accuracy", x = "Run") + theme(text = element_text(size=20)) +   ylim(0, 1) +  geom_hline(yintercept = 0.75 , color='red',size=1) 


#dev.off()

########################################
#tobepasted = c(dirname(opt$input4),"/","accuracies_selected1",".txt")
#towrite =paste(tobepasted,sep="/",collapse="")

#input_accurary <- read.table(towrite, header=FALSE,stringsAsFactors=FALSE)

#dataframe_accurary <- data.frame(numComparison=input_accurary$V1,accuracy=input_accurary$V2 )
#png(paste(opt$out,"_accuracies1.png",sep=""),width=1000,height=1000  )
#ggplot(dataframe_accurary, aes(x=feature, y=run, color=feature)) +  geom_point(size=4) +  theme(legend.position="none",axis.text.x = element_text(angle = 90),text = element_text(size=20)) # geom_text(aes(label=ifelse(score>=0.06 && run==3,as.character(gene),'')),hjust=0,vjust=0) +
#ggplot(dataframe_accurary, aes(x = numComparison , y = accuracy)) + geom_point() + geom_line() + labs(y = "Accuracy", x = "Run") + theme(text = element_text(size=20)) +   ylim(0, 1) +  geom_hline(yintercept = 0.75 , color='red',size=1) 

#dev.off()
