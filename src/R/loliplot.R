library(optparse)
library(ggplot2)
library(dplyr)
hrbrthemes::import_roboto_condensed()
library(hrbrthemes)


option_list = list(
  make_option(c("-q", "--input"), type="character", default=NULL, help="input", metavar="character")


)


parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options
#ASALB  Group: MB157  Other 0.42  0.58
#Other Group: SUM225CWN  Other 0.473 0.527

mypalette = c("BASALA"="#00BFC4",  "BASALB" = "#F8766D")
label = c("BASALA",  "BASALB" )

dataframe <- read.table(opt$input, header=FALSE,stringsAsFactors=FALSE,sep="\t")

dataframe$V2 <-gsub('Group: ', '', dataframe$V2 )

dataframe$V2 <- make.unique(dataframe$V2, sep = ".")

dataframe$V1[dataframe$V1=="Other"] <- "BASALA"
dataframe$V3[dataframe$V3=="Other"] <- "BASALA"

dataframe$V6<- "BASALA"
dataframe$V7<- "BASALB"
# Plot
#V1               V2     V3    V4    V5
#1 BASALB     Group: MB157  Other 0.420 0.580
groupLabeledAs <- ifelse(dataframe$V1 == 'BASALB', '#F8766D', '#00BFC4')
dataframe$V8<-groupLabeledAs

dataframe$V9<-abs(dataframe$V5-dataframe$V4)

reordered_tric <- dataframe[order(dataframe$V9),c(2,8)]
#ggplot(my.df, aes(x=reorder(Attribute, my.order),y=value)) + 

print(reordered_tric$V8)

p <- ggplot(dataframe) +
  geom_segment( aes(x=reorder(V2, V9), xend=V2, y=V4, yend=V5), color="grey" ,size=1)+
  geom_point( aes(x=reorder(V2, V9), y=V4, color=V7), size=2) + 
  geom_point( aes(x=reorder(V2, V9), y=V5,color=V6), size=2) +
  scale_color_manual(name = "",values = mypalette,labels=label )+ 
  #geom_hline(yintercept=0.4, linetype="dashed",color = "yellow", size=1) +
  #geom_hline(yintercept=0.6, linetype="dashed",color = "yellow", size=1) +
  labs( x ="", y = "Probability of being :")  +
  #labels="Forr legend to rename stuffs "
  theme(    legend.position = "bottom",axis.text.x = element_text(size=12,angle = 45),axis.text.y = element_text(size=12,colour=reordered_tric$V8))+
  coord_flip() 

tobepasted = c(dirname(opt$input),"/","LOLIPLOT_",basename(opt$input),".png")
towrite =paste(tobepasted,sep="/",collapse="")
#10
png(towrite,width=10,height=25,units = 'cm', res = 300)
p 
dev.off()