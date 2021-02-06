library(optparse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

hrbrthemes::import_roboto_condensed()
library(hrbrthemes)


option_list = list(

  make_option(c("-q", "--input"), type="character", default=NULL, help="input", metavar="character")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options


dataframe <- read.table(opt$input, header=TRUE,stringsAsFactors=FALSE,sep=",")
dataframe$X3 <-gsub('Group: ', '', dataframe$X3 )

head(dataframe)
mypalette = c("BASALA" = "#00BFC4", "BASALB"="#F8766D")
#scale_fill_manual("legend", values = c("BASALB-LIKE" = "#F8766D", "BASALA-LIKE" = "#00BFC4")
#dataframe$X0[dataframe$X0=="BASALB"]<- "BASALB-LIKE"
#dataframe$X0[dataframe$X0=="BASALA"]<- "BASALA-LIKE"
head(dataframe$X3)
#:::::::::::::::::::::::::::::::::::::::::::::::::
# linear trend + confidence interval
p3 <- ggplot(dataframe, aes(x=X1, y=X2,color=factor(X0))) +
  geom_point(alpha=0.5) +
  labs(x= "tSNE-1",y="tSNE-2",color='Cell Lines') +  
  geom_label_repel(data = subset(dataframe, X0=="BASALB" & X1 < 0),size = 2.4,aes(label = X3),
                  box.padding   = 2.9, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  theme(legend.position="none",text=element_text(size=10)) +   scale_colour_manual(values = mypalette)

tobepasted = c(dirname(opt$input),"/","Rtsne_scatter",basename(opt$input),".png")
towrite =paste(tobepasted,sep="/",collapse="")

png(towrite,,width=8,height=7,units = 'cm', res = 300)
p3 
dev.off()
