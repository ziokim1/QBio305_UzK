### Load packages----
library(tidyverse)

#set working directory
setwd("~/Downloads/Group2_data/admixture")

cv <- read.table("cross_validation.txt")

# Analyze the cross-validation results Then, add a K-cluster column indicating the number of K you test and select only two columns of interest, CV and K.

cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- select(cv, V4,K)
CV$K <- as.numeric(sub("K=", "", CV$K))

# Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

# Do a graph showing the cross validation results. Then select the optimal number of clusters regarding :
# the lowest cross validation error
# when the cross-validation error decrease the most

graph_title=""
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))

#Save the graph
#ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
#dev.off()

# using the palette function and the colors parameter. You can define a custom palette
# with these colors. Here's how you can do it:
# Define a custom color palette
my_colors <- c("red", "yellow", "green", "cyan", "blue")

# Set the custom palette
palette(my_colors)

#load admxiture function
source("admixFun.R") 

#list Q files and sort for K
files <- list.files("~/Downloads/Group2_data/admixture", full = TRUE, pattern = "Q")
files <- files[order(as.numeric(sub('.*a.thaliana_admix.(\\d+)\\.Q$', '\\1', files)))]

#population file
pop <- scan("~/Downloads/Group2_data/admixture/pop.txt",what="df",na="")
table(pop) 

# possible K
Kall <- 1:10

## read Qs
allQ <- list()
for(K in Kall)
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))

#make smaller line and change type to solid line
plotMulti(allQ,Kall=3:7,reorder=1,pop,fast=T,lwd=1,lty=1)

#writeLines(capture.output(sessionInfo()), "sessionInfo_admixture.txt")
