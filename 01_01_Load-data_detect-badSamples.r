###### ------  Load packages --------------------------------------------------------------------------
######

library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(expss)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(car)
library(emmeans)

###### ------  Load data --------------------------------------------------------------------------
######

data<-read.csv("methylation.bed", header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings = ("None"))  #Here we also turn all Nones into NA with the commande na.string = ("None")

###### ------  control step --------------------------------------------------------------------------
###### We're going to create R objects to check for mistakes.
###### sink() function saves the results in a .txt in the working directory.
###### dim(data [1]) Check the row numbers. In our case, each row is one C
###### str() function shows the column class.
###### The data has a class that isn't what we need.
###### Take a look.


control.01<-dim(data)[1]  ## save
control.03<-head(data[10:12])

str(data) ## Quick check

###### We're going to change data class to work on this data later.

row.names(data)<-c(1:dim(data)[1])
data[, 4:dim(data)[2]]<-lapply(data[,4:dim(data)[2]], as.numeric, na.rm=T)
data[, 1:3]<-lapply(data[,1:3], as.factor)
str(data) ## Quick check

---NAs introduced by coercion'

###### ------  Quantile by sample_called and read coverage ------------------------------------------------------
###### Now, we're going to calculate quantiles for samples_called.
###### This mean that we want to know how many samples have each C with at least one count or read per sample.
###### If we have 40 samples, and the C (one row) is 40, that means that all samples
###### presented that C.

data.quantile.samples_called<-quantile(data$samples_called, c(0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1), na.rm=TRUE)
data.quantile.samples_called

###### ------  control step --------------------------------------------------------------------------------------
###### We're going to create R objects to check for mistakes.
###### You can compare all control objects to find mistakes in the data set


control.04<-dim(data)[1]
control.06<-head(data[10:12])

save(data, file="01.1_data.Rdata")  #  You can load this data with the following command load("01.1_data.Rdata")


###### ------   Descriptive step: total reads for each sample ------------------------------------------------------
######  We'll work only with the total reads for each samples.

data.total<-select(data,-(ends_with("methylated")))
head(data.total)
nas<-sum(is.na(data.total))  # calculate number of NA
nas.perc<-(nas*100)/(dim(data.total)[1]*dim(data.total)[2]-4) # percentage of NA
nas.perc  ### % of NA in all data


###### ------  From wide to stacked data set  ----------------------------------------------------------------------
###### We need the original table with other format to work faster.
###### Also, we're going to calculate quantiles for read coverage


data.stacked<-melt(data.total, id=c("chr", "pos", "context", "samples_called")) # change to stacked Data

head(data.stacked)

data.quantile.reads<-quantile(data.stacked$value,  c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,
0.95,0.96,0.97,0.98,0.981,0.985,0.989,0.99,0.991,0.995,0.998,0.999,0.9995,0.9999,0.99995, 1), na.rm=TRUE)

data.quantile.reads

###### ------   Samples summary--------------------------------------------------------------------------------------
###### In order to detect bad samples, we're going to calculate: the *sum* of all cytosine reads for each sample.
###### The *average* of the sum and the *50% of the average*.
###### The bad samples are all those samples below the 50% of sum average.
###### You can see all this info in the plot

summ<-ddply(data.stacked, c("variable"), summarise, Sum=sum(value, na.rm=TRUE)) ## sum all reads per sample
head(summ)
average<-mean(summ$Sum)
mean.av<-average*0.5
bad.samples<-filter(summ, Sum<=c(mean.av))
bad.samples

###### ------ Plot ---------------------------------------------------------------------------------------------------
###### And we make a pretty plot and then, we save it.
###### The black line is the average of the sum.
###### The red line is the threshold to detect bad samples

### NOT GENERAL STEP--- CHECK your experimental design and modify

summ<- summ %>% mutate(Treat= str_replace_all(variable, "_total", ""))
summ<-separate(summ, Treat, c("Acc", "Treat"))
summ<- arrange(summ, Acc, Treat)
summ<- summ %>% mutate(Exp= str_replace_all(variable, "_total", ""))
summ$Exp <- factor(summ$Exp, levels=unique(summ$Exp))

ggplot(summ, aes(x=Exp, y=Sum)) +
    geom_bar(stat="identity", width=0.75, color="pink", fill="pink", size=0.1)+
    geom_hline(yintercept=c(average), color="black")+
    geom_hline(yintercept=c(mean.av), linetype="dashed", color="red")+
    geom_text(aes(6,c(mean.av), label=c(mean.av), vjust=0.1))+
    theme(axis.text.x=element_text(angle=90,hjust=1))

tiff("Total reads for sample_02.tiff")
  ggplot(summ, aes(x=variable, y=Sum)) +
    geom_bar(stat="identity", width=0.75, color="pink", fill="pink", size=0.1)+
    geom_hline(yintercept=c(average), color="black")+
    geom_hline(yintercept=c(mean.av), linetype="dashed", color="red")+
    geom_text(aes(6,c(mean.av), label=c(mean.av), vjust=0.1))+
    theme(axis.text.x=element_text(angle=90,hjust=1))
dev.off()

###### ------ Report --------------------------------------------------------------------------------------------------
###### In the text file you will find the data structure from bed file and the data structure after changing the class.
###### Also, all info for read coverage and distribution

sink("01.1_Report: Control uploaded data and detection of bad samples")
cat("#########   CONTROL  #######\n")
cat("#########   Before change data class  #######\n")
        cat("----- Number of Cs from bed file\n")
                control.01
        cat("-----One_sample_data_class\n")
                control.03
cat("#########   After change data class  #######\n")
cat("#########   Check that observations remain the same  #######\n")
        cat("----- Number of Cs from bed file\n")
                control.04
        cat("----- Data_class\n")
                str(data)
        cat("----- One sample data class\n")
                control.06
        cat("----- Quantiles for samples called\n")
                data.quantile.samples_called
        cat("---- NA % original bed file \n")
                nas.perc
        cat("---- Quantiles for read coverage in bed file\n")
                data.quantile.reads
cat("#########   BAD SAMPLES  #######\n")
cat("#########   Before change data class  #######\n")
        cat("Samples below sum average: BAD SAMPLES\n")
                bad.samples
sink()
