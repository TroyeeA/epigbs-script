###### ------PCA -------------------------------------------------------------------------------
###### For each context, we're going to make two PCAs. One with missing data (NA) and the other without NA.
###### load("02.1_data.Rdata") if you start from here
###### For each context, we're going to make two PCAs. One with missing data (NA) and the other without NA.


###### ------With NA-------------------------------------------------------------------------------
###### First, we're going to  for total and methylated samples.
###### Then, we calculate methylation proportion: meth/total
###### We need to impute NA and then take only the complete observacion with imputed data.


###### separate by context
        data.3.80.total.CHH<-select(filter(data.3.80, context == "CHH"), -(ends_with("methylated")), -context)
        data.3.80.meth.CHH<-select(filter(data.3.80, context == "CHH"), -(ends_with("total")), -context)

        data.3.80.total.CHG<-select(filter(data.3.80, context == "CHG"), -(ends_with("methylated")), -context)
        data.3.80.meth.CHG<-select(filter(data.3.80, context == "CHG"), -(ends_with("total")), -context)

        data.3.80.total.CG<-select(filter(data.3.80, context == "CG"), -(ends_with("methylated")), -context)
        data.3.80.meth.CG<-select(filter(data.3.80, context == "CG"), -(ends_with("total")), -context)

###### calculate methylation proportion by context

         new.3.80.CHH<-data.3.80.meth.CHH[,3:dim(data.3.80.meth.CHH)[2]]/data.3.80.total.CHH[,3:dim(data.3.80.total.CHH)[2]]
         new.3.80.CHG<-data.3.80.meth.CHG[,3:dim(data.3.80.meth.CHG)[2]]/data.3.80.total.CHG[,3:dim(data.3.80.total.CHG)[2]]
         new.3.80.CG<-data.3.80.meth.CG[,3:dim(data.3.80.meth.CG)[2]]/data.3.80.total.CG[,3:dim(data.3.80.total.CG)[2]]

###### impute NA
         imp.3.80.CHH<- imputePCA(new.3.80.CHH, ncp=10)
         imp.3.80.CHG<- imputePCA(new.3.80.CHG, ncp=10)
         imp.3.80.CG<- imputePCA(new.3.80.CG, ncp=10)

###### extract imputed observations
        CHH.80<-imp.3.80.CHH$completeObs
        CHG.80<-imp.3.80.CHG$completeObs
        CG.80<-imp.3.80.CG$completeObs
###### transpose the data. We need samples in rows and C's in columns
        CHH.80.transpose<-as.data.frame(t(as.matrix(CHH.80)))
        CHG.80.transpose<-as.data.frame(t(as.matrix(CHG.80)))
        CG.80.transpose<-as.data.frame(t(as.matrix(CG.80)))

###### NOT GENERAL PART: extract names to plot with colors
        y.CHH.80<-rownames(CHH.80.transpose)
        label.CHH.80<-do.call(rbind, strsplit(y.CHH.80, '_'))
        CHH.80$Treat<-label.CHH.80[,2]
        CHH.80$Sample<-label.CHH.80[,1]

        y.CHG.80<-rownames(CHG.80.transpose)
        label.CHG.80<-do.call(rbind, strsplit(y.CHG.80, '_'))
        CHG.80$Treat<-label.CHG.80[,2]
        CHG.80$Sample<-label.CHG.80[,1]

        y.CG.80<-rownames(CG.80.transpose)
        label.CG.80<-do.call(rbind, strsplit(y.CG.80, '_'))
        CG.80$Treat<-label.CG.80[,2]
        CG.80$Sample<-label.CG.80[,1]
###### do PCA
        pca.CHH.80 <- PCA(CHH.80.transpose, graph=FALSE)
        pca.CHG.80 <- PCA(CHG.80.transpose, graph=FALSE)
        pca.CG.80 <- PCA(CG.80.transpose, graph=FALSE)

#####  Save Whitout labels
        plot.CHH.80<-fviz_pca_ind(pca.CHH.80, geom.ind =c("point"),
        col.ind=CHH.80$Treat, repel=TRUE, title="CHH_80")

        tiff("CHH_80.tiff")
        print(plot.CHH.80)
        dev.off()

        plot.CHG.80<-fviz_pca_ind(pca.CHG.80, geom.ind =c("point"),
        col.ind=CHG.80$Treat, repel=TRUE, title="CHG_80")

        tiff("CHG_80.tiff")
        print(plot.CHG.80)
        dev.off()

        plot.CG.80<-fviz_pca_ind(pca.CG.80, geom.ind =c("point"),
        col.ind=CG.80$Treat, repel=TRUE, title="CG_80")

        tiff("CG_80.tiff")
        print(plot.CG.80)
        dev.off()

#####  Save Whit labels
        plot.CHH.80<-fviz_pca_ind(pca.CHH.80, geom.ind =c("point", "text"),
        col.ind=CHH.80$Treat, repel=TRUE, title="CHH_80")

        tiff("CHH_80_Labels.tiff")
        print(plot.CHH.80)
        dev.off()

        plot.CHG.80<-fviz_pca_ind(pca.CHG.80, geom.ind =c("point", "text"),
        col.ind=CHG.80$Treat, repel=TRUE, title="CHG_80")

        tiff("CHG_80_Labels.tiff")
        print(plot.CHG.80)
        dev.off()

        plot.CG.80<-fviz_pca_ind(pca.CG.80, geom.ind =c("point", "text"),
        col.ind=CG.80$Treat, repel=TRUE, title="CG_80")

        tiff("CG_80_Labels.tiff")
        print(plot.CG.80)
        dev.off()

###### ------WithOUT NA-------------------------------------------------------------------------------
###### First, we need to select Cs without NA.
###### We need to separete by context for total and methylated samples.
###### Then, we calculate methylation proportion: meth/total
###### Plot.

data.4.total<-select(data.3.80, -(ends_with("methylated")))
data.4.total$N_NA<-count_row_if(NA, data.4.total[,4:dim(data.4.total)[2]]) # we add a column
data.4.80<-data.3.80
data.4.80$N_NA<-data.4.total$N_NA
data.4.80<-select(filter(data.4.80, N_NA ==0), -N_NA)
head(data.4.80)
dim(data.4.80)[1]

###### separate by context

        data.4.80.total.CHH<-select(filter(data.4.80, context == "CHH"), -(ends_with("methylated")), -chr, -pos, -context)
        data.4.80.meth.CHH<-select(filter(data.4.80, context == "CHH"), -(ends_with("total")),  -chr, -pos, -context)

        data.4.80.total.CHG<-select(filter(data.4.80, context == "CHG"), -(ends_with("methylated")),  -chr, -pos, -context)
        data.4.80.meth.CHG<-select(filter(data.4.80, context == "CHG"), -(ends_with("total")),  -chr, -pos, -context)

        data.4.80.total.CG<-select(filter(data.4.80, context == "CG"), -(ends_with("methylated")),  -chr, -pos, -context)
        data.4.80.meth.CG<-select(filter(data.4.80, context == "CG"), -(ends_with("total")),  -chr, -pos, -context)

###### calculate methylation proportion by context
                 new.4.80.CHH<-data.4.80.meth.CHH/data.4.80.total.CHH
                 new.4.80.CHG<-data.4.80.meth.CHG/data.4.80.total.CHG
                 new.4.80.CG<-data.4.80.meth.CG/data.4.80.total.CG

###### transpose the data. We need samples in rows and C's in columns
                CHH.80.transpose<-as.data.frame(t(as.matrix(new.4.80.CHH)))
                CHG.80.transpose<-as.data.frame(t(as.matrix(new.4.80.CHG)))
                CG.80.transpose<-as.data.frame(t(as.matrix(new.4.80.CG)))

###### NOT GENERAL PART: extract names to plot with colors
                y.CHH.80<-rownames(CHH.80.transpose)
                label.CHH.80<-do.call(rbind, strsplit(y.CHH.80, '_'))
                CHH.80$Treat<-label.CHH.80[,2]
                CHH.80$Sample<-label.CHH.80[,1]

                y.CHG.80<-rownames(CHG.80.transpose)
                label.CHG.80<-do.call(rbind, strsplit(y.CHG.80, '_'))
                CHG.80$Treat<-label.CHG.80[,2]
                CHG.80$Sample<-label.CHG.80[,1]

                y.CG.80<-rownames(CG.80.transpose)
                label.CG.80<-do.call(rbind, strsplit(y.CG.80, '_'))
                CG.80$Treat<-label.CG.80[,2]
                CG.80$Sample<-label.CG.80[,1]
###### do PCA
                pca.CHH.80 <- PCA(CHH.80.transpose, graph=FALSE)
                pca.CHG.80 <- PCA(CHG.80.transpose, graph=FALSE)
                pca.CG.80 <- PCA(CG.80.transpose, graph=FALSE)
###### Save Whitout labels
                plot.CHH.80<-fviz_pca_ind(pca.CHH.80, geom.ind =c("point"),
                col.ind=CHH.80$Treat, repel=TRUE, title="CHH_80")

                tiff("CHH_80_NA0.tiff")
                print(plot.CHH.80)
                dev.off()

                plot.CHG.80<-fviz_pca_ind(pca.CHG.80, geom.ind =c("point"),
                col.ind=CHG.80$Treat, repel=TRUE, title="CHG_80")

                tiff("CHG_80_NA0.tiff")
                print(plot.CHG.80)
                dev.off()

                plot.CG.80<-fviz_pca_ind(pca.CG.80, geom.ind =c("point"),
                col.ind=CG.80$Treat, repel=TRUE, title="CG_80")

                tiff("CG_80_NA0.tiff")
                print(plot.CG.80)
                dev.off()

###### SAve With labels
                plot.CHH.80<-fviz_pca_ind(pca.CHH.80, geom.ind =c("point", "text"),
                col.ind=CHH.80$Treat, repel=TRUE, title="CHH_80")

                tiff("CHH_80_NA0_Labels.tiff")
                print(plot.CHH.80)
                dev.off()

                plot.CHG.80<-fviz_pca_ind(pca.CHG.80, geom.ind =c("point", "text"),
                col.ind=CHG.80$Treat, repel=TRUE, title="CHG_80")

                tiff("CHG_80_NA0_Labels.tiff")
                print(plot.CHG.80)
                dev.off()

                plot.CG.80<-fviz_pca_ind(pca.CG.80, geom.ind =c("point", "text"),
                col.ind=CG.80$Treat, repel=TRUE, title="CG_80")

                tiff("CG_80_NA0_Labels.tiff")
                print(plot.CG.80)
                dev.off()

###### ------ Histograms and bar plots ------------------------------------------------------------------
###### ------ Histogram for methylation proportion at each context
##### Stack methylated proportions to make the quantile distribution

                meth.CHH<-melt(new.3.80.CHH)
                meth.CHG<-melt(new.3.80.CHG)
                meth.CG<-melt(new.3.80.CG)

                qua.meth.CHH<-quantile(meth.CHH$value, c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,
                0.9,0.95,0.96,0.97,0.98, 0.99, 0.991, 0.995, 0.998,0.999,0.9991, 0.9995, 0.9999, 1), na.rm=TRUE)
                qua.meth.CHG<-quantile(meth.CHG$value, c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,
                0.9,0.95,0.96,0.97,0.98, 0.99, 0.991, 0.995, 0.998,0.999,0.9991, 0.9995, 0.9999, 1), na.rm=TRUE)
                qua.meth.CG<-quantile(meth.CG$value, c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,
                0.9,0.95,0.96,0.97,0.98, 0.99, 0.991, 0.995, 0.998,0.999,0.9991, 0.9995, 0.9999, 1), na.rm=TRUE)

#### Make the Histograms and save
                histo.CHH<-hist(meth.CHH$value)
                histo.CHG<-hist(meth.CHG$value)
                histo.CG<-hist(meth.CG$value)

                tiff("Histogram for CHH context.tiff")
                plot(histo.CHH)
                dev.off()


                tiff("Histogram for CHG context.tiff")
                plot(histo.CHG)
                dev.off()


                tiff("Histogram for CG context.tiff")
                plot(histo.CG)
                dev.off()

###### ------ Bar plot for mean methylated proportion at each context
##### Calculate Mean and SEM for methylated proportion at each context

                summ.CHH<-ddply(meth.CHH, c("variable"), summarise, Mean=mean(value, na.rm=TRUE), SE=sd(value, na.rm=TRUE),
                 SEM= sd(value, na.rm=TRUE)/sqrt(length(value)))
                summ.CHH<- summ.CHH %>% mutate(Treat= str_replace_all(variable, "_methylated", ""))
                summ.CHH<-separate(summ.CHH, Treat, c("Acc", "Treat"))
                summ.CHH<- arrange(summ.CHH, Acc, Treat)
                summ.CHH<- summ.CHH %>% mutate(Exp= str_replace_all(variable, "_methylated", ""))
                summ.CHH$Exp <- factor(summ.CHH$Exp, levels=unique(summ.CHH$Exp))
                head(summ.CHH)


                summ.CHG<-ddply(meth.CHG, c("variable"), summarise, Mean=mean(value, na.rm=TRUE), SE=sd(value, na.rm=TRUE),
                SEM= sd(value, na.rm=TRUE)/sqrt(length(value)))
                summ.CHG<- summ.CHG %>% mutate(Treat= str_replace_all(variable, "_methylated", ""))
                summ.CHG<-separate(summ.CHG, Treat, c("Acc", "Treat"))
                summ.CHG<- arrange(summ.CHG, Acc, Treat)
                summ.CHG<- summ.CHG %>% mutate(Exp= str_replace_all(variable, "_methylated", ""))
                summ.CHG$Exp <- factor(summ.CHG$Exp, levels=unique(summ.CHG$Exp))
                head(summ.CHG)


                summ.CG<-ddply(meth.CG, c("variable"), summarise, Mean=mean(value, na.rm=TRUE), SE=sd(value, na.rm=TRUE),
                SEM= sd(value, na.rm=TRUE)/sqrt(length(value)))
                summ.CG<- summ.CG %>% mutate(Treat= str_replace_all(variable, "_methylated", ""))
                summ.CG<-separate(summ.CG, Treat, c("Acc", "Treat"))
                summ.CG<- arrange(summ.CG, Acc, Treat)
                summ.CG<- summ.CG %>% mutate(Exp= str_replace_all(variable, "_methylated", ""))
                summ.CG$Exp <- factor(summ.CG$Exp, levels=unique(summ.CG$Exp))
                head(summ.CG)

##### Plot the bar

                limits.sd<-aes(ymax=Mean+SE, ymin=Mean-SE)
                limits.sem<-aes(ymax=Mean+SEM, ymin=Mean-SEM)

                tiff("Methylation proportion for CHH context.tiff")
                 ggplot(summ.CHH, aes(x=Exp, y=Mean)) +
                             geom_bar(stat="identity", width=0.75, color="pink", fill="pink", size=0.1)+
                             geom_errorbar(limits.sem, width=0.25)+
                             theme(axis.text.x=element_text(angle=90,hjust=1))
                dev.off()

                tiff("Methylation proportion for CHG context.tiff")
                ggplot(summ.CHG, aes(x=Exp, y=Mean)) +
                    geom_bar(stat="identity", width=0.75, color="pink", fill="pink", size=0.1)+
                    geom_errorbar(limits.sem, width=0.25)+
                    theme(axis.text.x=element_text(angle=90,hjust=1))
                dev.off()

                tiff("Methylation proportion for CG context.tiff")
                ggplot(summ.CG, aes(x=Exp, y=Mean)) +
                    geom_bar(stat="identity", width=0.75, color="pink", fill="pink", size=0.1)+
                    geom_errorbar(limits.sem, width=0.25)+
                    theme(axis.text.x=element_text(angle=90,hjust=1))
                dev.off()
