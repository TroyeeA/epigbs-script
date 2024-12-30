#!/home/anneaa/miniconda3/envs/assay/bin/Rscript

#-------------------- Parse arguments
argsis<-commandArgs(trailingOnly=TRUE)
directory_is=argsis[1]		#smps/gene_region/		
file_is=argsis[2]		#B15_d5D32_gene_inputDSS.cpgf

				#directory_is="/faststorage/project/spider2/assay_study_anne/methylation/smps/gene_region/"
				#file_is="S29_d5D32_gene_inputDSS.cpgf"


#-------------------- Setup
library(DSS)
require(bsseq)



#-------------------- Parse file name
	# Get the name that is shared for all sample files:
patternis=paste(sapply(strsplit(file_is,split="_"),"[[",2),sapply(strsplit(file_is,split="_"),"[[",3),sapply(strsplit(file_is,split="_"),"[[",4),sep="_")
	#patternis="d5D32_gene_inputDSS.cpg"
	# Adjust to your needs
files_are <- dir(directory_is,pattern=paste0(patternis,"$"),full.names=T) # "$" makes sure the string ends
	# split for output naming:
outname=sapply(strsplit(patternis,"_input"),"[[",1)
	# d5D32_gene

#-------------------- Import files into list of dataframes
samplelist=vector()
varlist=list()
for (each in files_are){
	message(each)
	sample_name=sapply(strsplit(each,"//|/|_"),"[[",12) #Gets the sample name (uniq) --> Adjust to your needs
	nam <- paste("data_", sample_name, sep = "")	# assigns a name to a variable for use in next step
	assign(nam,read.table(each,header=T)) # this assigns the read-in table to a variable with the name that is inside the nam variable - equivalent to data_POPtreat
	varlist[[nam]]<- get(nam)	# get the data now assigned to nam and add to list varlist
	samplelist <- append(samplelist,sample_name) # add also the sample name to a vector
	rm(nam)
}

#-------------------- Make design matrix

poplist<-sapply(strsplit(samplelist,"[0-9]"),"[[",1) # split samplenames to pop
templist<-sapply(strsplit(samplelist,"[A-Z]"),"[[",2) # split samplenames to temp
design<-as.data.frame(cbind(poplist,templist)); colnames(design)<- c("population","Temperature") # make design matrix


#-------------------- Make BSobject needed for DSS
BSobj <- makeBSseqData( varlist, samplelist)
BSobj
'
An object of type BSseq with
  177134 methylation loci
  20 samples
has not been smoothed
All assays are in-memory
'


#-------------------- Fit the model in DSS

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~population+Temperature)
	# Here you Could add the interaction: +pop:temp


#-------------------- Test the effect of each term in the model and write out data
		#### First term:
term="population"  # adjust the term to fit your data
test_1 = DMLtest.multiFactor(DMLfit, term=term) 
head(test_1)
		##       chr     pos     stat        pvals         fdrs
		## 1273 chr1 2930315 5.280301 1.289720e-07 0.0006448599
		## 4706 chr1 3321251 5.037839 4.708164e-07 0.0011770409

	# write all results to file
write.table(test_1, file=paste0(directory_is,"/DSS_result_addit_",outname,"_",term,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	# pick out results that are significant and below FDR 0.1
test_1_sub <- test_1[ which(test_1$pvals < 0.05 & test_1$fdrs < 0.1),]
write.table(test_1_sub, file=paste0(directory_is,"/DSS_result_addit_",outname,"_",term,"_significant.txt"),row.names=F,col.names=T,quote=F,sep="\t")


		#### Second term:
term="Temperature"  # adjust the term to fit your data
test_2 = DMLtest.multiFactor(DMLfit, term=term)
head(test_2)
	# write all results to file
write.table(test_2, file=paste0(directory_is,"/DSS_result_addit_",outname,"_",term,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	# pick out results that are significant and below FDR 0.1
test_2_sub <- test_2[ which(test_2$pvals < 0.05 & test_2$fdrs < 0.1),]
write.table(test_2_sub, file=paste0(directory_is,"/DSS_result_addit_",outname,"_",term,"_significant.txt"),row.names=F,col.names=T,quote=F,sep="\t")


		## Combine terms table
require("plyr")
test_1
test_2

### OBS fix this ;)

merged_table <- join_all(list(effect_of_temp,effect_of_pop,effect_of_interact), by ="rownamesis" , type = 'left')
colnames(merged_table)<-c("region_name","stat_T","pvalue_T","padj_T","stat_P","pvalue_P","padj_P","stat_I","pvalue_I","padj_I")

effect_of_temp <- as.data.frame(read.table(effect_of_temp_file,header=T,stringsAsFactors=T));
effect_of_temp <- effect_of_temp[,4:6];
effect_of_temp<-cbind(rownames(effect_of_temp),effect_of_temp);
colnames(effect_of_temp)[1]<-"rownamesis"


#Maybe call DMRs:
	#**Step 4*. DMRs for multifactor design can be called using {callDMR} function:

#callDMR(test_pop, p.threshold=0.05)
	##      chr   start     end length nCG   areaStat
	## 33  chr1 2793724 2793907    184   5  12.619968
	## 413 chr1 3309867 3310133    267   7 -12.093850




