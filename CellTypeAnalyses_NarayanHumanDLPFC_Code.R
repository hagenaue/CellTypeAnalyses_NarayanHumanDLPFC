#Narayan data download and analysis (GSE21138)
#Megan Hagenauer, 1/2017-7/2017
#Affymetrix Human Genome U133 Plus 2.0 Array

library(GEOquery)
gse <- getGEO("GSE21138", GSEMatrix = FALSE)
#lots of warnings. I wonder why.

head(Meta(gse))
str(Meta(gse))
Meta(GSMList(gse)$GSM528831)$characteristics_ch1
# [1] "brain region: BA46"                                                                                                  
# [2] "stage of illness [short doi=<5 yrs; intermediate doi=7-18yrs; long doi=>28 yrs]: short duration of illness - control"
# [3] "Sex: M"                                                                                                              
# [4] "age: 38"                                                                                                             
# [5] "tissue ph: 6.42"                                                                                                     
# [6] "pmi (hrs): 46"                                                                                                       
# [7] "type of drug: NA"                                                                                                    
# [8] "drug dose (chlorpromazine equivalents): NA" 

sub("Sex: ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[3])
as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[4]))
as.numeric(sub("pmi (hrs): ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[6], fixed = T))
as.numeric(sub("tissue ph: ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[5]))
sub("stage of illness [short doi=<5 yrs; intermediate doi=7-18yrs; long doi=>28 yrs]: ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[2], fixed=T)
sub("brain region: ","", Meta(GSMList(gse)$GSM528831)$characteristics_ch1[1])

SampleID<-as.matrix(names(GSMList(gse)))
Gender<-matrix("a", nrow=59, ncol=1)
Age<-matrix(0, nrow=59, ncol=1)
BrainpH<-matrix(0, nrow=59, ncol=1)
PMI<-matrix(0, nrow=59, ncol=1)
Diagnosis2<-matrix("a", nrow=59, ncol=1)
Tissue<-matrix("a", nrow=59, ncol=1)

length(GSMList(gse))
#[1] 59

for(i in c(1:59)){
  Gender[i]<-sub("Sex: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[3])
  Age[i]<-as.numeric(sub("age: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[4]))
  PMI[i]<-as.numeric(sub("pmi (hrs): ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[6], fixed = T))
  BrainpH[i]<-as.numeric(sub("tissue ph: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[5]))
  Diagnosis2[i]<-sub("stage of illness [short doi=<5 yrs; intermediate doi=7-18yrs; long doi=>28 yrs]: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2], fixed=T)
  Tissue[i]<-sub("brain region: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1])
}
Gender<-relevel(as.factor(Gender), ref="M")

Narayan_SampleCharacteristics<-data.frame(SampleID, Gender, Age, BrainpH, PMI, Diagnosis2, Tissue, stringsAsFactors=F)

head(Narayan_SampleCharacteristics)

Diagnosis<-c(rep("Control", 29), rep("Schiz", 30))
Diagnosis<-relevel(as.factor(Diagnosis), ref="Control")

Narayan_SampleCharacteristics<-data.frame(SampleID, Gender, Age, BrainpH, PMI, Diagnosis, Diagnosis2, Tissue, stringsAsFactors=F)

head(Narayan_SampleCharacteristics)

#Looks good.

setwd("~/Documents/Microarray Gen/Narayan_GSE21138")

write.csv(Narayan_SampleCharacteristics, "Narayan_SampleCharacteristics.csv")

#Now let's re-run RMA on their data:

library(org.Hs.eg.db)
library(plyr)
library(affy)

#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html

#cdf and the chip.

# hgu133plus2hsentrezg.db_19.0.2.tar.gz
# hgu133plus2hsentrezgcdf_19.0.0.tar.gz
# hgu133plus2hsentrezgprobe_19.0.0.tar.gz

#Already installed:
#install.packages(pkgs = c("hgu133plus2hsentrezg.db", "hgu133plus2hsentrezgcdf", "hgu133plus2hsentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

#Changed working directory to where the .cel files are located
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/GSE21138_RAW")

data2<-ReadAffy(cdfname ="hgu133plus2hsentrezg")
str(data2)
data2
# AffyBatch object
# size of arrays=1164x1164 features (41 kb)
# cdf=hgu133plus2hsentrezg (19764 affyids)
# number of samples=59
# number of genes=19764
# annotation=hgu133plus2hsentrezg
# notes=

eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDFplus2.txt")
RMAExpression_customCDFplus2<-read.delim("data_customCDFplus2.txt", sep="\t")
str(RMAExpression_customCDFplus2)
#'data.frame':	19764 obs. of  60 variables:
write.csv(RMAExpression_customCDFplus2, "RMAExpression_customCDFplus2.csv")


ScanDate<-protocolData(data2)$ScanDate
#Yep, there are definitely different scan dates here.
library(reshape2)
ScanDate_Split<-colsplit(ScanDate, pattern=" ", c("ScanDate", "ScanTime"))
table(ScanDate_Split$ScanDate)
04/27/06 04/28/06 05/17/05 05/19/05 06/14/06 06/22/05 
13       14       12       16        1        3
#Guess I should see if this matters...They are dates that are pretty close together, so unlikely to represent separate dissections.

#####################################################
## This is an alternative version of the analysis that I ran to see if controlling for RNA degradation improved the results (especially since the PMI is so long for these analyses). Controlling for RIN improves our analyses of the Illumina data, but those results are more sensitive to the reliability of individual probes, whereas Affy uses probesets.
source("https://bioconductor.org/biocLite.R")
biocLite("AffyRNADegradation")
library(AffyRNADegradation)

tongs <- GetTongs(data2, chip.idx = 4)
PlotTongs(tongs)
tongs <- GetTongs(data2, chip.idx = 5)
PlotTongs(tongs)

rna.deg<- RNADegradation(data2, location.type = "index")
RNADegradPerSample<-d(rna.deg)
str(RNADegradPerSample)

data3<-afbatch(rna.deg)
eset2 <- rma(data3)
write.exprs(eset2,file="data_customCDFplus2_RNADegCNTRL.txt")
RMAExpression_customCDFplus2_RNADegCNTRL<-read.delim("data_customCDFplus2_RNADegCNTRL.txt", sep="\t")
str(RMAExpression_customCDFplus2_RNADegCNTRL)
#'data.frame':	19764 obs. of  60 variables:
write.csv(RMAExpression_customCDFplus2_RNADegCNTRL, "RMAExpression_customCDFplus2_RNADegCNTRL.csv")
RMAExpression_customCDFplus2<-RMAExpression_customCDFplus2_RNADegCNTRL

#Alright - after snooping at the PCA vs. RNA degradation plots, it looks like they were actually aggravated by this correction. I suspect it is because we are using a custom .cdf and the correction that they are applying is based on the index of the probe within the probeset - i.e. probably references the original probeset and not the custom .cdf mapping of probe to transcript.
#Here's some quotes from the users manual to back up my suspicion:
# "Instead of using the probe index within the probeset as argument of the degra- dation degree, one can use the actual probe locations within the transcript. We have pre-computed the distance of each probe to the 3’ end of its target transcript for all Affymetrix 3’ expression arrays. These probe location files are available under the URL http://www.izbi.uni-leipzig.de/downloads_ links/programs/rna_integrity.php."
#"It is possible to use custom probe locations, for example if one wishes to analyze custom built microarrays or if one relies on alternative probe annotations."
#My thought: this sort of problem is likely to bias the results for any particular probeset, but the average amount of degradation across probesets per individual sample is likely to be unaltered - I think I should probably just include this as a covariate in my analyses instead of using their corrected Affybatch values.

#####################################################


head(RMAExpression_customCDFplus2)
RMAExpression_EntrezIDplus2<-sub("_at", "", RMAExpression_customCDFplus2[,1])
head(RMAExpression_EntrezIDplus2)
RMAExpression_customCDFAnnotationplus2<-data.frame(RMAExpression_customCDFplus2[,1], RMAExpression_EntrezIDplus2, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotationplus2)<-c("ProbesetID", "EntrezGeneID")

library(org.Hs.eg.db)
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 59887 

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2plus2<-join(RMAExpression_customCDFAnnotationplus2, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")


sum(is.na(RMAExpression_customCDFAnnotation2plus2[,3])==F)
#[1] 19528
dim(RMAExpression_customCDFAnnotation2plus2)
#[1] 19764     3
#So almost all of the results have gene symbols.

write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")

SignalSortedNoNA3<-as.matrix(RMAExpression_customCDFplus2[,-1])

#Double-checking that the Sample info and Expression data are in the same order:
cbind(SampleID, colnames(SignalSortedNoNA3))
#Yep, they are both in ascending numeric order.

Narayan_SampleCharacteristics<-data.frame(Narayan_SampleCharacteristics, ScanDate_Split)


####Quality Control################

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",]

RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]

str(RMAExpression_customCDFplus2)

setwd("~/Documents/Microarray Gen/Narayan_GSE21138")
png("Boxplot_RMAExpression_customCDFplus2.png", width=2000, height=400)
boxplot(SignalSortedNoNA3)
dev.off()
#Looks quantile normalized and the signal values now have an appropriate range. Hurray!

png("XIST_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]~Gender, col=2)
dev.off()
#Ooh - there are a couple of gender switches.
#let's double check that with another gene

png("RPS4Y1_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1",][169,]~Gender, col=2)
dev.off()
#At least one subject still looks mistaken here.

#Let's double check one more time:

png("DDX3Y_vs_Gender_customCDFplus2.png")
boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y",][175,]~Gender, col=2)
dev.off()
#Yep, there is still one mistaken subject

png("RPS4Y1vsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="RPS4Y1",][169,]~SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,], col=as.numeric(Gender)+1)
dev.off()

png("DDX3YvsXIST_GenderCheck.png")
plot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="DDX3Y",][175,]~SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,], col=as.numeric(Gender)+1)
dev.off()


#Interesting. Two subjects with really high XIST still have really high Y chromosome expression. XIST is expressed in males during spermatogenesis...but in the brain?
#Or could it be a female with a y chromosome? XXY is klinefelter syndrome
#Wikipedia: "The primary feature is sterility.[1] Often symptoms may be subtle and many people do not realize they are affected. Sometimes symptoms are more prominent and may include weaker muscles, greater height, poor coordination, less body hair, smaller genitals, breast growth, and less interest in sex.[2] Often it is only at puberty that these symptoms are noticed.[3] Intelligence is usually normal; however, reading difficulties and problems with speech are more common."
#"Klinefelter syndrome is one of the most common chromosomal disorders, occurring in 1:500 to 1:1000 live male births."
#About 64% of affected individuals are never recognized.[35] 

#Klinefelter’s syndrome (karyotype 47,XXY) and schizophrenia-spectrum pathology
#SOPHIE VAN RIJN, ANDRÉ ALEMAN, HANNA SWAAB, RENÉ KAHN
#The British Journal of Psychiatry Oct 2006, 189 (5) 459-461; DOI: 10.1192/bjp.bp.105.008961
#"Klinefelter’s syndrome, characterised by a 47,XXY chromosomal pattern, has largely been associated with physical abnormalities. Here, we report high levels of schizophrenia-spectrum pathology in 32 men with this syndrome in comparison with 26 healthy controls. This may have implications for treatment of have implications for treatment of Klinefelter’s syndrome and suggests that the X chromosome may be involved in the aetiology of schizophrenia."

Narayan_SampleCharacteristics[Narayan_SampleCharacteristics$Gender=="M" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]>7,]

#Ah - but the subjects do not have Schiz:
#SampleID Gender Age BrainpH PMI Diagnosis                                 Diagnosis2 Tissue
#9  GSM528839      M  38    6.19  44   Control intermediate duration of illness - control   BA46
#10 GSM528840      M  42    6.21  26   Control intermediate duration of illness - control   BA46


Narayan_SampleCharacteristics[Narayan_SampleCharacteristics$Gender=="F" & SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]<7,]

#SampleID Gender Age BrainpH PMI Diagnosis                                       Diagnosis2
#50 GSM528880      F  33    6.43  20     Schiz intermediate duration of illness - schizophrenia
#Tissue
#50   BA46

#Note - none of these are the outliers that were actually removed in the original study. 
#Interesting. 

#Let's check PCA and see if there are any other obvious problem samples:

################################
# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/PCA")


pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, RMAExpression_customCDFAnnotation2plus2)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

# #Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()


# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()
#Yes, there are three outliers identified through PCA too.


# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()


SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-cbind(SubjectFactorVariables, SubjectContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#######################################

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(SignalSortedNoNA3), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=1000, height=400)
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()
#Outliers: GSM528866, GSM528873
#The same outliers are present after correcting for RNADegradation

#Unfortunately, those don't overlap with the ones that failed gender check, so our n is going to drop by 5.

#Outlier removal: 

SampleCharacteristics_NoOutliers<-Narayan_SampleCharacteristics[(SampleID%in%c("GSM528866", "GSM528873", "GSM528880", "GSM528840", "GSM528839"))==F,]
dim(SampleCharacteristics_NoOutliers)
#[1] 54 10

SignalSortedNoNA3NoOutliers<-SignalSortedNoNA3[,(SampleID%in%c("GSM528866", "GSM528873", "GSM528880", "GSM528840", "GSM528839"))==F]
dim(SignalSortedNoNA3NoOutliers)
#[1] 19764    54

RNADegradPerSampleNoOutliers<-RNADegradPerSample[(SampleID%in%c("GSM528866", "GSM528873", "GSM528880", "GSM528840", "GSM528839"))==F]

png("RNADegradPerSampleNoOutliersVsPMI.png")
plot(RNADegradPerSampleNoOutliers~PMI)
dev.off()
#Interesting - the correlation with PMI is pretty loose.

#Redefining the variables to only include good data:
SignalSortedNoNA3<-SignalSortedNoNA3NoOutliers

Gender<-SampleCharacteristics_NoOutliers$Gender
Age<-SampleCharacteristics_NoOutliers$Age
PMI<-SampleCharacteristics_NoOutliers$PMI
BrainpH<-SampleCharacteristics_NoOutliers$BrainpH
Diagnosis<-SampleCharacteristics_NoOutliers$Diagnosis
ScanDate<-SampleCharacteristics_NoOutliers$ScanDate
ScanDate<-as.factor(ScanDate)

#############################################

#Characterizing subjects and looking for relationships between subject variables:
#Changed Working directory
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/SubjVarVsSubjVar")

SubjectFactorVariables<-cbind(Diagnosis, Gender, ScanDate)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "ScanDate")

SubjectContinuousVariables<-cbind(BrainpH, PMI, Age, RNADegradPerSampleNoOutliers)
colnames(SubjectContinuousVariables)<-c("BrainpH", "PMI", "Age", "RNADegradPerSample")

for (i in 1:length(SubjectContinuousVariables[1,])){
  png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
  hist(SubjectContinuousVariables[, i], col=i+1)
  dev.off()		
}
#Wide age distribution, but many of the subjects are younger (~30 years old)
#wide pH distribution (5.6-6.8), but typically low (~6.3), two moderate outliers with 5.6 and 5.8
#PMI is almost always longer than 40 hrs


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}


#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}


#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("14 Cross Tabs Between Subject Factors.txt")
out<-c(
  capture.output(
    
    summary(Diagnosis),
    summary(Gender),
    
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(ContingencyTable)
        print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
      }		
    }
  )
)
cat(out, file="14 Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)


library(car)

StatisticalRelationshipsIV<-file("14 Statistical Relationships between Subject Variables.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMI+Diagnosis+Gender+Age))
    
  ),
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(chisq.test(ContingencyTable))		
      }		
    }
  )
  
)
cat(out, file="14 Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
        if(chisq.test(ContingencyTable)$p.value<0.05){
          print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))	
        }else{}		
      }		
    }
  )
)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)

############################################


#Re-Ran PCA, then looked at relationships with subject variables.

setwd("~/Documents/Microarray Gen/Narayan_GSE21138/PCAVsSubjVar")

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
    RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}



#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}




#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
    summary.lm(lm(PC1noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSampleNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC2noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSampleNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC3noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSampleNoOutliers))
  ),
  
  capture.output(
    summary.lm(lm(PC4noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSampleNoOutliers))
  ),
  
  
  #Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between PCA and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  )
  
)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)



#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  )
)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)

#Holy crap - Scan date really really  matters.
#RNADegradation matters quite a bit too.
#RNADegradation *still* matters even after correcting for it within the AffyBatch - perhaps even more.

#####################

#Running Cell Type analyses:
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/BrainInABlender")

install.packages("devtools")

library("devtools")

install_github("hagenaue/BrainInABlender")

library("BrainInABlender")

temp<-data.frame(as.character(RMAExpression_customCDFAnnotation2plus2[,3]), SignalSortedNoNA3, stringsAsFactors=F)
write.csv(temp, "UserInput.csv")
str(temp)
table(table(RMAExpression_customCDFAnnotation2plus2[,3]))
SirUnMixALotOutput<-Sir_UnMixALot(userInput=temp, dataColumns=c(2:55), geneColumn=1, species="human")

PublicationSpecific_CellTypeIndex<-SirUnMixALotOutput$PublicationSpecific_CellTypeIndex
AveragePrimary_CellTypeIndex<-SirUnMixALotOutput$AveragePrimary_CellTypeIndex

str(PublicationSpecific_CellTypeIndex)

#Wow- the output is really beautiful. Beautiful clusters.


temp<-cbind(as.matrix(SubjectContinuousVariables), t(PublicationSpecific_CellTypeIndex), t(AveragePrimary_CellTypeIndex)) 
str(temp)
SubjectContinuousVariables<-temp

########################################################################
#Then I re-ran the subjvar and PCA code.
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/CellTypeVsSubjVar")


row.names(AveragePrimary_CellTypeIndex)
for(i in c(1:10)){
  print(" 
        ")
  print(row.names(AveragePrimary_CellTypeIndex)[i])
  print(summary.lm(lm(AveragePrimary_CellTypeIndex[i,]~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSampleNoOutliers)))
}


setwd("~/Documents/Microarray Gen/Narayan_GSE21138/PCA")

###############################

setwd("~/Documents/Microarray Gen/Narayan_GSE21138/LM4_Basic")
#Later, I re-ran LM4_Basic because for some reason the results of LM4_Basic were bizarrely stronger than any other model
#setwd("~/Documents/Microarray Gen/Narayan_GSE21138/RNADegCNTRL/LM4_BasicReRun")
#...but the results look exactly the same.

#I need to output results from a basic linear model for comparison:

BrainpHcentered<-BrainpH-median(BrainpH)
Agecentered<-Age-median(Age)
PMIcentered<-PMI-median(PMI)

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 11)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 11)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+ScanDate+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}

GeneByCellTypeSubjVar2_Pvalues2<-cbind(GeneByCellTypeSubjVar2_Pvalues, RMAExpression_customCDFAnnotation2plus2)
GeneByCellTypeSubjVar2_Betas2<-cbind(GeneByCellTypeSubjVar2_Betas, RMAExpression_customCDFAnnotation2plus2)

write.csv(GeneByCellTypeSubjVar2_Pvalues2, "GeneByCellTypeSubjVar2_Pvalues.csv")
write.csv(GeneByCellTypeSubjVar2_Betas2, "GeneByCellTypeSubjVar2_Betas.csv")


GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 11)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")

for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+ScanDate+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  
}

GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")



for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), "png", sep="."))	
  hist(GeneByCellTypeSubjVar2_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
  abline(a=(length(GeneByCellTypeSubjVar2_Pvalues[,1])/100), b=0)
  dev.off()		
}	


GeneByCellTypeSubjVar2_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar2_Pvalues[,1]), length(GeneByCellTypeSubjVar2_Pvalues[1,]))
colnames(GeneByCellTypeSubjVar2_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar2_Pvalues)
row.names(GeneByCellTypeSubjVar2_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar2_Pvalues)


library(multtest)

for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
  
  #Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
  TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar2_Pvalues[,i], proc=c("BH"))
  GeneByCellTypeSubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
  
}

GeneByCellTypeSubjVar2_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar2_PvaluesAdj, RMAExpression_customCDFAnnotation2plus2)
write.csv(GeneByCellTypeSubjVar2_PvaluesAdj2, "GeneByCellTypeSubjVar2_PvaluesAdj.csv")

GeneByCellTypeSubjVar2DF<-as.data.frame(cbind(GeneByCellTypeSubjVar2_Betas, GeneByCellTypeSubjVar2_Pvalues, GeneByCellTypeSubjVar2_PvaluesAdj))

temp<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2DF)
write.csv(temp, "GeneByCellTypeSubjVar2DF.csv" )

#Batch correction seemed to help pull out Schiz signatures. :)
#RNA degradation also really matters.
#Controlling for RNA degradation seemed to help draw out recognizable Schiz signatures too. :)

########################

setwd("~/Documents/Microarray Gen/Narayan_GSE21138/LM4_Prevelent")

row.names(AveragePrimary_CellTypeIndex)
Astrocyte_All<-AveragePrimary_CellTypeIndex[1,]
Oligodendrocyte_All<-AveragePrimary_CellTypeIndex[8,]
Microglia_All<-AveragePrimary_CellTypeIndex[3,]
Endothelial_All<-AveragePrimary_CellTypeIndex[2,]
RBC_All<-AveragePrimary_CellTypeIndex[10,]
Neuron_All<-AveragePrimary_CellTypeIndex[5,]
Neuron_Interneuron<-AveragePrimary_CellTypeIndex[6,]
Neuron_Projection<-AveragePrimary_CellTypeIndex[7,]
Oligodendrocyte_Immature<-AveragePrimary_CellTypeIndex[9,]
Mural_All<-AveragePrimary_CellTypeIndex[4,]

summary.lm(lm(PC1noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+ScanDate+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC2noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+ScanDate+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC3noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+ScanDate+RNADegradPerSampleNoOutliers))

summary.lm(lm(PC4noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+ScanDate+RNADegradPerSampleNoOutliers))

setwd("~/Documents/Microarray Gen/Narayan_GSE21138/LM4_Prevelent")

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)

GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")


for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+ScanDate+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}



GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")





###############
setwd("~/Documents/Microarray Gen/Narayan_GSE21138/LM4_Everything")
#I don't know if the sample size for this dataset can actually handle this.

GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 20)
GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 20)
colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC","ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")
row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
head(GeneByCellTypeSubjVar2_Pvalues)

GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 20)
colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "ScanDate2", "ScanDate3", "ScanDate4", "ScanDate5", "RNADeg")

for(i in c(1:length(SignalSortedNoNA3[,1]))){
  
  temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Endothelial_All+Microglia_All+Mural_All+Neuron_All+Neuron_Projection+Neuron_Interneuron+Oligodendrocyte_All+RBC_All+ScanDate+RNADegradPerSampleNoOutliers))
  
  GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
  GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
  GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
  
}




GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)

write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")





###############
