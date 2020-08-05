URL <- "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6418/files"
SDRF.file <- "E-MTAB-6418.sdrf.txt"
Data.file <- "E-MTAB-6418.raw.1.zip"
download.file(paste(URL,SDRF.file,sep="/"), SDRF.file)
download.file(paste(URL,Data.file,sep="/"), Data.file)
unzip(Data.file)

#Read the sdrf file, this is an updated version I added birthweight centiles and LGA/AGA/SGA classes too.
SDRF <- read.delim("withcentiles.sdrf.txt", check.names=FALSE,stringsAsFactors=FALSE)

#Have a quick look at birthweight classes
SDRF[,c("Array Data File", "Class")]

#Set the birthweight class as a variable
Class <- SDRF[, "Class"]

#Set the different levels of class
levels <- c("LGA", "AGA", "SGA")
Class <- factor(Class, levels=levels)

str(SDRF)

#Extracting data from idat files
library(illuminaio)
#ie. read first file.
idat <- readIDAT("3998423016_A_Grn.idat")
names(idat)

#using limma instead
library(limma)
idatfiles = dir(pattern="idat")

#Accessed the HT-12 Human v4 bgx file from Ilumina Site https://support.illumina.com/downloads/humanht-12_v4_product_file.html
bgxfile = dir(pattern="bgx")

raw <- read.idat(idatfiles, bgxfile)

#plotting the probes signal intensity against negative control probes
pdf("boxplots_preNorm.pdf")
par(mfrow=c(1,2))
boxplot(log2(raw$E[raw$genes$Status=="regular",]),range=0, xlab="Arrays", ylab="log2 intensities", main= "Regular probes")
boxplot(log2(raw$E[raw$genes$Status=="negative",]),range=0, xlab="Arrays", ylab="log2 intensities", main= "Negative control probes")
dev.off()

#check for the proportion of probes on a microarray that correspond to the expressed genes 
#based on shi et al. 2010
propexp <- propexpr(raw)
propexp
dim(propexp) <- c(12,9) #set to this so I can see all 108 files clearly in View.
View(propexp)

mean(propexp)
summary(propexp)
summary(propexpr(raw))
#the proportion of probes varies between 0.3776 and 0.5869 (38% - 58%) mean is 53%
#lowest is 3998443071_F_Grn 33.776

#read the SDRF as a targets file
library(affy)
targets <- read.AnnotatedDataFrame(file="withcentiles.sdrf.txt", header=TRUE)
targets <- pData(targets)