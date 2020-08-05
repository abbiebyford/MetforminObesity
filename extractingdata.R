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
idatfiles = dir(pattern="idat")

#Accessed the HT-12 Human v4 bgx file from Ilumina Site https://support.illumina.com/downloads/humanht-12_v4_product_file.html
bgxfile = dir(pattern="bgx")