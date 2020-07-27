URL <- "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6418/files"
SDRF.file <- "E-MTAB-6418.sdrf.txt"
Data.file <- "E-MTAB-6418.raw.1.zip"
download.file(paste(URL,SDRF.file,sep="/"), SDRF.file)
download.file(paste(URL,Data.file,sep="/"), Data.file)
unzip(Data.file)
