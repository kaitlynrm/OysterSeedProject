

path <- "Documents/robertslab/gannet.fish.washington.edu/Atumefaciens/20190225_cpg_oe" #defining master directory
folder.list <- dir(path) #making list of subdirectories w/in master directory

d <- read.table(paste(path,folder.list[1], "ID_CpG",sep="/"), sep="\t", stringsAsFactors=FALSE) #reading in ID_CpG file in first subdirectory
colnames(d) <- c("ID", "CL_1") #renaming columns of file just loaded in 
for(i in 2:length(folder.list)){
  dtemp <- read.table(paste(path,folder.list[i], "ID_CpG",sep="/"), sep="\t", stringsAsFactors=FALSE) #starting loop for ID_CpG files in all subdirectories
  colnames(dtemp)[1] <- "ID" #rename col 1 
  colnames(dtemp)[2] <- sub(".*ANACfill.","",folder.list[i]) #sub all characters before (and including) ANACfill. with NOTHING ("") in folder (i)
  colnames(dtemp)[2] <- sub("_GENE_analysis","",colnames(dtemp)[2]) #sub GENE... with nothing in folder (i)
  d <- merge(d,dtemp, by = "ID") #merge file in subdirectory with previously created file
}

write.csv(d, "Documents/robertslab/gannet.fish.washington.edu/Atumefaciens/20190225_cpg_oe/CpG_oe-gene_analysis")
