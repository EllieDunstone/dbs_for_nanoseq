# Ellie Dunstone, 2023

# This script annotates dinucleotides and runs of consecutive subs of any length using a function written by Iñigo Martincorena.
# This function takes a data frame of mutations with the following columns: 
# sampleID, chr = chromosome, pos = mutated site, ref = reference base, mut = mutant base.


# NB. For users at the Wellcome Sanger Institute, the wrapper script (dinuc_nanoseq_caller_wrapper.sh) can be used to pull mutation files from the compute farm, filter them, and format them correctly to be used in this script. That script also acts as a wrapper for this script.

library(tidyverse)
library(vcfR)

#list file names in the area
files <- list.files(path = "/Users/ed4/Documents/temp/dbs_call_test/", pattern = "*.tsv$")

#read in as tables
for (file in files) {
  temp <- read.table(paste0("/Users/ed4/Documents/temp/dbs_call_test/", file), sep="\t", header = TRUE)
  assign(file, temp)
}

#add sample ID as column
tables <- c()
for (file in files) {
  new <- mutate(get(file), sampleID = paste0(strsplit(file, "_")[[1]][1], "_", strsplit(file, "_")[[1]][2]))
  assign(paste0("sampleID_", file), new)
  tables <- cbind(tables, paste0("sampleID_", file))
}

#concatenate tables
mut_table <- data.frame(matrix(ncol=ncol(new))) #init table
colnames(mut_table) <- colnames(new)

for (table in tables) {
  mut_table <- rbind(mut_table, get(table))
  print(table)
}

mut_table <- mut_table[-1,] #drop blank first row


#grab cols you need
mut_table <- select(mut_table, c(sampleID, chrom, chromStart, context, call))
mut_table <- mutate(mut_table, "context" = substr(context, 2, 2))
colnames(mut_table) <- c("sampleID", "chr", "pos", "ref", "mut")
  

#DBS calling algorithm from Inigo Martincorena (2022) - define function

annotruns = function(mutations) {
  mutations = mutations[order(mutations$sampleID, mutations$chr, mutations$pos), ]
  d = mutations$pos-(1:nrow(mutations))
  runs = rle(d)
  rmpos = rep(0,nrow(mutations))
  runstarts = cumsum(runs$length)-runs$length+1
  for (h in 1:length(runs$length)) {
    if (runs$length[h]>1) { # Adjacent mutations
      mutcluster = runstarts[h]:(runstarts[h]+runs$lengths[h]-1)
      rmpos[mutcluster[-1]] = 1 # Removing all the affected rows except the first one (which we will edit to capture the complex event)
      mutations[mutcluster[1],"ref"] = paste(mutations[mutcluster,"ref"],collapse="")
      mutations[mutcluster[1],"mut"] = paste(mutations[mutcluster,"mut"],collapse="")
    }
  }
  mutations = mutations[!rmpos,]
  return(mutations)
}

# run DBS calling
collapsed_table <- annotruns(mut_table) #this now shows the dbs separately

#filter for just the dbs
dbs_table <- filter(collapsed_table, nchar(ref) == 2)

#filter for longer variants
mnv_table <- filter(collapsed_table, nchar(ref) > 2)

#write out
write_csv(dbs_table, "/Users/ed4/Documents/temp/dbs_call_test/dbs_calls.csv")
write_csv(mnv_table, "/Users/ed4/Documents/temp/dbs_call_test/mnv_calls.csv")

#make input for SigProfilerMatrixGenerator
dbs_table_txt <- dbs_table %>% mutate(Project="Dinucs") %>% mutate(ID=".") %>% mutate(Genome="GRCh37") %>% mutate(mut_type="DNP") %>% mutate(pos_end=pos+1) %>% mutate(Type="SOMATIC")
colnames(dbs_table_txt) <- c("Sample", "chrom", "pos_start", "ref", "alt", "Project", "ID", "Genome", "mut_type", "pos_end", "Type")
dbs_table_txt_final <- dbs_table_txt[,c(6,1,7,8,9,2,3,10,4,5,11)]

#output
write_tsv(dbs_table_txt_final, "/Users/ed4/Documents/temp/dbs_call_test/dbs_calls_text.txt")
