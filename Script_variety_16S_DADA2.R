
#                              Script DADA2  16S V3V4 Rice varieties microbiome 


#Load libraries

library(dada2)
library(vegan)
library(Biostrings)



setwd("path")

#starting the dada2 pipeline
path <- "path"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)



# Forward and reverse fastq filenames have format: SAMPLENAME-16S_1.fastq and SAMPLENAME-16S_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

#Filter and trim: Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter:  Macrogen 16S primers are 341F :  CCTACGGGNGGCWGCAG (17 nt) and  805R GACTACHVGGGTATCTAATCC (21 nt)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(17, 21), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
# On Windows set multithread=FALSE;


head(out) 
#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

#Remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted length. 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 382:432]

saveRDS(seqtab2, file="seqtab216SVarCry4help.RDS")


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
View(seqtab.nochim)


## [1]  20 232
sum(seqtab.nochim)/sum(seqtab2)
## [1] 0.964263

saveRDS(seqtab.nochim, file="seqtab.nochim-16SVarCry4help.RDS")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Stat_assemblyDADA2_16S_DivCry4help.csv")


#Assign taxonomy (download first silva_nr_v132_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v132.fa.gz from silva website)
taxa <- addSpecies(taxa, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_species_assignment_v138.fa.gz")


#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save seqtab and taxa as RDS files (easier to reload later)
saveRDS(taxa, file="taxa16S_DivVarieteCry4help_silva138.RDS")

#load metadata file
samdf <- read.csv2("Metadata_Varieties.csv")
samdf

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE (name check)

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Variety", "Condition", "Group") # a changer avec le nom de tes colonnes de m?tadonn?es
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE


# extract info in a single table (without metadata)
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa6<-as.data.frame(taxa)
taxa6$ID <- rownames(taxa)
all_data5  <- merge(taxa6,seqtab6,by='ID')
write.csv(all_data5, "alldata-VarieteyCry4help-16S_silva138.csv")


#phyloseq object without phylogenetic tree

ps_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_16S

