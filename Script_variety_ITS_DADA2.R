#               Script DADA2 pipeline ITS2 Rice variety microbiome


#Load libraries

library(phyloseq)
library(DECIPHER)
library(dada2)
library(vegan)
library(Biostrings)


setwd("path")

#starting the dada2 pipeline
path <- "path"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


#DETECT AND REMOVE PRIMERS
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

#Filter:
#Macrogen ITS2 primers: ITS3 /  ITS4 : 3F - (20 bp) GCATCGATGAAGAACGCAGC; 4R - (20 bp)TCCTCCGCTTATTGATATGC

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(20, 20), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=TRUE)
# On Windows set multithread=FALSE; trimleft of 20 nt to remove primers; barcodes already removed in Macrogen data


head(out) 
#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, pool=FALSE, multithread=TRUE)

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

saveRDS(seqtab, "seqtabITSPathowater.RDS")

## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

#Remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted length. 
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 260:433]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## [1]  20 232
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.964263

saveRDS(seqtab.nochim, "seqtabnochim_ITS.RDS")

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


seqtab.nochim <- readRDS("seqtabnochim_ITS.RDS")

#Assign taxonomy (UNITE database)
taxa <- assignTaxonomy(seqtab.nochim, "////sh_general_release_dynamic_s_all_19.02.2025_dev.fasta", multithread=TRUE)


#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, "taxa_ITS_Pathowater.RDS")

View(sample.names)

#load metadata file
samdf <- read.csv2("C:///metadata_its_geneticdiv.csv")

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE (name check)

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Variety", "Group", "Origin", "Status", "Type") # a changer avec le nom de tes colonnes de m?tadonn?es
samdf <- samdf[rownames(seqtab.nochim), keep.cols]


# check correspondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

# extract info in a single table (without metadata)
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa7<-as.data.frame(taxa)
taxa7$ID <- rownames(taxa)
all_data6  <- merge(taxa7,seqtab6,by='ID')
write.csv(all_data6, "Alldata_ITS.csv")


#phyloseq object without phylogenetic tree

ps_its <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_its


