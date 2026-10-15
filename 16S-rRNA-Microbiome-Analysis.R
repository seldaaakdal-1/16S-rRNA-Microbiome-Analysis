#################### 16S-rRNA-Microbiome-Analysis ##############################

# Installs and loads the packages required for 16S rRNA amplicon sequencing analysis.

install.packages(c("ggplot2", "vegan", "ape", "tidyverse", "reshape2", "phangorn", "knitr"))
BiocManager::install(c("phyloseq","dada2", "DESeq2", "ggtree", "DECIPHER", "microbiomeMarker"))

library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan) 
library(DESeq2)
library(ape)
library(tidyverse)
library(reshape2)
library(phangorn)
library(DECIPHER)
library(ggtree) 
library(microbiomeMarker)
library(knitr)

# Imports sample identifiers from a text file to define the input FASTQ files.

samples <- scan(file="samples.txt", what = "character")
samples

# Generates file paths for forward and reverse paired-end FASTQ reads.

forward.reads <- paste0(samples, "_L001_R1_001.fastq")
forward.reads

reverse.reads <- paste0(samples, "_L001_R2_001.fastq")
reverse.reads

# # Visualizes the quality scores of raw forward and reverse sequencing reads.

plotQualityProfile(forward.reads[1:4])
plotQualityProfile(reverse.reads[1:4])

# Defines output file names for quality-filtered forward and reverse reads.

flt.forward.reads <- paste0(samples, "_L001_R1_001_filtered.fastq")
flt.reverse.reads <- paste0(samples, "_L001_R2_001_filtered.fastq")

# Filters low-quality reads, removes PhiX contamination, and trims reads based on predefined quality thresholds.

out <- filterAndTrim(forward.reads, flt.forward.reads, reverse.reads, flt.reverse.reads,
                     maxEE = c(2,2), truncLen= c(240,160), rm.phix = TRUE)

out

# Evaluates the quality profiles of the filtered reverse reads.

plotQualityProfile(flt.reverse.reads[1:4])

# Learns sequencing error models from the filtered forward and reverse reads.

errF<- learnErrors(flt.forward.reads)
errR <- learnErrors(flt.reverse.reads)

# Visualizes the estimated sequencing error rates across quality scores.

plotErrors(errF, nominalQ = T)

# Infers exact amplicon sequence variants (ASVs) by correcting sequencing errors in the filtered reads.

dadaF <- dada(flt.forward.reads, err = errF, multithread = TRUE)
dadaR <- dada(flt.reverse.reads, err = errR, multithread = TRUE)

# Merges paired-end reads to reconstruct full-length amplicon sequences.

mergers <- mergePairs(dadaF, flt.forward.reads, dadaR, flt.reverse.reads)
mergers

# Creates an ASV count table containing the number of reads assigned to each unique sequence in each sample.

seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

# Removes chimeric sequences to improve the accuracy of downstream microbial community analyses. 

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)

# Extracts and ranks ASVs detected in the mock community sample.

unqs.mock <- seqtab.nochim["Mock_S280_L001_R1_001_filtered.fastq",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

# Compares inferred mock community sequences with expected reference sequences to evaluate sequence recovery.

mock.ref <- readDNAStringSet("HMP_MOCK.v35.fasta")
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# Retains the first 19 samples for downstream analyses.

seqtab.nochim <- seqtab.nochim[1:19,]

# Assigns taxonomic classifications to ASVs using the SILVA reference database.

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", multithread = T)

# Refines taxonomic assignments by matching ASVs against the SILVA species-level reference database.

taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa")

View(taxa)

# Saves the chimera-free ASV count table and taxonomic assignments as RDS files.

saveRDS(seqtab.nochim, "raw.seqtab.nochim.rds")
saveRDS(taxa, "raw.taxa.rds")

# Generates unique identifiers for the inferred ASVs.

asv_seqs <- colnames(seqtab.nochim)

asv_ids <- paste0(">ASV", seq(ncol(seqtab.nochim)))

asv_ids

# Exports ASV sequences in FASTA format for phylogenetic analysis and downstream applications.

fasta <- c(rbind(asv_ids, asv_seqs))

fasta

write(fasta, "ASVs.fa")

# Updates ASV identifiers in the count table and taxonomy table for consistency.

colnames(seqtab.nochim) <- asv_ids

rownames(taxa) <- asv_ids

# Assigns sample identifiers to the ASV table for compatibility with the metadata.

rownames(seqtab.nochim) <- c("F3D0", "F3D141", "F3D142", "F3D143", "F3D144", "F3D145", "F3D146", "F3D147", "F3D148", "F3D149", "F3D150", "F3D1", "F3D2", "F3D3", "F3D5", "F3D6", "F3D7", "F3D8", "F3D9")

# Saves the finalized ASV table and taxonomy data for subsequent microbiome analyses.

saveRDS(seqtab.nochim, "seqtab.nochim.rds")
saveRDS(taxa, "taxa.rds")

# Imports ASV sequences from the FASTA file.

seqs <- readDNAStringSet("ASVs.fa")

# Performs multiple sequence alignment of ASV sequences.

alignment <- AlignSeqs(seqs, anchor=NA)

# Converts aligned sequences into a phylogenetic data structure.

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# Calculates maximum-likelihood genetic distances between aligned sequences.

dm <- dist.ml(phang.align)

# Constructs a neighbor-joining phylogenetic tree based on the genetic distance matrix.

treeNJ <- NJ(dm)

# Roots the phylogenetic tree at its midpoint.

tree <- midpoint(treeNJ)

# Visualizes the phylogenetic relationships among ASVs using a circular tree layout.

ggtree(tree, layout="circular") + 
  geom_tiplab2(size=3, align=TRUE, linesize=0.3) +  
  ggtitle("Circular Phylogenetic Tree") +
  theme_minimal()

# Imports the ASV count table and taxonomic classification data.

seqtab.nochim <- readRDS("seqtab.nochim.rds")

taxa <- readRDS("taxa.rds")

# Removes FASTA-style prefix characters from ASV identifiers to ensure consistent naming.

colnames(seqtab.nochim) <- gsub(">", "", colnames(seqtab.nochim))

rownames(taxa) <- gsub(">","",rownames(taxa))

# Creates sample metadata containing sampling time points and corresponding collection days.

metadata <- data.frame(time=c("early", rep("late",10), rep("early", 8)),
                       day= c(0, 141:150,1:3, 5:9))

# Matches metadata row names to the sample identifiers in the ASV table.

rownames(metadata) <- rownames(seqtab.nochim) 

# Combines the ASV table, taxonomy, sample metadata, and phylogenetic tree into a phyloseq object.

physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                   tax_table(taxa),
                   sample_data(metadata),
                   phyloseq::phy_tree(tree))
physeq
physeq@otu_table

# Removes ASVs with a total abundance of 10 reads or fewer. 

physeq <- prune_taxa(taxa_sums(physeq) > 10, physeq)

# Excludes samples containing 1,000 reads or fewer.

physeq <- prune_samples(sample_sums(physeq) > 1000, physeq)

physeq

# Visualizes alpha diversity metrics across early and late sampling time points.

plot_richness(physeq = physeq, color = "time", title= "Alpha Diversity",)+
  geom_point(size=3)+
  theme(plot.title = element_text(hjust=0.5))

# Calculates alpha diversity indices to assess microbial richness and diversity within individual samples.

alpha_div <- estimate_richness(physeq)
alpha_div

# Exports alpha diversity measurements to a CSV file.

write.csv(alpha_div, "alpha_div.csv")

# Calculates Bray???Curtis dissimilarities to compare microbial community composition between samples.

bray_dist <-phyloseq::distance(physeq, method = "bray")

# Performs Principal Coordinates Analysis (PCoA) using Bray???Curtis distances.

ordination_bray <- ordinate(physeq, method = "PCoA", distance = bray_dist)

# Visualizes sample-level community differences using Bray???Curtis PCoA with group confidence ellipses.

plot_ordination(physeq, ordination_bray, color = "time", title = "PCoA Bray")+
  geom_point(size=3)+
  theme_bw()+
  stat_ellipse(aes(group=time), type = "t", size=0.5, linetype = "dashed")

# Calculates weighted UniFrac distances, accounting for phylogenetic relationships and relative microbial abundances.

wunifrac_dist <-phyloseq::distance(physeq, method = "wunifrac")

# Performs PCoA based on weighted UniFrac distances.

ordination_wunifrac <- ordinate(physeq, method = "PCoA", distance = wunifrac_dist)

# Visualizes phylogeny-aware beta diversity using weighted UniFrac PCoA.

plot_ordination(physeq, ordination_wunifrac, color = "time", title = "Weighted Unifrac PCoA")+
  geom_point(size=3)+
  theme_bw()+
  stat_ellipse(aes(group=time), type = "t", size=0.5, linetype = "dashed")

# Calculates unweighted UniFrac distances based on microbial presence or absence and phylogenetic relationships.

unifrac_unweighted <- phyloseq::distance(physeq, method = "uunifrac")

# Performs PCoA using unweighted UniFrac distances.

uunifrac_ordination <- ordinate(physeq, method = "PCoA", distance = unifrac_unweighted)

# Visualizes differences in microbial community membership using unweighted UniFrac PCoA.

plot_ordination(physeq, uunifrac_ordination, color = "time", title = "Unweighted Unifrac PCoA") +
  geom_point(size=4) +
  theme_bw() + coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = time), type = "t", size=0.5, linetype="dashed")

# Tests whether microbial community composition differs between time groups using PERMANOVA.

adonis_result <- adonis2(bray_dist ~ sample_data(physeq)$time, permutations=999)
adonis_result

# Evaluates differences in community composition between groups using ANOSIM.

anosim_result <- anosim(bray_dist, sample_data(physeq)$time, permutations=999)
anosim_result

# Assesses differences in within-group and between-group community dissimilarities using MRPP.

mrpp_result <- mrpp(bray_dist, sample_data(physeq)$time, permutations = 999)
mrpp_result

# Performs PERMANOVA on weighted UniFrac distances to test for differences between time groups.

weighted.unifrac.adonis <- adonis2(wunifrac_dist ~ sample_data(physeq)$time)
weighted.unifrac.adonis

# Applies ANOSIM to weighted UniFrac distances to evaluate group separation.

weighted.unifrac.anosim <- anosim(wunifrac_dist, sample_data(physeq)$time)
weighted.unifrac.anosim

# Applies MRPP to weighted UniFrac distances to assess differences between groups.

weighted.unifrac.mrpp <- mrpp(wunifrac_dist, sample_data(physeq)$time)
weighted.unifrac.mrpp

# Performs PERMANOVA on unweighted UniFrac distances to test for differences between time groups.

unweighted.unifrac.adonis <- adonis2(unifrac_unweighted ~ sample_data(physeq)$time)
unweighted.unifrac.adonis

# Applies ANOSIM to unweighted UniFrac distances to evaluate differences in microbial community composition between time groups.

unweighted.unifrac.anosim <- anosim(unifrac_unweighted, sample_data(physeq)$time)
unweighted.unifrac.anosim

# Applies MRPP to unweighted UniFrac distances to assess differences in microbial community composition between time groups.

unweighted.unifrac.mrpp <- mrpp(unifrac_unweighted, sample_data(physeq)$time)
unweighted.unifrac.mrpp

# Converts raw ASV counts into relative abundances for each sample.

ps <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

# Visualizes the relative abundance of microbial taxa at the phylum level across time groups.

plot_bar(ps, fill="Phylum")+
  facet_wrap(~time, scales="free_x")

# Visualizes the relative abundance of microbial taxa at the family level across time groups.

plot_bar(ps, fill="Family")+
  facet_wrap(~time, scales="free_x")

# Converts the phyloseq object into a DESeq2 dataset using time as the explanatory variable.

deseq_obj <- phyloseq_to_deseq2(physeq, ~time)

# Fits a negative binomial model to identify differentially abundant ASVs between time groups.

deseq_obj <- DESeq(deseq_obj)
 
# Extracts differential abundance results while disabling outlier filtering.

res <- results(deseq_obj, cooksCutoff = FALSE)
res

# Sets the adjusted p-value significance threshold.

alpha <- 0.05

# Selects statistically significant ASVs.

sigtab <- res[which(res$padj < alpha), ]

# Combines differential abundance results with taxonomic annotations.

sigtab <- cbind(as(sigtab, "data.frame"), 
                as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

head(sigtab)

# Exports significant results to a CSV file.

write.csv(sigtab, "sigtab.csv")

# Calculates the maximum log2 fold change for each phylum.

x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))

# Sorts phyla by their maximum log2 fold-change values.

x <- sort(x, TRUE)

# Orders phyla according to their maximum log2 fold-change values.

sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels = names(x))

# Visualizes log2 fold changes by genus, with points colored by phylum.

ggplot(sigtab, aes(Genus, log2FoldChange, colour = Phylum))+
  geom_point(size=5)+
  theme(axis.text.x = element_text(angle= -90, hjust=0, vjust=0.5))

# LEfSe (Linear Discriminant Analysis Effect Size) Analysis

# Removes taxa with missing phylum classifications.

physeq_filtered <- subset_taxa(physeq, !is.na(Phylum))

# Identifies potential microbial biomarkers at the phylum level between time groups.

lefse <- run_lefse(physeq_filtered, group= "time", taxa_rank = "Phylum")

# Visualizes LEfSe results using a bar plot.

plot_ef_bar(lefse)
