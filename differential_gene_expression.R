#!/usr/bin/env Rscript

library(DESeq2)
library(readr)
library(dplyr)


suppressMessages(library(optparse))
option_list = list(
  make_option(c("-c", "--counts"), type="character", default="counts.csv", 
              help="Counts CSV file. Should have an initial 'locus_tag' column and then columns of the counts of reads mapping to each gene that are named by sample", metavar="counts.csv"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv", 
              help="Metadata CSV file. Should have an 'id' column with values that match the sample names in the Counts CSV file and a 'treatment' column with descriptions of the categories of different samples that will be compared. They will be compared in alphabetical order, so prefix with numbers if you want the differential gene expression log2 fold changes to be defined as 02_treatment_B versus 01_treatment_A, for example.", metavar="metadata.csv"),
  make_option(c("-g", "--genes"), type="character", default="genes.csv", 
              help="Gene information CSV file. Can be created by breseq CONVERT-REFERENCE -f CSV. Only lines that are marked as CDS, rRNA, tRNA, and ncRNA will be used.", metavar="genes.csv")
)

usage_string = paste(
  "differential_gene_expression.R -c counts.csv -m metadata.csv -g genes.csv",
  sep = ""
) 

opt_parser = OptionParser(usage=usage_string, option_list=option_list);
opt = parse_args(opt_parser);

# Load counts, the first column must be 'locus_tag'
counts = read.csv(opt$counts)
row.names(counts) = as.character(counts$locus_tag)
counts = counts %>% select(-locus_tag)

# Load metadata, it must have 'id' and 'treatment' columns 
metadata = read.csv(opt$metadata)
if (exists("metadata$biorep")) {
  metadata$biorep = as.factor(metadata$biorep)
}
metadata$treatment = factor(metadata$treatment)

# Load genes file
genedata = read.csv(opt$genes)
genedata = genedata %>% filter(feat_type %in% c("CDS", "tRNA", "rRNA", "ncRNA"))
genedata = genedata %>% rename(locus_tag = feat_id)


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design= ~ treatment)
dds <- DESeq(dds)

res <- results(dds)

res = as.data.frame(res)
res$locus_tag = row.names(res)

# Add gene information
res = res %>% left_join(genedata, by="locus_tag")

res = res %>% select(locus_tag, everything())

res$normalizedBaseMean = res$baseMean / (res$end - res$start + 1)

write.csv(res, "differential_expression.csv", row.names=F)
write_tsv(res, "differential_expression.tsv")


#Transform counts by both methods
vsd <- vst(dds, blind=FALSE)

vsd.to.write = assay(vsd)
row.names(vsd.to.write) = row.names(counts)
write.csv(vsd.to.write, "vst.normalized.counts.csv", row.names=T)

rld <- rlog(dds, blind=FALSE)

rld.to.write = assay(vsd)
row.names(rld.to.write) = row.names(counts)
write.csv(rld.to.write, "rlog.normalized.counts.csv", row.names=T)



