library(ggplot2)
library(tidyverse)

transformed.counts = read.csv("vst.normalized.counts.csv")
metadata = read.csv("EW_EL_strain_table.csv")
gene.info = read.csv("gene_info.csv")


#"tidy" the data
transformed.counts$locus_tag = transformed.counts$X
transformed.counts = transformed.counts %>% select(-X)
transformed.counts.tidy = transformed.counts %>% gather(id, "log2exp", -locus_tag)

# add strain data for grouping
transformed.counts.tidy = transformed.counts.tidy  %>% left_join(metadata %>% select(id, strain, genotype), by="id")

#add gene names to make selection easier
transformed.counts.tidy = transformed.counts.tidy  %>% left_join(gene.info %>% select(gene, locus_tag), by="locus_tag")

#create a table of means across replicates
transformed.counts.means.tidy = transformed.counts.tidy %>% group_by(strain, locus_tag) %>% summarize(mean.log2exp = mean(log2exp))
transformed.counts.means.tidy = transformed.counts.means.tidy  %>% left_join(metadata %>% select(strain, genotype), by="strain")
transformed.counts.means.tidy = transformed.counts.means.tidy  %>% left_join(gene.info %>% select(gene, locus_tag), by="locus_tag")

####### Repeat this part for any graph
## You can select different strains/genes for graphing

## fabA and fabB
locus_tags.of.interest = c("ECB_00958", "ECB_02248")

strains.of.interest = levels(transformed.counts.tidy$strain)
#uncomment to set to a subset
#strains.of.interest = c("DS1086", "DS1087", "DS1088", "DS1089", "DS1090")

counts.to.graph = transformed.counts.tidy
counts.to.graph = counts.to.graph  %>% filter(locus_tag %in% locus_tags.of.interest)
counts.to.graph = counts.to.graph  %>% filter(strain%in% strains.of.interest )

mean.counts.to.graph = transformed.counts.means.tidy
mean.counts.to.graph = mean.counts.to.graph  %>% filter(locus_tag %in% locus_tags.of.interest)
mean.counts.to.graph = mean.counts.to.graph  %>% filter(strain%in% strains.of.interest )

counts.to.graph$locus_tag.gene = paste0(counts.to.graph$locus_tag, " (", counts.to.graph$gene, ")")
mean.counts.to.graph$locus_tag.gene = paste0(mean.counts.to.graph$locus_tag, " (", mean.counts.to.graph$gene, ")")


ggplot() + geom_jitter(data=counts.to.graph, aes(x=strain, y=log2exp, color=locus_tag.gene), width=0.1) + geom_errorbar(data=mean.counts.to.graph, aes(x=strain, ymin=mean.log2exp, ymax=mean.log2exp, color=locus_tag.gene), width=.5)

