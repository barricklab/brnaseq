# brnaseq
Basic analysis of bacterial RNAseq data for differential gene expression

## Materials and Software

* Data files
 * RNA-Seq data in FASTQ format (Ex: dataset1.fastq, dataset2.fastq)
* Genome sequence files
 * Genome sequence in FASTA format (Ex: REL606.fna)
 * Genome sequence gene annotations in GFF3 format (Ex: REL606.gff3)
* Adaptor filtering software
 * [Download Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 * Use this [protocol](https://barricklab.org/twiki/bin/view/Lab/ProtocolsTrimmomaticCommands).
* Read mapping software
 * Bowtie2
  * [Download Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2)
  * Download an executable for your platform
  * Add this directory to your $PATH or move bowtie, and bowtie-* star executable to your $PATH
* _breseq_ and included _gdtools_
  * [Download and install breseq](http://barricklab.org/breseq)
* htseq python module and scripts
  * Install using pip
* R statistics package
  * [Download and install R](http://cran.r-project.org/)
* Bioconductor R modules
  * library(deseq2)
* <code>brnaseq</code> script
  * [Download from Github](https://github.com/barricklab/barricklab/blob/master/brnaseq).

## Procedure

### Create <code>genomediff</code> metadata files

The files should contain information about the reads and the

```
Example

```

### Trim adaptor sequences from reads

Use this [protocol](https://barricklab.org/twiki/bin/view/Lab/ProtocolsTrimmomaticCommands).

### Run brnaseq

```
brnaseq
```

### Analyze differential gene expression

```
differential_gene_expression.R
```

For an explanation of the methods [[http://bioconductor.org/packages/devel/bioc/html/DESeq.html][Manual and Instructions]]

### Graph data for a specific gene

```
graph_gene_counts.R
```

### Examine enrichment in certain gene categories

MetaCyc

### View reads in IGV

And convert to BAM format (assumes single-end data):  %BR%
<code>$samtools faidx REL606.fna </code> %BR%
<code>$samtools import REL606.fna datasetX.sam datasetX.unsorted.bam </code> %BR%
<code>$samtools sort datasetX.unsorted.bam datasetX </code> %BR%
<code>$samtools index datasetX.bam </code> %BR%

Now you can use IGV to view them.

## Citations

If you use this pipeline, you should cite:
1. trimmomatic
2. bowtie2
3. htseq
4. DESeq2
