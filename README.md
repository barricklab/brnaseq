# _brnaseq_
Basic analysis of bacterial RNAseq data for differential gene expression.


## Materials
* Read files
 * RNA-Seq data in FASTQ format (Ex: dataset1.fastq, dataset2.fastq)
* Reference files
 * Genome sequence in GenBank or GFF3 format (Ex: REL606.gbk)

## Software
* Adaptor filtering (Trimmomatic)
 * [Download Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 * Use this [protocol](https://barricklab.org/twiki/bin/view/Lab/ProtocolsTrimmomaticCommands).
* Read mapping (Bowtie2)
  * [Download Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2)
  * Download an executable for your platform
  * Add this directory to your $PATH or move bowtie, and bowtie-* star executable to your $PATH
* Sequence conversion and runfule generation (_breseq_ and _gdtools_)
  * [Download and install breseq](http://barricklab.org/breseq)
* Read counting (htseq)
  * Python3
  * Install using pip
* Differential gene expression (DEseq2)
  * [Download and install R](http://cran.r-project.org/)
  * Bioconductor R modules
    * library(deseq2)
* These _brnaseq_ scripts
  * Put them into your `$PATH`.
  * [Download from Github](https://github.com/barricklab/barricklab/blob/master/brnaseq).

## Procedure

### Create <code>genomediff</code> metadata files

The files should contain information about the reads and the references used in the analysis.

```
Example

```

### Trim adaptor sequences from reads

Use this [protocol](https://barricklab.org/twiki/bin/view/Lab/ProtocolsTrimmomaticCommands).

```
gdtools RUNFILE --mode  trimmomatic-PE-unique --preserve-pairs
```

### Run brnaseq

```
brnaseq
```

```
gdtools RUNFILE --preserve-pairs --executable brnaseq --options "-j 12 -k"
```

### Analyze differential gene expression

```
differential_gene_expression.R
```

For an explanation of the methods: [Manual and Instructions](http://bioconductor.org/packages/devel/bioc/html/DESeq.html)

### Graph data for a specific gene

```
graph_gene_counts.R
```

### Examine enrichment in certain gene categories

* MetaCyc
* [clusterProfiler R package](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

### View reads in IGV

You should always inspect your data at the nucleotide level to see if there are
indications that an analysis has gone awry. For example, are all of the forward
reads on the correct strand? Are there the mutations you expect in a certain strain there?

And convert to BAM format (assumes single-end data):
```
samtools faidx REL606.fna
samtools import REL606.fna datasetX.sam datasetX.unsorted.bam
samtools sort --threads 8 -o datasetX.bam datasetX.unsorted.bam
samtools index datasetX.bam
```

Now you can use IGV to view them!

### Plotting coverage

Extract R1, sort and index; Index FASTA reference
```
samtools faidx reference.fna
samtools view -hbf 64 aligned.paired.sam > unsorted_R1.bam
samtools sort --threads 8 -o sorted_R1.bam unsorted_R1.bam
samtools index datasetX.bam
```

Plot to table file
```
breseq BAM2COV -b sorted_R1.bam -f reference.fna -t -p 0 <seq_id:start-end>
```

## Citations

If you use this pipeline, you should cite:
1. trimmomatic
2. bowtie2
3. htseq
4. DESeq2
