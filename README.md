![logo](QBRC.jpg)
# FASTR
## Introduction
FASTR (**F**astq **A**lignment-based **S**or**T**ing of sc**R**NA seq reads) is a preprocessing tool designed to convert scRNA seq results into a format suitable for mutation or SNV calling at the single cell level for the purpose of lineage tracing. The output is incidentally suitable for use with general high-throughput sequencing analysis tools to examine the transcriptome of each individual cell. The main advantage of scRNA seq is that it allows researchers to sequence individual cells, whereas more traditional techniques sequence the aggregate genetic material from a population of cells. However, this means that scRNA seq results contain a mix of reads from perhaps thousands of diffrent cells. FASTR seperates sequecing reads from scRNAseq by their cell of origin and performs preliminary QC using STAR alignment, producing high-quality inputs for finding rare mutations in individual cells.
## Getting started
### Installation
FASTR is written in python can be installed from its [GitHub page](https://github.com/zzhu33/test/blob/master/FASTR.zip). 
#### System requirements
FASTR requires a linux x86-64 operating system with basic utilities (split and gzip; tested on RHEL 6, kernel 3.10.0-693).

Hardware requirements are dependent on reference genome, CPU, and input size:
  - free drive space: 30x the size of compressed input fastqs, or 3x the size of uncompressed fastqs.
  - CPU: no minimum requirement, but >16 core system with single core performance comparable or better than Intel e5-2680 is   recommended.
  - memory: same as STAR (31 GB + 150 MB per logical processor for human genome index hg38) 64 GB recommended.
### Dependencies
[STAR](https://github.com/alexdobin/STAR) (tested using version 2.6.1b)

python 3.6.4 or 3.7.4
  
**python packages**

numpy, pandas


A STAR reference index is required. Refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for instructions on generating a genome index.

## Tutorial
This tutorial will guide the user in processing a sample scRNA seq result from 10x Genomic. 
### Inputs
The example inputs are [sample fastqs](http://cf.10xgenomics.com/samples/cell-exp/1.2.0/hgmm_100/hgmm_100_fastqs.tar) from 10x Genomics. The tarball needs to be extracted, although the individual files can be left as .fastq.gz files, or extracted fully as .fastq files. These files should appears similar to those shown in Fig.1.

**Fig.1:** Input fastq files.

![example_fastq](input_fastq.PNG)

Note that since FASTR is desinged to be easily integrated into other pipelines, filenames of inputs need to be stored in a tab-delimited text file, the included file is `exampleInputNames.txt`. The format of the input file is as follows:
```
lane1_r1.fastq  lane1_r2.fastq  lane1_I1.fastq
lane2_r1.fastq  lane2_r2.fastq  lane2_I1.fastq
...
```

The example reference genome is the [human GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) index. Users can use the following STAR command to generate the index. Refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for details.
```
STAR --runThreadN <logical cores> --runMode genomeGenerate --genomeDir <output path> --genomeFastaFiles <input fasta paths> --sjdbGTFfile <annotation files path> -- sjdbOverhang <readLength - 1>
```
In the example, the index is placed in `/home/STAR_indices/hg38/STAR`, other reference genomes such as mm10 can also be used.
### Running FASTR
Usage:
```
python3 FASTR.py [options]* --f <input_file> --i <input_directory> --r <output_directory> --ind <STAR_index_directory_path>
```
notes: 
  - syntax is not strictly enforced and inputs can be in any order that is suitable to the user.

**Example run:**

command:
```
$home/<name> python3 FASTR.py --cc 1000 --sc 100000 --f exampleInputNames.txt --i /home/FASTR_testing/fastqs --r /home/FASTR_testing/test_run --ind /home/STAR_indices/hg38/STAR
```

### Results
Results will be written to `/home/FASTR_testing/test_run` along with intermediate files and summary log. However, all intermediate files will be deleted. If the output directory does not exist, it will be automatically created.  The results will be in `/home/FASTR_testing/test_run/results`. By default, the results are compressed fastq files. There will be directories corresponding to each sample index that met the sample index read cuttoff, and numerous fastqs inside each directory, one for each cell barcode that met the specified cuttoff (Fig.2). A log ,`summary.txt`, is also provided, which includes run parameters and run times, as well as other information about the run. 

Note that 10x Genomics include four different sample indices in runs with only one sample. Here, each of the four sample indeces are treated as a seperate sample. Users can add `--ig True` option to the command or omit the I1 read filenames in the input file to 
**Fig.2** Output fastq files

![example output](output_fastq.PNG)

Note that there are four directories in the output (Fig.2, left), each with a different sample index and contains close to 100 cells (Fig.2, right). However, the [summary](http://cf.10xgenomics.com/samples/cell-exp/1.2.0/hgmm_100/hgmm_100_web_summary.html) provided by 10x Genomics indicates that there should be only ~100 valid cells present in the data.

A closer examination of the results would reveal that results from each sample appear to have nearly identical sets of cell barcodes. In this case, it can be useful to use the `--ig` option to check if the different sample indeces are truly what they seem:
```
$home/<name> python3 FASTR.py --ig True --cc 1000 --sc 100000 --n hgmm_100_S1_L001_I1_001.fastq --i /home/FASTR_testing/fastqs --r /home/FASTR_testing/test_run_pooled
```
The new outputs consists of 109 cells. Running FASTR with `--ig True` assumes that all reads are from the same sample and ignores any sample index. The summary log indicates that 88.52% of aligned reads were selected. This means that reads with the same cell and UMI barcodes but diferrent sample indeces aligned to the same positions. In this case, it appears safe to pool reads from different sample indeces. Conversely, if ~20% of aligned reads were selected, then it would be counterproductive to pool the reads.

It may seem odd that only 61.51% of reads are aligned to the reference. However, this sample data is a mix of mouse and human cells, and this tutorial only aligned to a human reference genome. When using the mouse mm10 reference index, 44.17% of reads successfully align. Note that there is overlap due to substantial similarities between the mouse and human genomes. FASTR currently only aligns to one genome for each run. FASTR still identified ~100 cells instead of ~50 due to genomes similarities, although fewer reads are expected in cells of mouse origin. The converse is true when aligning using the mouse mm10 reference index.



