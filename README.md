![logo](QBRC.jpg)
# <name>
## Introduction
<name> is a preprocessing tool designed to convert scRNA seq results into a format compatible with analysis tools designed for general high-throughput sequencing data such as DNA seq. The main advantage of scRNA seq is that it allows researchers to sequence individual cells, whereas more traditional techniques sequence the aggregate genetic material from a population of cells. However, this means that scRNA seq results contain a mix of reads from perhaps thousands of diffrent cells. <name> seperates sequecing reads from scRNAseq by their cell of origin and performs a preliminary filter using bowtie2 alignment, producing high-quality inputs that is compatible with existing custom analysis tools.
## Getting started
### Installation
<name> is written in python can be installed from its [GitHub page](https://github.com/zzhu33/test).
#### System requirements
<name> requires a linux x86-64 operating system (tested on RHEL 6, kernel 3.10.0-693).

Hardware requirements are dependent on CPU and input size:
  - free drive space: 30x the size of compressed input fastqs, or 3x the size of uncompressed fastqs.
  - CPU: no minimum requirement, but >16 core system with single core performance comparable or better than Intel e5-2680 is   recommended.
  - memory: >500 MB per logical processor, 1 GB per logical processor recommended.
### Dependencies
bowtie2 2.3.2 or later

python 3.6.4 or 3.7.4
  
**python packages**

numpy, pandas


It is recommended to create a `bowtie2` directory in the installation directory and then dowload and extract the [human GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz) Bowtie2 reference index to the directory. This is the default reference index for alignment. Users can also specify any reference index explicitly with `--ind` using Bowtie2 syntax. Ex: `--ind /home/bowtie_inds/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index`.

## Tutorial
This tutorial will guide the user in processing a sample scRNA seq result from 10x Genomic. 
### Inputs
The example inputs are [sample fastqs](http://cf.10xgenomics.com/samples/cell-exp/1.2.0/hgmm_100/hgmm_100_fastqs.tar) from 10x Genomics. The tarball needs to be extracted, although the individual files can be left as .fastq.gz files, or extracted fully as .fastq files. 

The example reference index is the [human GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz) provided by [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Please refer to the Bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) for its installation instructions, if needed. This is the default Bowtie2 reference index and should be placed in `<name>/bowtie2`, as shown in the installation instructions. However, any index from any accessible directory can be used.
### Running <name>
Usage:
```
python3 <name>.py [options] --n <input_name> --i <input_directory> --r <output_directory> --ind <bowtie2_index>
```
note: syntax is not strictly enforced and inputs can be in any order that is suitable to the user.

**Example run:**

install numpy and pandas if needed:
```
pip install numpy
pip install pandas
```
run <name>:
```
$home/<name> python3 <name>.py --cc 1000 --sc 10000 --n hgmm_100_S1_L001_I1_001.fastq --i /home/<name>_testing/fastqs --r /home/<name>_testing/test_run 
```

If Bowtie2 index was not placed in the default location, add `--ind <path>/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index` to the run command.
### Results
Results will be written to `/home/<name>_testing/test_run` along with intermediate files and summary log. However, all intermediate files except sorted sams will be deleted. If the output directory does not exist, it will be automatically created.  The results will be in `/home/<name>_testing/test_run/results`. By default, the results are compressed fastq files. There will be directories corresponding to each sample index that met the sample index read cuttoff (10000 reads, `--sc 10000`), and numerous fastqs inside each directory, one for each cell barcode that met the specified cuttoff (1000 reads, `--cc 1000`).
A log ,`summary.txt`, is also provided, which includes run parameters and run times, as well as other information about the run. 
