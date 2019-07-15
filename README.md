![logo](QBRC.jpg)
# <name>
## introduction
<name> is a preprocessing tool designed to convert scRNA seq results into a format compatible with analysis tools designed for general high-throughput sequencing data such as DNA seq. The main advantage of scRNA seq is that it allows researchers to sequence individual cells, whereas more traditional techniques sequence the aggregate genetic material from a population of cells. However, this means that scRNA seq results contain a mix of reads from perhaps thousands of diffrent cells. <name> seperates sequecing reads from scRNAseq by their cell of origin and performs a preliminary filter using bowtie2 alignment, producing high-quality inputs that is compatible with existing custom analysis tools.
## Getting started
### Installation
<name> is written in python can be installed from ...
#### System requirements
<name> requires a linux x86-64 installation (tested on RHEL 6, kernel 3.10.0-693)
Hardware requirements are dependent on CPU and input size
  - free drive space: 30x the size of compressed input fastqs, or 3x the size of uncompressed fastqs
  - CPU: no minimum requirement, but >16 core system with single core performance comparable to Intel e5-2680 is recommended
  - memory: >500 MB per logical processor, 1 GB per logical processor recommended
### Dependencies
python 3.6.4 or 3.7.4
bowtie2 2.3.2 or later
  
**python packages**

numpy, pandas

## Tutorial
This tutorial will guide the user in processing a sample scRNA seq
