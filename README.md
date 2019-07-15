![logo](QBRC.jpg)
# <name>
## introduction
<name> is a preprocessing tool designed to convert scRNAseq results into a format compatible with analysis tools designed for general high-throughput sequencing data such as DNAseq. The main advantage of scRNAseq is that it allows researchers to sequence individual cells, whereas more traditional techniques sequence the aggregate genetic material from a population of cells. However, this means that scRNAseq results contain a mix of reads from perhaps thousands of diffrent cells. <name> seperates sequecing reads from scRNAseq by their cell of origin and performs a preliminary filter using bowtie2 alignment, producing high-quality inputs that is compatible with existing custom analysis tools.
## Getting started
  
