# RNAseq_analysis
* Purpose: 

> Useful pipeline parts and functions for RNAseq analysis.

**RNAseq tecnhologies:**

|Type|Description|Expected Features (#Raw reads,corr,#variables in final test)|  
|:---------|:-----------------|:------------------------|
|**Conventional RNAseq**|mRNA sequencing. Classic|13-30 million reads. 80-95% Pearson's correlantion by pair of samples after normalization. Low-expressed tags ~ 3-8k genes|
|**Single-Cell RNAseq**| mRNA seq from a single cell| 2-8 million reads. 90% Pearson's correlation by pairs. ~ #13k low-expressed tags.|
|**small RNAseq**|Tipically microRNAs sequencing|4-15 million reads, depending on soMir analysis. Unknown but many zeros as Single-cell.|

-----------------------------
**ToDo List**

*Quaility Control*:

* [X] *Plot* : Cross-species contamination plot by using the Fast_screen output. **Babraham Institute** <http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/>
* [ ] *Plot* : Improved FastQC reports somehow. **Babraham Institute** <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
* [ ] *Plot* : Snapshot of alignment given a Exon ID or Gene ID (EnsEMBL, UCSC, genenames, etc).

*Data exploration*:
* [ ] *R function* : Principal Component Analysis.

*Differential gene expression*:
* [x] *R function* : **Batch effect removal** from gene expression matrix by Generalized Linear Models.
* [ ] *R function* : **limma-voom** implementation.
* [ ] *R function* : **Gene Set test** statistics.
* [ ] *Plot* : **Gene Set test** barcode plot improvements: worm and Enrichment Score.

--------------------------------

**Dependencies**
* `edgeR` (R/Bioconductor)
* `limma` (R/Bioconductor)
* `ggplot2` (R)
* `GSEABase` (R/Bioconductor)
