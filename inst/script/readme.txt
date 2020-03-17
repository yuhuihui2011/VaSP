The extdata directory contains:

1. Three .bam files and the corresponding index files;
2. Six directories of expression information of 6 samples and 33 transcripts.

Source

The raw RNA-seq data were from the project of variation in transcriptional 
responses to salt stress in rice (SRA Accession: 
SRP106054<https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP106054>)


Processing

The raw RNA-seq data were screened and trimmed using Trimmomatic 
(Bolger et al., 2014) and RNA-seq mapping, transcript assembly, and 
quantification were conducted with HISAT, StringTie, and Ballgown by 
following the method described by Pertea et al. (Pertea et al., 2016). The 3 
.bam files and the corresponding index files were ouput from HISAT2 mapping and 
the 6 directories of expression information were ouput from StringTie. 
The dataset rice.bg was a subset of Ballgown with 33 transcripts and 6 samples.
See ?rice.bg for details.


References

Bolger, A.M., Lohse, M., and Usadel, B. (2014). Trimmomatic: a flexible 
trimmer for Illumina sequence data. Bioinformatics 30, 2114-2120.
 
Pertea, M., Kim, D., Pertea, G.M., Leek, J.T., and Salzberg, S.L. (2016). 
Transcript-level expression analysis of RNA-seq experiments with HISAT, 
StringTie and Ballgown. Nat Protoc 11, 1650-1667.
