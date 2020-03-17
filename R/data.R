#' Rice Ballgown Object
#'
#' Small ballgown object created with a subset of rice RNAseq data,
#' for demonstration purposes
#'
#' @format a ballgown object with 33 transcripts and 6 samples
#' @source The raw RNA-seq data were from the project of variation in 
#' transcriptional responses to salt stress in rice (SRA Accession: 
#' \href{https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP106054}{SRP106054})
#' @details The raw RNA-seq data were screened and trimmed using Trimmomatic 
#' (Bolger et al., 2014) and RNA-seq mapping, transcript assembly, and 
#' quantification were conducted with HISAT, StringTie, and Ballgown by 
#' following the method described by Pertea et al. (Pertea et al., 2016). 
#' The rice.bg is a subset ballgown object with 33 transcripts and 6 samples.
#' #> library(vasp)
#' #> ?ballgown
#' #> path<-system.file('extdata', package='vasp')
#' #> rice.bg<-ballgown(samples = list.dirs(path = path,recursive = F))
#' @references 
#' Bolger, A.M., Lohse, M., and Usadel, B. (2014). Trimmomatic: a flexible 
#' trimmer for Illumina sequence data. Bioinformatics 30, 2114-2120.
#' 
#' Pertea, M., Kim, D., Pertea, G.M., Leek, J.T., and Salzberg, S.L. (2016). 
#' Transcript-level expression analysis of RNA-seq experiments with HISAT, 
#' StringTie and Ballgown. Nat Protoc 11, 1650-1667.
#' @examples
#' data(rice.bg)
#' rice.bg
#' # ballgown instance with 33 transcripts and 6 samples
"rice.bg"
