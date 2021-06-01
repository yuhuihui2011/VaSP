<img src="README_files/VaSP_logo_s.jpg" align='right' alt="logo" width="120" 
 style="vertical-align:middle;margin:20px" />

# VaSP: Quantification and Visualization of Variations of Splicing in Population 

<!-- badges: start -->
[ Release ](http://bioconductor.org/packages/release/bioc/html/VaSP.html) ![in Bioc](http://bioconductor.org/shields/years-in-bioc/VaSP.svg) ![platform](http://bioconductor.org/shields/availability/release/VaSP.svg) [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/VaSP.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/VaSP) ![update](http://bioconductor.org/shields/lastcommit/release/bioc/VaSP.svg) 
<!-- badges: end -->

<br>

## Table of contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Data input](#data-input)
- [Quick start](#quick-start)
- [Features](#features)
- [Contributors](#contributors)
- [License](#license)
- [Citation](#citation)


## Introduction

**VaSP** is an R package for discovery of genome-wide variable alternative splicing events from short-read RNA-seq data and visualizations of gene splicing information for publication-quality multi-panel figures.

![](README_files/VaSP.png)

**Figure 1. Overview of VaSP**. **(A)**. The workflow and functions of [VaSP](https://github.com/yuhuihui2011/VaSP). The input is an R data object ballgown (see `?ballgown`) produced by a standard RNA-seq data analysis protocol, including mapping with HISAT, assembling with StringTie, and collecting expression information with R package [Ballgown](https://github.com/alyssafrazee/ballgown). VaSP calculates the Single Splicing Strength (3S) scores for all splicing junctions in the genome (`?spliceGenome`) or in a particular gene (`?spliceGene`), identifies genotype-specific splicing (GSS) events (`?BMfinder`), and displays differential splicing information (`?splicePlot`). The 3S scores can be also used for other analyses, such as differential splicing analysis or splicing QTL identification. **(B)**. VaSP estimates 3S scores based on junction-read counts normalized by gene-level read coverage. In this example, VaSP calculates the splicing scores of four introns in a gene X with two transcript isoforms. Only the fourth intron is a full usage intron excised by both the two isoforms and the other three are alternative donor site (AltD) sites or Intron Retention (IntronR), respectively. **(C)**. Visualization of splicing information in gene MSTRG.183 (LOC_Os01g03070), whole gene without splicing scores. **(D)**. Visualization of differential splicing region of the gene MSTRG.183 with splicing score displaying. In C and D, the y-axes are read depths and the arcs (lines between exons) indicate exon-exon junctions (introns). The dotted arcs indicate no junction-reads spanning the intron (3S = 0) and solid arcs indicate 3S > 0. The transcripts labeled beginning with ‘LOC_Os’ indicate annotated transcripts by reference genome annotation and the ones beginning with “MSTRG” are transcripts assembled by StringTie. ([Yu et al., 2021](#citation))

## Installation

Start R (>= 4.0) and run:

```{r,eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VaSP")
vignette('VaSP')
```

If you use an older version of R (>= 3.5), enter:

```{r,eval=FALSE}
BiocManager::install("yuhuihui2011/VaSP", build_vignettes=TRUE)
vignette('VaSP')
```

## Data input

Users need to follow the manual of R package Ballgown (<https://github.com/alyssafrazee/ballgown>) to creat a ballgown object as an input for the VaSP package. See `?ballgown` for detailed information on creating Ballgown objects. The object can be stored in a `.RDate` file by `save()` . Here is an example of constructing rice.bg object from HISAT2+StringTie output

```{r,eval=FALSE}
library(VaSP)
?ballgown
path<-system.file('extdata', package='VaSP')
rice.bg<-ballgown(samples = list.dirs(path = path,recursive = F) )
```

## Quick start

Calculate 3S (Single Splicing Strength) scores, find GSS (genotype-specific splicing) events and display the splicing information.

-   Calculating 3S scores:

<!-- -->

    library(VaSP)
    #> Loading required package: ballgown
    #> 
    #> Attaching package: 'ballgown'
    #> The following object is masked from 'package:base':
    #> 
    #>     structure
    data(rice.bg)
    ?rice.bg
    rice.bg
    #> ballgown instance with 33 transcripts and 6 samples
    score<-spliceGene(rice.bg, gene="MSTRG.183", junc.type = "score")
    tail(round(score,2),2)
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 58          0       0.02       0.22       0.23       0.23       0.23
    #> 59          0       0.00       0.00       0.12       0.12       0.12

-   Discovering GSS:

<!-- -->

    gss <- BMfinder(score, cores = 1) 
    #> Warning in BMfinder(score, cores = 1): sample size < 30, the result should be
    #> further tested!
    #> ---BMfinder: 16 features * 6 samples
    #> ---Completed:  a total of 4 bimodal distrubition features found!
    gss
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 55          2          2          1          1          1          1
    #> 57          1          1          2          2          2          2
    #> 58          1          1          2          2          2          2
    #> 59          1          1          1          2          2          2

-   Extracing intron information

<!-- -->

    gss_intron<-structure(rice.bg)$intron
    (gss_intron<-gss_intron[gss_intron$id%in%rownames(gss)])
    #> GRanges object with 4 ranges and 2 metadata columns:
    #>       seqnames          ranges strand |        id transcripts
    #>          <Rle>       <IRanges>  <Rle> | <integer> <character>
    #>   [1]     Chr1 1179011-1179226      - |        55       15:16
    #>   [2]     Chr1 1179011-1179134      - |        57          17
    #>   [3]     Chr1 1179011-1179110      - |        58          18
    #>   [4]     Chr1 1179011-1179106      - |        59          19
    #>   -------
    #>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
    range(gss_intron)
    #> GRanges object with 1 range and 0 metadata columns:
    #>       seqnames          ranges strand
    #>          <Rle>       <IRanges>  <Rle>
    #>   [1]     Chr1 1179011-1179226      -
    #>   -------
    #>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

-   Showing the splicing information

<!-- -->

    splicePlot(rice.bg,gene='MSTRG.183',samples = sampleNames(rice.bg)[c(1,3,5)],
               start = 1179000, end = 1179300)
    #> [1] "yes"
![](README_files/splicePlot-1.png)

## Features

Detailed usage examples are available in the [Vignette](README_files/VaSP.md).

## Contributors

* [Huihui Yu](https://github.com/yuhuihui2011)
* [Qian Du](https://github.com/purod)
* Chi Zhang
  
## License

The code is freely available under the GPL (>= 2.0) license

## Citation

Yu, H., Du, Q., Campbell, M., Yu, B., Walia, H. and Zhang, C. (2021), 
Genome‐wide discovery of natural variation in pre‐mRNA splicing and prioritising
causal alternative splicing to salt stress response in rice. ***New Phytol***. 230: 1273-1287. 
https://doi.org/10.1111/nph.17189
<br />
<br />
