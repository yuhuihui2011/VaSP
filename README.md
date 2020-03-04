# vasp
Quantification and Visulization of Variations of Splicing in Population
## Introduction

**VaSP** is an R package for discovery of genome-wide variable splicing events from short-read RNA-seq data. Based on R package [Ballgown](https://github.com/alyssafrazee/ballgown), VaSP calculates Single Splicing Strength (3S) score of any intron by the junction count normalized by the gene-level average read coverage (score=junction count/gene-level average read coverage). The 3S scores can be used for further analysis, such as differential splicing analysis between two sample groups and sQTL (splicing Quantitative Trait Locus) identification in a large population. The VaSP package provides a function to find large-effect differential splicing events without the need of genotypic information in an inbred plant population, so called genotype-specific splicing (GSS). Integrated with functions from R package [Sushi](https://github.com/dphansti/Sushi), VaSP package also provides function to visualization of gene splicing information for publication-quality multi-panel figures.

## Installation

Start R and run:
```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("vasp")
```

## Data input

Users need to follow the manual of R package Ballgown (<https://github.com/alyssafrazee/ballgown>) to creat a ballgown object as an input for the VaSP package. See `?ballgown` for detailed information on creating Ballgown objects. The object can be stored in a `.RDate` file by `save()` .

```{r,eval=FALSE}
library(vasp)
?ballgown
path<-system.file('extdata', package='vasp')
rice.bg<-ballgown(samples = list.dirs(path = path,recursive = F) )
```

## Quick start

Calculate 3S (Single Splicing Strength) scores, find GSS (genotype-specific splicing) events and display the splicing information.

-   Calculating 3S scores:

<!-- -->

    library(vasp)
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
![](https://github.com/yuhuihui2011/vasp/blob/master/figure/splicePlot-1.png)

## Features

Detailed usage examples are available in the [Vignette](https://github.com/yuhuihui2011/vasp/blob/master/vignettes/vasp.md).

## Contributors

* [Huihui Yu](https://github.com/yuhuihui2011)
* [Qian Du](https://github.com/purod)
* Chi Zhang
  
## License

The code is freely available under the GPL (>= 2.0) license
