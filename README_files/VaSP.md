<img src="VaSP_logo_s.jpg" align='right' alt="logo" width="120" 
 style="vertical-align:middle;margin:20px" />
 
# VaSP: Quantification and Visulization of <br>Variations of Splicing in Population
*by [Huihui Yu](https://github.com/yuhuihui2011), [Qian Du](https://github.com/purod) and Chi Zhang*

## Table of contents
- [1. Introduction](#1-introduction)
- [2. Citation](#2-citation)
- [3. Installation](#3-installation)
- [4. Data input](#4-data-input)
- [5. Quick start](#5-quick-start)
- [6. Functions](#6-functions)
	- [6.1 getDepth](#61-getdepth)
	- [6.2 getGeneinfo](#62-getgeneinfo)
	- [6.3 spliceGene](#63-splicegene)
	- [6.4 spliceGenome](#64-splicegenome)
	- [6.5 BMfinder](#65-bmfinder)
	- [6.6 splicePlot](#66-spliceplot)

## 1. Introduction
---------------

**VaSP** is an R package for discovery of genome-wide variable alternative splicing events from short-read RNA-seq data and visualizations of gene splicing information for publication-quality multi-panel figures.

![](VaSP.png)

**Figure 1. Overview of VaSP**. **(A)**. The workflow and functions of [VaSP](https://github.com/yuhuihui2011/vasp). The input is an R data object ballgown (see `?ballgown`) produced by a standard RNA-seq data analysis protocol, including mapping with HISAT, assembling with StringTie, and collecting expression information with R package [Ballgown](https://github.com/alyssafrazee/ballgown). VaSP calculates the Single Splicing Strength (3S) scores for all splicing junctions in the genome (`?spliceGenome`) or in a particular gene (`?spliceGene`), identifies genotype-specific splicing (GSS) events (`?BMfinder`), and displays differential splicing information (`?splicePlot`). The 3S scores can be also used for other analyses, such as differential splicing analysis or splicing QTL identification. **(B)**. VaSP estimates 3S scores based on junction-read counts normalized by gene-level read coverage. In this example, VaSP calculates the splicing scores of four introns in a gene X with two transcript isoforms. Only the fourth intron is a full usage intron excised by both the two isoforms and the other three are alternative donor site (AltD) sites or Intron Retention (IntronR), respectively. **(C)**. Visualization of splicing information in gene MSTRG.183 (LOC_Os01g03070), whole gene without splicing scores. **(D)**. Visualization of differential splicing region of the gene MSTRG.183 with splicing score displaying. In C and D, the y-axes are read depths and the arcs (lines between exons) indicate exon-exon junctions (introns). The dotted arcs indicate no junction-reads spanning the intron (3S = 0) and solid arcs indicate 3S > 0. The transcripts labeled beginning with ‘LOC_Os’ indicate annotated transcripts by reference genome annotation and the ones beginning with “MSTRG” are transcripts assembled by StringTie. ([Yu et al., 2021](#2-citation))

## 2. Citation
---------------

Yu, H., Du, Q., Campbell, M., Yu, B., Walia, H. and Zhang, C. (2021), 
Genome‐wide discovery of natural variation in pre‐mRNA splicing and prioritising
causal alternative splicing to salt stress response in rice. ***New Phytol***.
https://doi.org/10.1111/nph.17189

## 3. Installation
---------------

Start R and run:

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BBiocManager::install("VaSP",build_vignettes=TRUE)
    vignette('VaSP')

## 4. Data input
-------------

Users need to follow the manual of R package Ballgown
(<https://github.com/alyssafrazee/ballgown>) to create a ballgown object
as an input for the VaSP package. See `?ballgown` for detailed
information on creating Ballgown objects. The object can be stored in a
`.RDate` file by `save()` . Here is an example of constructing rice.bg object 
from HISAT2+StringTie output

    library(VaSP)
    ?ballgown
    path<-system.file('extdata', package='VaSP')
    rice.bg<-ballgown(samples = list.dirs(path = path,recursive = F) )

## 5. Quick start
--------------

Calculate 3S (Single Splicing Strength) scores, find GSS
(genotype-specific splicing) events and display the splicing
information.

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

![](splicePlot-1.png)

## 6. Functions
------------

Currently, there are 6 functions in VaSP:  
***getDepth***: Get read depth from a BAM file (in bedgraph format)  
***getGeneinfo***: Get gene informaton from a ballgown object  
***spliceGene***: Calculate 3S scores for one gene  
***spliceGenome***: Calculate genome-wide splicing scores  
***BMfinder***: Discover bimodal distrubition features  
***splicePlot***: Visualization of read coverage, splicing information
and gene information in a gene region

### 6.1 getDepth

Get read depth from a BAM file (in bedgraph format) and return a
data.frame in bedgraph file format which can be used as input for
`plotBedgraph` in the **SuShi** package.

    path <- system.file("extdata", package = "VaSP")
    bam_files <- list.files(path, "*.bam$")
    bam_files
    #> [1] "Sample_027.bam" "Sample_102.bam" "Sample_237.bam"

    depth <- getDepth(file.path(path, bam_files[1]), "Chr1", start = 1171800, end = 1179400)
    head(depth)
    #>   chrom   start    stop value
    #> 1  Chr1 1171799 1171859     0
    #> 2  Chr1 1171859 1171899     1
    #> 3  Chr1 1171899 1171902     2
    #> 4  Chr1 1171902 1171903     4
    #> 5  Chr1 1171903 1171909     5
    #> 6  Chr1 1171909 1171911     6

    library(Sushi)
    #> Loading required package: zoo
    #> 
    #> Attaching package: 'zoo'
    #> The following objects are masked from 'package:base':
    #> 
    #>     as.Date, as.Date.numeric
    #> Loading required package: biomaRt
    par(mar=c(3,5,1,1))
    plotBedgraph(depth, "Chr1", chromstart = 1171800, chromend = 1179400, yaxt = "s")
    mtext("Depth", side = 2, line = 2.5, cex = 1.2, font = 2)
    labelgenome("Chr1", 1171800, 1179400, side = 1, scipen = 20, n = 5, scale = "Kb")

![](plotBedgraph-1.png)

### 6.2 getGeneinfo

Get gene informaton from a ballgown object by genes or by genomic
regions and return a data.frame in bed-like file format that can be used
as input for `plotGenes` in the **SuShi** package

    data(rice.bg)
    unique(geneIDs(rice.bg))
    #>  [1] "MSTRG.177" "MSTRG.178" "MSTRG.179" "MSTRG.180" "MSTRG.181" "MSTRG.182"
    #>  [7] "MSTRG.183" "MSTRG.184" "MSTRG.185" "MSTRG.186"

    gene_id <- c("MSTRG.181", "MSTRG.182", "MSTRG.183")
    geneinfo <- getGeneinfo(genes = gene_id, rice.bg)
    trans <- table(geneinfo$name)  # show how many exons each transcript has
    trans
    #> 
    #> LOC_Os01g03050.1 LOC_Os01g03060.1 LOC_Os01g03060.2 LOC_Os01g03060.3 
    #>                5                2                3                3 
    #> LOC_Os01g03070.1 LOC_Os01g03070.2      MSTRG.181.1      MSTRG.182.2 
    #>               14               14                6                3 
    #>      MSTRG.183.3      MSTRG.183.4      MSTRG.183.5 
    #>               14               14               14


    chrom = geneinfo$chrom[1]
    chromstart = min(geneinfo$start) - 1e3
    chromend = max(geneinfo$stop) + 1e3
    color = rep(SushiColors(2)(length(trans)), trans)

    par(mar=c(3,1,1,1))
    p<-plotGenes(geneinfo, chrom, chromstart, chromend, col = color, bheight = 0.2, 
              bentline = FALSE, plotgenetype = "arrow", labeloffset = 0.5)
    #> [1] "yes"
    labelgenome(chrom, chromstart , chromend, side = 1, n = 5, scale = "Kb")

![](plotGenes-1.png)

### 6.3 spliceGene

Calculate 3S Scores from ballgown object for a given gene. This function
can only calculate one gene. Please use function `spliceGenome` to
obtain genome-wide 3S scores.

    data(rice.bg)
    rice.bg
    #> ballgown instance with 33 transcripts and 6 samples
    head(geneIDs(rice.bg))
    #>           1           2           3           4           5           6 
    #> "MSTRG.177" "MSTRG.177" "MSTRG.177" "MSTRG.178" "MSTRG.178" "MSTRG.179"

    score <- spliceGene(rice.bg, "MSTRG.183", junc.type = "score")
    count <- spliceGene(rice.bg, "MSTRG.183", junc.type = "count")

    ## compare
    tail(score)
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 54   1.073545  1.1034944  0.9799037  1.0581254  1.0581254  1.0581254
    #> 55   1.073545  1.1034944  0.1224880  0.1594436  0.1594436  0.1594436
    #> 56   0.788727  1.3148018  1.1803386  1.3190330  1.3190330  1.3190330
    #> 57   0.000000  0.0000000  0.5790340  0.4058563  0.4058563  0.4058563
    #> 58   0.000000  0.0234786  0.2227054  0.2319179  0.2319179  0.2319179
    #> 59   0.000000  0.0000000  0.0000000  0.1159589  0.1159589  0.1159589
    tail(count)
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 54         49         47         88         73         73         73
    #> 55         49         47         11         11         11         11
    #> 56         36         56        106         91         91         91
    #> 57          0          0         52         28         28         28
    #> 58          0          1         20         16         16         16
    #> 59          0          0          0          8          8          8

    ## get intron structrue
    intron <- structure(rice.bg)$intron
    intron[intron$id %in% rownames(score)]
    #> GRanges object with 16 ranges and 2 metadata columns:
    #>        seqnames          ranges strand |        id transcripts
    #>           <Rle>       <IRanges>  <Rle> | <integer> <character>
    #>    [1]     Chr1 1172688-1173213      - |        43       15:19
    #>    [2]     Chr1 1173416-1173795      - |        44       15:19
    #>    [3]     Chr1 1173877-1173965      - |        45       15:19
    #>    [4]     Chr1 1174170-1174670      - |        46       15:19
    #>    [5]     Chr1 1175175-1176065      - |        48       15:19
    #>    ...      ...             ...    ... .       ...         ...
    #>   [12]     Chr1 1179011-1179226      - |        55       15:16
    #>   [13]     Chr1 1174881-1174952      - |        56       16:19
    #>   [14]     Chr1 1179011-1179134      - |        57          17
    #>   [15]     Chr1 1179011-1179110      - |        58          18
    #>   [16]     Chr1 1179011-1179106      - |        59          19
    #>   -------
    #>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

### 6.4 spliceGenome

Calculate 3S scores from ballgown objects for all genes and return a
list of two elelments: "score' is matrix of intron 3S scores with intron
rows and sample columns and "intron" is a `GRanges` object of intron
structure.

    data(rice.bg)
    rice.bg
    #> ballgown instance with 33 transcripts and 6 samples

    splice <- spliceGenome(rice.bg, gene.select = NA, intron.select = NA)
    #> ---Calculate gene-level read coverage:
    #>   10 genes selected.
    #> ---Extract intron-level read ucount:
    #>   81 introns in 9 genes selected.
    names(splice)
    #> [1] "score"  "intron"

    head(splice$score)
    #>   Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 1  1.4852698  1.6132865   1.934796   1.834662   1.834662   1.834662
    #> 2  1.1580070  1.5164893   1.924559   1.726741   1.726741   1.726741
    #> 3  1.2838773  1.7423494   2.129300   2.309516   2.309516   2.309516
    #> 4  1.8377067  1.3551606   1.904085   2.050505   2.050505   2.050505
    #> 5  0.9817885  0.9679719   1.289864   1.446146   1.446146   1.446146
    #> 6  1.1076589  1.1615663   1.637923   1.661988   1.661988   1.661988
    splice$intron
    #> GRanges object with 81 ranges and 3 metadata columns:
    #>      seqnames          ranges strand |        id   transcripts     gene_id
    #>         <Rle>       <IRanges>  <Rle> | <integer>   <character> <character>
    #>    1     Chr1 1146321-1146479      - |         1           1:2   MSTRG.177
    #>    2     Chr1 1146612-1147484      - |         2           1:2   MSTRG.177
    #>    3     Chr1 1147563-1148021      - |         3           1:3   MSTRG.177
    #>    4     Chr1 1148200-1148268      - |         4           1:3   MSTRG.177
    #>    5     Chr1 1148442-1148530      - |         5           1:3   MSTRG.177
    #>   ..      ...             ...    ... .       ...           ...         ...
    #>   77     Chr1 1187243-1187436      - |        77 c(28, 29, 32)   MSTRG.186
    #>   78     Chr1 1189347-1190773      - |        78         29:30   MSTRG.186
    #>   79     Chr1 1190863-1190978      - |        79         29:30   MSTRG.186
    #>   80     Chr1 1189347-1189914      - |        80         31:33   MSTRG.186
    #>   81     Chr1 1189961-1190037      - |        81            33   MSTRG.186
    #>   -------
    #>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

### 6.5 BMfinder

Find bimodal distrubition features and divide the samples into 2 groups
by k-means clustering and return a matrix with feature rows and sample
columns.

    data(rice.bg)
    score <- spliceGene(rice.bg, "MSTRG.183", junc.type = "score")
    score <- round(score, 2)
    as <- BMfinder(score, cores = 1)  # 4 bimodal distrubition features found
    #> Warning in BMfinder(score, cores = 1): sample size < 30, the result should be
    #> further tested!
    #> ---BMfinder: 16 features * 6 samples
    #> ---Completed:  a total of 4 bimodal distrubition features found!

    ## compare
    as
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 55          2          2          1          1          1          1
    #> 57          1          1          2          2          2          2
    #> 58          1          1          2          2          2          2
    #> 59          1          1          1          2          2          2
    score[rownames(score) %in% rownames(as), ]
    #>    Sample_027 Sample_042 Sample_102 Sample_137 Sample_237 Sample_272
    #> 55       1.07       1.10       0.12       0.16       0.16       0.16
    #> 57       0.00       0.00       0.58       0.41       0.41       0.41
    #> 58       0.00       0.02       0.22       0.23       0.23       0.23
    #> 59       0.00       0.00       0.00       0.12       0.12       0.12

### 6.6 splicePlot

Visualization of read coverage, splicing information and gene
information in a gene region. This function is a wrapper of `getDepth`,
`getGeneinfo`, `spliceGene`, `plotBedgraph` and `plotGenes`.

    data(rice.bg)
    rice.bg
    #> ballgown instance with 33 transcripts and 6 samples
    samples <- paste("Sample", c("027", "102", "237"), sep = "_")
    bam.dir <- system.file("extdata", package = "VaSP")

    ## plot the whole gene region without junction lables
    splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", junc.text = FALSE, bheight = 0.2)
    #> [1] "yes"

![](unnamed-chunk-9-1.png)


    ## plot the alternative splicing region with junction splicing scores
    splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", start = 1179000)
    #> [1] "yes"

![](unnamed-chunk-9-2.png)

If the bam files are provided (`bam.dir` is not NA), the read depth for
each sample is plotted. Otherwise (`bam.dir=NA`), the coserved exons of
the samples are displayed by rectangles (an example is the figure in
**4. Quick start**). And by default (`junc.type = 'score'`,
`junc.text = TRUE`), the junctions (represented by arcs) are labeled
with splicing scores. You can change the argument `junc.text = FALSE` to
unlabel the junctions or change the argument `junc.type = 'count'` to
label with junction read counts.

    splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", junc.type = 'count', start = 1179000)
    #> [1] "yes"

![](unnamed-chunk-10-1.png)

There are other more options to modify the plot, please see the function
`?splicePlot` for details.
