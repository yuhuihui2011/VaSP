## ----setup, include = FALSE---------------------------------------------------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("vasp")

## ----eval=FALSE---------------------------------------------------------------
#  library(vasp)
#  ?ballgown
#  path<-system.file('extdata', package='vasp')
#  rice.bg<-ballgown(samples = list.dirs(path = path,recursive = F) )

## -----------------------------------------------------------------------------
library(vasp)
data(rice.bg)
?rice.bg
rice.bg
score<-spliceGene(rice.bg, gene="MSTRG.183", junc.type = "score")
tail(round(score,2),2)

## -----------------------------------------------------------------------------
gss <- BMfinder(score, cores = 1) 
gss

## -----------------------------------------------------------------------------
gss_intron<-structure(rice.bg)$intron
(gss_intron<-gss_intron[gss_intron$id%in%rownames(gss)])
range(gss_intron)

## ----splicePlot, fig.width=8, fig.height=5------------------------------------
splicePlot(rice.bg,gene='MSTRG.183',samples = sampleNames(rice.bg)[c(1,3,5)],
           start = 1179000, end = 1179300)

## ----plotBedgraph, fig.height=4, fig.width=8----------------------------------
path <- system.file("extdata", package = "vasp")
bam_files <- list.files(path, "*.bam$")
bam_files

depth <- getDepth(file.path(path, bam_files[1]), "Chr1", 
                  start = 1171800, end = 1179400)
head(depth)

library(Sushi)
par(mar=c(3,5,1,1))
plotBedgraph(depth, "Chr1", chromstart = 1171800, 
             chromend = 1179400, yaxt = "s")
mtext("Depth", side = 2, line = 2.5, cex = 1.2, font = 2)
labelgenome("Chr1", 1171800, 1179400, side = 1, 
            scipen = 20, n = 5, scale = "Kb")

## ----plotGenes, fig.height=5,fig.width=8--------------------------------------
data(rice.bg)
unique(geneIDs(rice.bg))

gene_id <- c("MSTRG.181", "MSTRG.182", "MSTRG.183")
geneinfo <- getGeneinfo(genes = gene_id, rice.bg)
trans <- table(geneinfo$name)  # show how many exons each transcript has
trans

library(Sushi)

chrom            = geneinfo$chrom[1]
chromstart       = min(geneinfo$start) - 1e3
chromend         = max(geneinfo$stop) + 1e3
color            = rep(SushiColors(2)(length(trans)), trans)

par(mar=c(3,1,1,1))
p<-plotGenes(geneinfo, chrom, chromstart, chromend, col = color, bheight = 0.2, 
          bentline = FALSE, plotgenetype = "arrow", labeloffset = 0.5)
labelgenome(chrom, chromstart , chromend, side = 1, n = 5, scale = "Kb")

## -----------------------------------------------------------------------------
data(rice.bg)
rice.bg
head(geneIDs(rice.bg))

score <- spliceGene(rice.bg, "MSTRG.183", junc.type = "score")
count <- spliceGene(rice.bg, "MSTRG.183", junc.type = "count")

## compare
tail(score)
tail(count)

## get intron structrue
intron <- structure(rice.bg)$intron
intron[intron$id %in% rownames(score)]

## -----------------------------------------------------------------------------
data(rice.bg)
rice.bg

splice <- spliceGenome(rice.bg, gene.select = NA, intron.select = NA)
names(splice)

head(splice$score)
splice$intron

## -----------------------------------------------------------------------------
data(rice.bg)
score <- spliceGene(rice.bg, "MSTRG.183", junc.type = "score")
score <- round(score, 2)
as <- BMfinder(score, cores = 1)  # 4 bimodal distrubition features found

## compare
as
score[rownames(score) %in% rownames(as), ]

## ---- fig.width=8, fig.height=4.4---------------------------------------------
data(rice.bg)
rice.bg
samples <- paste("Sample", c("027", "102", "237"), sep = "_")
bam.dir <- system.file("extdata", package = "vasp")

## plot the whole gene region without junction lables
splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", 
           junc.text = FALSE, bheight = 0.2)

## plot the alternative splicing region with junction splicing scores
splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", start = 1179000)

## ---- fig.width=8, fig.height=4.4---------------------------------------------
splicePlot(rice.bg, samples, bam.dir, gene = "MSTRG.183", 
           junc.type = 'count', start = 1179000)

## -----------------------------------------------------------------------------
sessionInfo()

