#' Gene Splicing Plot
#'
#' Visualization of read coverage, splicing information and gene information in
#' a gene region. This function is a wrapper of \code{\link{getDepth}},
#' \code{\link{getGeneinfo}}, \code{\link{spliceGene}},
#' \code{\link[Sushi]{plotBedgraph}} and \code{\link[Sushi]{plotGenes}}.
#' @param bg ballgown object. See \code{\link[=ballgown-class]{ballgown}}.
#' @param gene string indicating a gene ID (must be in the 'bg')
#' @param samples names of the samples to be shown (must be in the 'bg' and
#' have bam files in the 'bam.dir')
#' @param bam.dir bam file directory of the samples.
#' If NA, instead of read depth, conserved exons are drawn.
#' @param start start position to be shown. If NA, start position of the gene
#' will be used.
#' @param end stop position to be shown. If NA, end position of the gene will
#' be used.
#' @param labels labels for samples (default: sample names).
#' If it is NA, neigher sample names nor gene names will be labeled
#' @param junc.type type of junction estimates to be shown ('score' for
#' junction score; count' for junction read count)
#' @param junc.text TRUE/FALSE indicating whether junction estimates should
#' be labeled
#' @param trans.select logical expression-like string, indicating transcript
#' rows to select from a matrix of transcript coverages: NA value keeps all
#' transcripts. See \code{\link{spliceGene}}
#' @param junc.select logical expression-like string, indicating junction rows
#' to select from a matrix of junction counts: NA value keeps all junctions.
#' See \code{\link{spliceGene}}
#' @param col a vector of length(samples) specifying colors of read depths.
#' @param transparency value between 0 and 1 indicating the degree of 
#' transparency of read depths.
#' @param scale scale of the \code{\link[Sushi]{labelgenome}} ('bp','Kb','Mb')
#' @param plotgenetype string specifying whether the genes should resemble a
#' 'box' or a 'arrow'. See \code{\link[Sushi]{plotGenes}}.
#' @param ... values to be passed to \code{\link[Sushi]{plotGenes}}.
#' @return see \code{\link[Sushi]{plotGenes}}.
#' @import GenomicRanges
#' @importFrom IRanges IRanges subsetByOverlaps mid
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsInUse seqlengths<-
#' @import Sushi graphics
#' @export splicePlot
#' @examples
#'
#' data(rice.bg)
#' rice.bg
#'
#' samples <- paste('Sample', c('027','102','237'),sep='_')
#' bam.dir <- system.file('extdata',package = 'VaSP')
#'
#' ## plot the whole gene region
#' splicePlot(rice.bg,samples,bam.dir,gene='MSTRG.183',bheight=0.2)
#'
#' ## plot the alternative splicing region
#' splicePlot(rice.bg,samples,bam.dir,gene='MSTRG.183',start=1179000)

splicePlot <- function(bg, gene, samples, bam.dir = NA, start = NA, end = NA, 
    labels = samples, junc.type = c("score", "count"), junc.text = TRUE, 
    trans.select = "rowMaxs(x)>=1", junc.select = "rowMaxs(x)>=5", 
    col = SushiColors(2)(length(samples) + 1)[-1], transparency=0.5,
    scale = "Kb", plotgenetype = "arrow", ...) {
    ### check samples and files
    sample_check <- samples[!samples %in% sampleNames(bg)]
    if (length(sample_check) > 0) {
        stop("no such sample(s) <", paste0(sample_check, collapse = "; "), 
            "> !")
    }
    if (!is.na(bam.dir)) {
        bam_files <- file.path(bam.dir, paste0(samples, ".bam"))
        file_check <- bam_files[!file.exists(bam_files)]
        if (length(file_check) > 0) {
            stop("no such bam file(s) <", paste0(file_check, collapse = "; "), 
                "> !")
        }
    }
    
    ### get gene information and transcript labels
    geneinfo <- getGeneinfo(genes = gene, bg = bg, samples = samples, 
                            trans.select = trans.select)
    trans <- GRanges(geneinfo)
    gr <- range(trans, ignore.strand = FALSE)
    chrom <- as.character(seqnames(gr))
    if (is.na(start)) 
        start <- start(gr)
    if (is.na(end)) 
        end <- end(gr)
    gr <- GRanges(chrom, IRanges::IRanges(start, end))
    trans <- subsetByOverlaps(trans, gr)
    start(trans) <- ifelse(start(trans) < start, start, start(trans))
    end(trans) <- ifelse(end(trans) > end, end, end(trans))
    gap <- gaps(trans)
    gap <- subsetByOverlaps(gap, gr, type = "within")
    trans.start <- min(start(trans))
    trans.end <- max(end(trans))
    geneinfo <- as.data.frame.array(as.data.frame(trans))[, c(1, 2, 3, 
        6, 7, 5), drop = FALSE]
    trans <- GenomicRanges::split(trans, trans$name)
    trans_labels <- names(trans)[order(unlist(width(range(trans))),
                                        decreasing = TRUE)]
    
    ### get junction score/count
    score <- spliceGene(bg = bg, gene = gene, samples = samples, 
                        junc.type = junc.type, trans.select = trans.select, 
                        junc.select = junc.select)
    intron <- ballgown::structure(bg)$intron
    intron <- intron[match(rownames(score), intron$id)]
    intron <- keepSeqlevels(intron, seqlevelsInUse(intron))
    intron <- subsetByOverlaps(intron, gr, type = "within")
    intron <- GenomicRanges::sort(intron)
    score <- score[match(intron$id, rownames(score)), , drop = FALSE]
    
    ### get read depth
    if (is.na(bam.dir)) {
        seqlengths(intron) <- trans.end
        exons <- gaps(intron)
        exons <- exons[strand(exons) == strand(intron)[1]]
        start(exons) <- ifelse(start(exons) < trans.start, trans.start, 
            start(exons))
        depth <- data.frame(chrom = chrom, start = start(exons) - 1, 
                            stop = end(exons), value = 4, 
                            stringsAsFactors = FALSE)
        for (i in seq_len(length(samples))) {
            assign(samples[i], depth)
        }
    } else {
        for (i in seq_len(length(samples))) {
            assign(samples[i], getDepth(bam_files[i], chrom, start, end))
        }
    }
    
    ### plot set layout
    def.par <- par(no.readonly = TRUE)
    sample_layout <- seq_len(length(samples))
    trans_layout <- rep(length(samples) + 1, 1.2 + length(trans_labels) * 
        0.3)
    layout(matrix(c(sample_layout, trans_layout)))
    # set lable position
    if (length(gap) == 0) {
        label_x <- mid(gr)
    } else {
        label_x <- mid(gap[which.max(width(gap))])
    }
    # plot depth
    par(mar = c(0.5, 3, 1, 0.5))
    for (i in seq_len(length(samples))) {
        height <- max(get(samples[i])[, 4], na.rm = TRUE)
        height <- ifelse(height < 5, 5, height)
        plotBedgraph(get(samples[i]), chrom, start - 1, end + 1, 
                    transparency = transparency, color = col[i], xaxt = "n", 
                    yaxt = "n", range = c(-height/2, height))
        if (!is.na(bam.dir)) {
            axis(side = 2, las = 2, at = pretty(c(0, height)))
            grid()
        }
        text(x = label_x, y = height, labels = labels[i], font = 2, xpd = TRUE, 
            cex = 1.2, adj = c(0.5, 0.5))
        if (length(intron) == 0) 
            next
        adj <- order(width(intron))[seq_len(ceiling(length(intron)/4))]
        adj <- ifelse(mean(adj%%2) > 0.5, 0, 1)
        for (j in seq_len(length(intron))) {
            if (is.na(score[j, i])) 
                next
            x1 <- start(intron)[j]
            x2 <- end(intron)[j]
            h <- height/3
            lwd <- 1.5
            arc <- function(x) h * sin(pi * (x - x1)/(x2 - x1) * (-1)^(j + 
                adj))
            curve(arc, from = x1, to = x2, lwd = lwd, add = TRUE, xname = "x", 
                type = "l", col = col[i], lty = ifelse(score[j, i] == 
                0, 2, 1))
            if (junc.text) {
                text((x1 + x2)/2, h * (-1)^(j + adj) * 1.5, round(score[j, 
                i], 2), xpd = TRUE, bg = 1)
            }
        }
    }
    # plot gene structure
    par(mar = c(2.5, 3, 0.5, 0.5))
    plotGenes(geneinfo, chrom, start - 1, end + 1, types = geneinfo$type, 
        bentline = FALSE, labeltext = FALSE, plotgenetype = plotgenetype, 
        ...)
    labelgenome(chrom, start - 1, end + 1, side = 1, scale = scale, 
        chromfont = 1.5, chromcex = 0.9, chromadjust = 0, scalefont = 1.5, 
        scalecex = 1.05, scaleadjust = 1, edgeblankfraction = 0.15)
    grid(ny = 0)
    if (!is.na(labels[1])) {
        text(label_x, y = seq_len(length(trans_labels)) + 0.5, 
            labels = trans_labels, font = 2, xpd = TRUE, cex = 1)
    }
    # reset to previous settings
    par(def.par)
}
