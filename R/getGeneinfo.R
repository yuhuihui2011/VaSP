#' Get Gene Informaton from a ballgown object
#'
#' Get gene informaton from a ballgown object by genes or by genomic regions
#'
#' @param genes a character vector specifying gene IDs in 'bg'. Any values other
#' than NA override genomic region (chrom, start, stop)
#' @param bg ballgown object
#' @param chrom chromosome of a region
#' @param start start postion
#' @param end stop postion
#' @param samples names of samples. The transcrpts in these samples are
#' subjected to 'trans.select'
#' @param trans.select logical expression-like string, indicating transcript
#' rows to select from a matrix of transcript coverages: NA value keeps all 
#' transcripts.
#' @return  a data.frame in bed-like file format that can be used as input for
#' \code{\link[Sushi]{plotGenes}} in the \bold{SuShi} package
#' @import ballgown
#' @export getGeneinfo
#' @seealso \code{\link{splicePlot}}; \code{\link[Sushi]{plotGenes}} in
#' \bold{Sushi} package

#' @examples
#' data(rice.bg)
#' unique(geneIDs(rice.bg))
#'
#' gene_id <- c('MSTRG.181', 'MSTRG.182', 'MSTRG.183')
#' geneinfo <- getGeneinfo(genes=gene_id,rice.bg)
#' trans <- table(geneinfo$name) # show how many exons each transcript has
#' trans
#'
#' library(Sushi)
#' chrom = geneinfo$chrom[1]
#' chromstart = min(geneinfo$start) - 1e3
#' chromend = max(geneinfo$stop) + 1e3
#' color = rep(SushiColors(2)(length(trans)), trans)
#'
#' par(mar=c(3,1,1,1))
#' plotGenes(geneinfo, chrom, chromstart, chromend, col = color, bheight = 0.2,
#'            bentline = FALSE, plotgenetype = 'arrow', labeloffset = 0.5)
#' labelgenome(chrom, chromstart , chromend, side = 1, n = 5, scale = 'Kb')

getGeneinfo <- function(genes = NA, bg, chrom, start, end, 
    samples = sampleNames(bg), trans.select = NA) {
    if (is.na(genes[1])) {
        gr <- GRanges(chrom, IRanges(start, end))
        trans <- subsetByOverlaps(structure(bg)$trans, gr)
        names(trans) <- transcriptNames(bg)[as.numeric(names(trans))]
    } else {
        index <- which(geneIDs(bg) %in% genes)
        if (length(index) == 0) 
            stop("none of the genes exist!")
        trans <- structure(bg)$trans[index]
        names(trans) <- transcriptNames(bg)[index]
    }
    if (!is.na(trans.select)) {
        cov <- subset(texpr(bg, "cov"), transcriptNames(bg) %in% names(trans), 
            paste0("cov.", samples))
        trans <- trans[eval(parse(text = trans.select), list(x = cov))]
    }
    geneinfo <- unlist(trans)
    geneinfo$name <- names(geneinfo)
    names(geneinfo) <- NULL
    geneinfo <- as.data.frame.array(as.data.frame(geneinfo))
    geneinfo$score <- "."
    geneinfo <- geneinfo[, c(1, 2, 3, 8:9, 5)]
    names(geneinfo) <- c("chrom", "start", "stop", "name", "score", "strand")
    geneinfo
}
