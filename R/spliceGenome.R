#' Calculate Genome-wide Splicing Scores
#'
#' Calculate splicing scores from ballgown objects for all genes.
#'
#' @param bg ballgown object
#' @param gene.select logical expression-like string, indicating genes to select
#' from a matrix of gene-level coverages: NA value keeps all genes.
#' e.g. gene.select = 'rowQuantiles(x,probs = 0.05)>=1' keeps the genes with the
#' read coverage greater than or equal to 1 in at least 95% of the samples
#' (0.05 quantile). Used to filter low expressed genes.
#' @param intron.select logical expression-like string, indicating introns to
#' select from a matrix of junction counts: NA value keeps all introns.
#' e.g. intron.select = 'rowQuantiles(x,probs = 0.95)>=5' keeps the introns
#' with the read count greater than or euqal to 5 in at least 5% of the samples
#' (0.95 quantile). Used to filter introns with very few junction reads
#' supporting.
#' @return  a list of two elelments: 'score' is matrix of intron splicing scores
#' with intron rows and sample columns and 'intron' is a
#' \code{\link[=GenomicRanges-class]{GRanges}} object of intron structure.
#' See \code{\link[ballgown]{structure}} in \bold{ballgown} package
#' @details score = junction count/gene-level per base read coverage.
#' Row functions for matrices in \code{\link{matrixStats}} package are useful
#' to select genes and introns.
#' @export spliceGenome
#' @seealso \code{\link{spliceGene}}, which calculates splicing scores in one
#' gene.
#' @examples
#' data(rice.bg)
#' rice.bg
#'
#' splice<-spliceGenome(rice.bg,gene.select=NA,intron.select=NA)
#' names(splice)
#'
#' head(splice$score)
#' splice$intron
#' 
#' @references 
#' Yu, H., Du, Q., Campbell, M., Yu, B., Walia, H. and Zhang, C. (2021), 
#' Genome‐wide discovery of natural variation in pre‐mRNA splicing and 
#' prioritising causal alternative splicing to salt stress response in rice. 
#' New Phytol. <https://doi.org/10.1111/nph.17189>

spliceGenome <- function(bg, gene.select = "rowQuantiles(x,probs = 0.05)>=1", 
    intron.select = "rowQuantiles(x,probs = 0.95)>=5") {
    options(stringsAsFactors = FALSE)
    cat("---Calculate gene-level read coverage:\n")
    t_cov <- texpr(bg, "cov")
    t_cov <- split.data.frame(t_cov, geneIDs(bg))
    g_cov <- lapply(t_cov, colSums)
    g_cov <- do.call("rbind", g_cov)
    if (!is.na(gene.select)) {
        g_cov <- subset(g_cov, subset = eval(parse(text = gene.select), 
            list(x = g_cov)))
    }
    cat("\t", nrow(g_cov), "genes selected.\n")
    
    cat("---Extract intron-level read ucount:\n")
    intron<-ballgown::structure(bg)$intron
    intron$gene_id<-geneIDs(bg)[vapply(intron$transcripts,function(x) {
        as.numeric(eval(parse(text=x))[1])
    }, 1)]
    intron<-intron[intron$gene_id%in%rownames(g_cov)]    
    ie <- iexpr(bg, "ucount")
    ie <- subset(ie,rownames(ie)%in%intron$id)
    if (!is.na(intron.select)) {
        ie <- subset(ie, eval(parse(text = intron.select), list(x = ie)))
    }
    intron<-intron[match(rownames(ie),intron$id)]
    g_cov <- g_cov[match(intron$gene_id, rownames(g_cov)), , drop = FALSE]
    cat("\t", nrow(g_cov), "introns in", length(unique(intron$gene_id)), 
        "genes selected.\n")
    score <- ifelse(g_cov == 0, NA, ie/g_cov)
    rownames(score) <- intron$id
    colnames(score) <- sampleNames(bg)
    list(score = score, intron = intron)
}
