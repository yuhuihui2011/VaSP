#' Calculate Splicing Scores for One Gene
#'
#' Calculate splicing Scores from ballgown object for a given gene.
#' This function can only calculate one gene. Please use function
#' \code{\link{spliceGenome}} to obtain genome-wide splicing scores.
#'
#' @param bg ballgown object
#' @param gene a character string specifying gene id
#' @param samples names of samples
#' @param junc.type type of junction estimate ('score' for junction score;
#' 'count' for junction read count)
#' @param trans.select logical expression-like string, indicating transcript
#' rows to select from a matrix of transcript coverages: NA value keeps all
#' transcripts. e.g. use trans.select='rowMaxs(x)>=1' to filter the transcrpts
#' with the maximium coverage among all the samples less than 1.
#' @param junc.select logical expression-like string, indicating junction rows
#' to select from a matrix of junction counts: NA value keeps all junctions.
#' e.g. use junc.select='rowMaxs(x)>=5' to filter the junctions with the
#' maximium read count among all the samples less than 5.
#' @return  a matrix of junction scores with intron rows and sample columns.
#' @seealso \code{\link{spliceGenome}}, which calculates splicing scores in
#' whole genome.
#' @details score = junction count/gene-level per base read coverage.
#' Row functions for matrices are useful to select transcripts and junctions.
#' See \code{\link{matrixStats}} package.
#' @import matrixStats
#' @export spliceGene
#' @examples
#' data(rice.bg)
#' rice.bg
#' head(geneIDs(rice.bg))
#'
#' score<-spliceGene(rice.bg,'MSTRG.183',junc.type='score')
#' count<-spliceGene(rice.bg,'MSTRG.183',junc.type='count')
#'
#' ## compare
#' tail(score)
#' tail(count)
#'
#' ## get intron structrue
#' intron<-structure(rice.bg)$intron
#' intron[intron$id%in%rownames(score)]
#' 
#' @references 
#' Yu, H., Du, Q., Campbell, M., Yu, B., Walia, H. and Zhang, C. (2021), 
#' Genome‐Wide Discovery of Natural Variation in Pre‐mRNA Splicing and 
#' Prioritizing Causal Alternative Splicing to Salt Stress Response in Rice. 
#' New Phytol. <https://doi.org/10.1111/nph.17189>

spliceGene <- function(bg, gene, samples = sampleNames(bg), 
                        junc.type = c("score", "count"), 
                        trans.select = "rowMaxs(x)>=1", 
                        junc.select = "rowMaxs(x)>=5") {
    index <- which(geneIDs(bg) == gene)
    if (length(index) == 0) 
        stop("No such gene ", gene)
    cov <- subset(texpr(bg, "cov"), geneIDs(bg) == gene, paste0("cov.", 
        samples))
    if (!is.na(trans.select)) {
        cov <- subset(cov, eval(parse(text = trans.select), list(x = cov)))
    }
    if (nrow(cov) == 0) {
        warning("No transcript in the gene retained!\n")
        return()
    }
    i2t <- indexes(bg)$i2t
    count <- iexpr(bg, "ucount")
    count <- subset(count, rownames(count) %in% i2t[i2t$t_id %in% index, 
        ]$i_id, paste0("ucount.", samples))
    if (!is.na(junc.select) & nrow(count) > 0) {
        count <- subset(count, eval(parse(text = junc.select), list(x = count)))
    }
    if (nrow(count) == 0) {
        warning("No junction retained!\n")
        return()
    }
    colnames(count) <- samples
    if (junc.type[1] == "count") {
        return(count)
    } else if (junc.type[1] == "score") {
        g_cov <- matrix(colSums(cov), nrow = nrow(count), ncol = ncol(count), 
            byrow = TRUE)
        score <- ifelse(g_cov == 0, NA, count/g_cov)
        dimnames(score) <- dimnames(count)
        return(score)
    } else {
        stop("juction must be \"score\" or \"count\"!")
    }
}
