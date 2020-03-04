#' Calculate Genome-wide Splicing Scores
#'
#' Calculate splicing scores from ballgown objects for all genes.
#'
#' @param bg ballgown object
#' @param gene.select logical expression-like string, indicating genes to select
#' from a matrix of gene-level coverages: NA value keeps all genes.
#' e.g. gene.select = "rowQuantiles(x,probs = 0.05)>=1" keeps the genes with the
#' read coverage greater than or equal to 1 in at least 95% of the samples
#' (0.05 quantile). Used to filter low expressed genes.
#' @param intron.select logical expression-like string, indicating introns to
#' select from a matrix of junction counts: NA value keeps all introns.
#' e.g. intron.select = "rowQuantiles(x,probs = 0.95)>=5" keeps the introns
#' with the read count greater than or euqal to 5 in at least 5% of the samples
#' (0.95 quantile). Used to filter introns with very few junction reads
#' supporting.
#' @return  a list of two elelments: "score' is matrix of intron splicing scores
#' with intron rows and sample columns and "intron" is a
#' \code{\link[=GenomicRanges-class]{GRanges}} object of intron structure.
#' See \code{\link[ballgown]{structure}} in \bold{ballgown} package
#' @details score = junction count/gene-level per base read coverage.
#' Row functions for matrices in \code{\link{matrixStats}} package are useful
#' to select genes and introns.
#' @import GenomicRanges
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

spliceGenome<-function(bg, gene.select = "rowQuantiles(x,probs = 0.05)>=1",
                       intron.select = "rowQuantiles(x,probs = 0.95)>=5"){
    options(stringsAsFactors = FALSE);
    cat("---Calculate gene-level read coverage:\n")
    t_cov<-texpr(bg,'cov');
    t_cov<-split.data.frame(t_cov,geneIDs(bg))
    g_cov<-lapply(t_cov,colSums);
    g_cov<-do.call('rbind',g_cov)
    if(!is.na(gene.select)){
        g_cov<-subset(g_cov, subset=eval(parse(text=gene.select),list(x=g_cov)))
    }
    cat('\t',nrow(g_cov),'genes selected.\n')

    cat("---Extract intron-level read ucount:\n")
    ie<-iexpr(bg,'ucount');
    index <- which(geneIDs(bg) %in% rownames(g_cov))
    i2t <- indexes(bg)$i2t
    i_id  <- sort(unique(i2t$i_id[i2t$t_id %in% index]))
    ie<-ie[i_id,,drop=FALSE]
    if (!is.na(intron.select)){
        ie<-subset(ie, subset=eval(parse(text=intron.select),list(x=ie)))
    }
    i_id <-  as.numeric(rownames(ie))
    g_id <-  geneIDs(bg)[i2t$t_id[match(i_id,i2t$i_id)]]
    g_cov<- g_cov[match(g_id,rownames(g_cov)),,drop=FALSE]
    cat('\t',nrow(g_cov),'introns in',length(unique(g_id)) ,'genes selected.\n')
    score<-ifelse(g_cov==0,NA,ie/g_cov)
    rownames(score)<-i_id
    colnames(score)<-sampleNames(bg)
    intron<-ballgown::structure(bg)$intron[i_id]
    names(intron)<-i_id
    intron$gene_id<-g_id
    list(score=score,intron=intron)
}
