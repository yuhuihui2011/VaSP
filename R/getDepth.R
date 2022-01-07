#' Get Read Depth
#'
#' Get read depth from a BAM file (in bedgraph format)
#'
#' @param x path to a BAM file
#' @param chrom chromosome of a region to be searched
#' @param start start position
#' @param end end position
#' @return  a data.frame in bedgraph file format
#' @importFrom IRanges IRanges
#' @export getDepth
#' @examples
#' path <- system.file('extdata',package='VaSP')
#' bam_files<-list.files(path,'bam$')
#' bam_files
#'
#' depth<-getDepth(file.path(path, bam_files[1]), 'Chr1',
#'                 start=1171800, end=1179400)
#' head(depth)
#'
#' # library(Sushi)
#' # plotBedgraph(depth,'Chr1',chromstart=1171800, chromend=1179400,yaxt='s')
#' # mtext('Depth',side=2,line=2.5,cex=1.2,font=2)
#' # labelgenome('Chr1',1171800,1179400,side=1,scipen=20,n=5,scale='Kb')

getDepth <- function(x, chrom, start, end) {
    gr <- GRanges(chrom, IRanges(start, end))
    depth <- GenomicAlignments::coverage(x = x, param = Rsamtools::ScanBamParam(
        which = gr))
    depth <- depth[names(depth) == seqnames(gr)][[1]]
    depth <- data.frame(chrom = chrom, start = start(depth) - 1, 
                        stop = end(depth), value = S4Vectors::runValue(depth), 
                        stringsAsFactors = FALSE)
    depth <- depth[depth$stop >= start & depth$start + 1 <= end, , drop = FALSE]
    depth$start[1] <- start - 1
    depth$stop[nrow(depth)] <- end
    rownames(depth) <- seq_len(nrow(depth))
    return(depth)
}
