#' Discover bimodal distrubition features
#'
#' Find bimodal distrubition features and
#' divide the samples into 2 groups by k-means clustering.
#'
#' @param x a numeric matrix with feature rows and sample columns, e.g.,
#' splicing score matrix from \code{\link{spliceGenome}} or
#' \code{\link{spliceGene}} function.
#' @param p.value p.value threshold for bimodal distrubition test
#' @param maf minor allele frequency threshold in k-means clustering
#' @param miss missing grouping rate threshold in k-means clustering
#' @param fold fold change threshold between the two groups
#' @param log whether the scores are to be logarithmic. If TRUE, all the scores
#' are log2 tranformed before k-means clustering: x = log2(x+1).
#' @param cores threads to be used. This value is passed
#'  to \code{\link[parallel]{mclapply}} in \bold{parallel} package
#' @return a matrix with feature rows and sample columns.
#' @details The matrix contains 1, 2 and NA, and
#' values of 'x' in group 2 are larger than group 1.
#' @import parallel
#' @export BMfinder
#' @examples
#' data(rice.bg)
#' score<-spliceGene(rice.bg,'MSTRG.183',junc.type='score')
#' score<-round(score,2)
#' as<-BMfinder(score,cores=1) # 4 bimodal distrubition features found
#'
#' ##compare
#' as
#' score[rownames(score)%in%rownames(as),]

BMfinder<-function(x,p.value=0.01,maf=0.05,miss=0.05,
                   fold=2,log=FALSE,cores=detectCores()-1){
    if (!methods::is(x,'matrix') || !is.numeric(x) ) {
        stop('x must be a numeric matrix!')
    }
    qnorm.value<-stats::qnorm(1-p.value)
    ind<-colnames(x)
    if (log) x<-log2(x+1)
    if (ncol(x)<30) {
        warning('sample size < 30, the result should be further tested!')
    }
    cat("---BMfinder:", nrow(x), "features *", ncol(x),"samples\n")
    res<-mclapply(seq_len(nrow(x)),function(i) {
        R <- x[i, ]
        if (all(R==R[1],na.rm=TRUE)) return(rep(NA,length(ind)))
        na.id<-which(is.na(R))
        if (length(na.id)>0) R[na.id]<- mean(R,na.rm=TRUE)
        R<-sort(R);
        Rcluster<-cluster::pam(R,2,cluster.only =TRUE);
        d <- split(R,Rcluster);
        d1<-d[[1]];d2<-d[[2]]
        m1<-mean(d1); m2<-mean(d2);
        s1<-sd(d1);   s2<-sd(d2);
        z1<-(d1-m2)/s2;z2<-(d2-m1)/s1;
        Rcluster[names(d1[abs(z1)<=qnorm.value])]<-NA;
        Rcluster[names(d2[abs(z2)<=qnorm.value])]<-NA;
        Rcluster<-Rcluster[ind]
        if (length(na.id)>0) Rcluster[na.id]<-NA;
        tab<-table(Rcluster)
        if (length(tab)<2 | min(tab)/sum(tab) < maf) return(rep(NA,length(ind)))
        fc<-mean(x[i,Rcluster==2],na.rm=TRUE)/mean(x[i,Rcluster==1],na.rm=TRUE)
        if ( fc >= fold) {
            return(Rcluster)
        }else{
            return(rep(NA,length(ind)))
        }
    }, mc.cores = cores)
    res<-do.call('rbind',res)
    dimnames(res)<-dimnames(x)
    res<-subset(res,rowSums(is.na(res)) <= length(ind)*miss)
    cat ("---Completed:  a total of",nrow(res),
         "bimodal distrubition features found!","\n");
    res;
}
