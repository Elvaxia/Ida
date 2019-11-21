# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Generate reads per base of genomic range
#' @param SampleInfo a data frame with columns named sample,Bam,Seq_depth,Group
#' @param GR a genomic range of unioned 3'utr frame
#' @param mapq a numeric for mapping quality,default is 255 for uniquelly mapped reads of STAR output,
#' Can be set to 30 when using Bowtie
#' @return A data frame of reads per base of given samples
#' @author Xia Xiao
#' @details function which generate reads per base,default for properly paired ,uniquelly mapped reads
#'and normalized by sequence depth
#' @export
GenerateReadsPerBase <- function(GR, SampleInfo, mapq=255, ...) {
  flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                      isNotPassingQualityControls = FALSE,
                      isUnmappedQuery = FALSE,
                      isPaired = TRUE,
                      isProperPair = TRUE,
                      isDuplicate = FALSE)
  p_param <- PileupParam(distinguish_strands=FALSE,
                         distinguish_nucleotides=FALSE,
                         min_mapq = 1L,
                         max_depth=2e4,
                         min_nucleotide_depth = 0L)
  sbp <- ScanBamParam(which=GR, flag=flag,mapqFilter = mapq)
  bamFile <- SampleInfo$Bam
  tab <- mclapply(FUN = function(x) as.data.frame(pileup(x, scanBamParam=sbp, pileupParam=p_param)),
                  bamFile,
                  mc.cores=5)
  #check if no reads in all sample
  checkreads <- do.call(c,lapply(tab, function(x)nrow(x)))
  if(sum(checkreads) == 0){
    tab.counts <- NULL
  }else{
    tab.tb <- lapply(tab, function(x)x[,2:3])
    tab.counts <- Reduce(function(x,y)merge(x,y,by='pos',all=T),tab.tb)
    tab.counts[is.na(tab.counts)] <- 0
    colnames(tab.counts)[-1] <- SampleInfo$sample
    tab.counts.n <- t(tab.counts[,-1]) * (1e6 /SampleInfo$Seq_depth)
    tab.counts$Total.reads <- rowSums(tab.counts[,-1])
    tab.counts$Normalized.reads <- colSums(tab.counts.n)
  }
  return(tab.counts)
}




#' Changepoint of normalized coverage
#' @param SampleInfo a data frame with columns named sample,Bam,Seq_depth,Group
#' @param GR a genomic range of unioned 3'utr frame
#' @return A data frame of candidate poly(A) sites
#' @author Xia Xiao
#' @details function which compute the change points of reads coverage of 3'utr frame
#' @export
DetectPAs <- function(GR, SampleInfo, cut.off=30, ...){
  cv <- suppressWarnings(GenerateReadsPerBase(GR=GR,
                             SampleInfo=SampleInfo
                             ))
  samples <- SampleInfo$sample
  if(is.null(cv)|mean(cv$Total.reads) <= cut.off){
    GR.res <- as.data.frame(matrix(NA,1,length(samples)*3 +1))
    colnames(GR.res) <-c(paste0(samples,'.L_Coverage'),paste0(samples,'.S_Coverage'),
                         paste0(samples,'.RUD'),'PA')
    GR.res$gene_id <- unique(GR$use.names)
    GR.res$strand <- unique(strand(GR))
    GR.res$chr <- unique(seqnames(GR))
    GR.res$utr_start <- min(start(GR))
    GR.res$utr_end <- max(end(GR))
  }else{
    GR.res <- ChangePoint.RUD(GR=GR,SampleInfo=SampleInfo,cv=cv)
  }
  return(GR.res)
}


#' Changepoint of normalized coverage
#' @param GR a genomic range of unioned 3'utr frame
#' @param SampleInfo a data frame with columns named sample,Bam,Seq_depth,Group
#' @param cv a data frame generate by GenerateReadsPerBase function
#' @return A data frame of candidate poly(A) sites with Cov_L and Cov_S
#' @author Xia Xiao
#' @details function which compute the change points of reads coverage of 3'utr frame
#' @export
ChangePoint.RUD <- function(GR,SampleInfo,cv){
  samples <- SampleInfo$sample
  Str <- unique(as.character(strand(GR)))
  n <- nrow(cv)
  #at most two PA sites,the cpt.meanvar will out put 3 sites,
  #the last one alway at the last element which need to be removed.
  PAs <- suppressWarnings(cpt.meanvar(cv$Normalized.reads,method="BinSeg",
                                      test.stat="Exponential",Q=2,
                                      pen.value = 0.5,penalty='Manual'))
  nPAs <- length(PAs@cpts)-1
  if(nPAs >=1){
    PA <- PAs@cpts[1:nPAs]
    PA.df <- as.data.frame(do.call(rbind,lapply(PA,function(PA){
      if(Str == '-'){
        L.cv <- colMeans(cv[1:PA,samples])
        S.cv <- colMeans(cv[PA:n,samples])
      }else{
        S.cv <- colMeans(cv[1:PA,samples])
        L.cv <- colMeans(cv[PA:n,samples])
      }
      RUD <- L.cv/S.cv
      res <- c(L.cv,S.cv,RUD)
      return(res)
    })))
  }else{
    PA <- 1
    PA.df <- as.data.frame(matrix(NA,1,length(samples)*3))
  }
  colnames(PA.df) <- c(paste0(samples,'.L_Coverage'),paste0(samples,'.S_Coverage'),
                       paste0(samples,'.RUD'))
  PA.df$PA <- cv$pos[PA]
  PA.df$gene_id <- unique(GR$use.names)
  PA.df$strand <- unique(strand(GR))
  PA.df$chr <- unique(seqnames(GR))
  PA.df$utr_start <- min(start(GR))
  PA.df$utr_end <- max(end(GR))
  return(PA.df)
}
