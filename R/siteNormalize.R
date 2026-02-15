normalize <- function(sitesGRanges, guitarTxdb, txType,overlapIndex,siteLengthIndex)
{ 
  sitesInformation <- data.frame()
  
  sitesPointsPositionNormalize <- list()
  
  sitesPointsPosTx<-list()
  
  sitesPointsPosTx<-end(sitesGRanges[[txType]])
  
  names(sitesPointsPosTx)<-seqnames(sitesGRanges[[txType]])
  #step 1
  startPointMat <- guitarTxdb[[txType]]$startPoint[names(sitesPointsPosTx), ]
  startPointDiffer <- startPointMat - sitesPointsPosTx
  # Vectorized: find last negative column per row using max.col
  neg_mask <- startPointDiffer < 0
  # Set non-negative entries to 0, negative to column index; max.col finds rightmost
  col_idx_mat <- col(startPointDiffer) * neg_mask
  sitesPointsComponet <- max.col(col_idx_mat, ties.method = "last")
  # step 2
  sitesPointsPositionComponet<-startPointDiffer[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)] * -1
  # step 3
  sitesPointsComponetMat<-startPointMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  # step 4
  sitesPointsComponetWidthAvg <- guitarTxdb[[txType]]$componentWidthAverage_pct[sitesPointsComponet]
  # step 5
  sitesPointsComponetStart_pct <- guitarTxdb[[txType]]$componentStartAverage_pct[sitesPointsComponet]
  # step 6
  componentWidthMat <- guitarTxdb[[txType]]$componentWidth[names(sitesPointsPosTx), ]
  sitesPointsComponetWidth <- componentWidthMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  # step 7
  sitesPointsPositionNormalize <- sitesPointsPositionComponet / sitesPointsComponetWidth * sitesPointsComponetWidthAvg + sitesPointsComponetStart_pct
  names(sitesPointsPositionNormalize) <- sitesGRanges[[txType]]$xHits
  #step 8 
  sitesComponet_pct <- guitarTxdb[[txType]]$componentWidthPtc[names(sitesPointsPosTx), ]
  # Vectorized: direct matrix indexing instead of per-element lapply
  sitesPointsComponet_pct <- sitesComponet_pct[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  sitesPointsWeight <-  sitesPointsComponetWidthAvg / (sitesGRanges[[txType]]$pointsOverlapTx^overlapIndex )/ sitesPointsComponet_pct * (sitesGRanges[[txType]]$sitesLength ^ siteLengthIndex)
  names(sitesPointsWeight) <- sitesGRanges[[txType]]$xHits 
  return(list(sitesPointsPositionNormalize,sitesPointsWeight))
}
