#' This function is used to extract peak features significantly higher in treatment
#'
#' @param peaks 
#' @param cutoff 
#' @param pcutoff 
#' @export
#' @return

DiffPeaks<-function(peaks, cutoff, pcutoff,control,treat){
  #'select the samples from control and treatment groups
  control<-grep(control,colnames(peaks))
  treat<-grep(treat,colnames(peaks))
  
  #'fold change between control and treatment
  peaks$fold<-rowMeans(peaks[,treat])/rowMeans(peaks[,control])
  
  #'calculate pvalues
  peaks$pvalue<-rep(1,nrow(peaks))
  for (i in 1:nrow(peaks)){
    #'excluding error for ttest
    if (sd(peaks[i,c(control,treat)])==0){
      next
    }
    
    #'calculating pvalues
    pvalue<-t.test(peaks[i,control],peaks[i,treat])
    peaks$pvalue[i]<-pvalue$p.value
  }
  
  #' selecting the peak features significantly higher in treatment
  index.fold<-which(peaks$fold>cutoff)
  index.sig<-which(peaks$pvalue[index.fold]<pcutoff)
  
  return(peaks[index.fold[index.sig],])
}