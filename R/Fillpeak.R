#' Fill missing values across samples
#'
#' @param xset 
#' @param ppm 
#' @param RT 
#'
#' @return a complete peak feature matrix
#' 
Fillpeak<-function(xset,ppm,btw,msfiles){
  ppm<-ppm/10^6
  xset.input<-xset

#'peak feature matrix  
  idsave<-matrix(rep(0,length(xset.input@groupidx)*length(msfiles)*4),nrow=length(xset.input@groupidx)*length(msfiles),ncol=4)

#' find peaks across groups  
  k<-1
  for (i in 1:length(xset.input@groupidx)){
    index<-unlist(xset.input@groupidx[i])
    sampleid<-xset.input@peaks[index,11]
#' the average mz and Rt retention time across samples    
    mz.value<-mean(xset.input@peaks[index,1])
    rt.value<-mean(xset.input@peaks[index,4])

#' find the peaks across samples
    for (j in 1:length(unlist(phenoData(xset.input)))){
      index<-which(sampleid==j)
      if (length(index)<1){
        idsave[k,]<-c(i,j,mz.value,rt.value)##save groupidx,sampleid
        k<-k+1
      }
    }
  }
  
#' filling missing values  
  index<-which(idsave[,1]==0)
  idsave<-idsave[-index,]
  
#' no value to fill  
  if(nrow(idsave)==0){
    return(xset.input)
  }
  
#' fill peaks  
  newpeak<-matrix(rep(0,nrow(idsave)*5),nrow=nrow(idsave),ncol=5)
  kk<-1
  msfile<-filepaths(xset.input)
  minmz<-min(xset.input@peaks[,1])
  maxmz<-max(xset.input@peaks[,1])
  minrt<-min(xset.input@peaks[,4])
  maxrt<-max(xset.input@peaks[,4])
  if (length(idsave)>0){
    for (k in 1:length(unlist(phenoData(xset.input)))){
      print(c('filling peakID...',k,'of...',length(unlist(phenoData(xset.input)))))
      index<-which(idsave[,2]==k)
      if (length(index)==0){next}
      xraw<-xcmsRaw(msfiles[k])
      for (n in 1:length(index)){
        mz.value<-idsave[index[n],3]
        mzmin<-max(minmz,mz.value-mz.value*ppm)
        mzmax<-min(maxmz,mz.value+mz.value*ppm)
        rt.value<-idsave[index[n],4]
        rtmin<-max(minrt,rt.value-btw)
        rtmax<-min(maxrt,rt.value+btw)
        
#' extract peak intensities        
        peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        intensity<-max(peak$intensity)
        newpeak[kk,]<-c(mz.value,rt.value,intensity,idsave[index[n],1],k)##mz, rt, intensity,groupidx,sampleid
        kk<-kk+1
      }
    }
    }
  
#' organize all peak features into a new peak matrix
  peakid<-length(xset.input@peaks[,1])
  for (i in 1:nrow(newpeak)){
    groupid<-newpeak[i,4]
    tempvalue<-unlist(xset.input@groupidx[groupid])
    peakid<-peakid+1
    xset.input@groupidx[groupid]<-list(c(tempvalue,peakid))
    }
  peak.combine<-matrix(0,ncol=11,nrow=nrow(newpeak))
  peak.combine[,1]<-newpeak[,1]
  peak.combine[,4]<-newpeak[,2]
  peak.combine[,9]<-newpeak[,3]
  peak.combine[,11]<-newpeak[,5]
  xset.input@peaks<-rbind(xset.input@peaks,peak.combine)
  return(xset.input)
}