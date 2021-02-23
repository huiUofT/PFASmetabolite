#' This function is used to extract peak features
#'
#' @param Intensitycut 
#' @param ppm 
#' @return peak features
#' @export
#' @import xcms
#' @description This function takes raw data and then extract all peaks to dataframe

Peakextract <- function(Intensitycut,ppm) {
#' the path to save results
  path<-getwd()
  
#' the path to save raw data
  path.data<-paste0(path,"/data")
  
#' ----------------------------------------------
#' detect peaks from mass spec raw files
#' ---------------------------------------------
  setwd(path.data)
  msfiles<-list.files()
  xset.raw<-xcmsSet(msfiles,method='centWave',ppm=ppm,peakwidth=c(10,30),snthresh=10,nSlaves=1)##peak width, the min and max range of chromatographic peaks in seconds
  xtest<-xcmsRaw(msfiles[1])
  mzrange<-xtest@mzrange

#' delete peak features with extreme m/z values
  xset<-xset.raw
  index<-which(xset.raw@peaks[,1]<mzrange[1]+1)
  if (length(index)>0){
    xset@peaks<-xset.raw@peaks[-index,]}
  index<-which(xset@peaks[,1]>mzrange[2]-1)
  if (length(index)>0){
    xset@peaks<-xset@peaks[-index,]}
  peaklist<-xset@peaks                                                                    
  len<-length(xset@peaks[,1])
  
#' group peaks across samples
  xset1<-group(xset,bw=60,minsamp=1,minfrac=1/length(msfiles),mzwid=0.001)
  
#' filling missing values
  xset2<-Fillpeak(xset1,10,20,msfiles)

#' peak ID for each group
  test<-unlist(xset2@groupidx)
  len<-length(xset2@groupidx)
  len2<-length(msfiles)
  
#' creating peak matrix
  Allpeak<-array(rep(0,len*(len2+2)),dim=c(len,(len2+2)))
  for (i in 1:len){
    print(c('PeakID...',i,'of...',len))
    temp<-unlist(xset2@groupidx[i])
    len3<-length(temp)
    for (j in 1:len3){
      index1<-xset2@peaks[temp[j],11]
#' m/z vlaues      
      Allpeak[i,1]<-xset2@peaks[temp[j],1]
      
#' retention time in minutes
      Allpeak[i,2]<-xset2@peaks[temp[j],4]/60

#‘ peak intensity
      Allpeak[i,index1+2]<-max(xset2@peaks[temp[j],9],Allpeak[i,index1+2])
    }}
  
#' replacing 0 values
  Allpeak[which(Allpeak==0)]<-100

#' extracting peak intensities across all samples
  colnames(Allpeak)<-c('mz','rt',msfiles)
  Allpeak<-data.frame(Allpeak)
  Allpeak$SampleID<-rep(1,nrow(Allpeak))
  
#' extracting peak intensities
  index.save<-NULL
  for (i in 1:nrow(Allpeak)){
    index<-which.max(Allpeak[i,3:ncol(Allpeak)])
    if (max(Allpeak[i,3:ncol(Allpeak)])>Intensitycut){
      index.save<-c(index.save,i)
      }
    Allpeak$SampleID[i]<-index[1]
  }
  if(length(index.save)>0){
    Allpeak<-Allpeak[index.save,]
  }

#' setup the working folder
  setwd(path)
  
  return(Allpeak)
}
#‘ devtools::document()