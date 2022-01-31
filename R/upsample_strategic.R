### =========================================================================
### upsample_thin function
### =========================================================================
#' Sample a supsample from a large set of points with a minimum distance
#' constraint
#'
#' Thin a spatial points object by number of points or minimum distance
#'
#' @param spdf SpatialPoints or SpatialPointsDataFrame object
#' @param lim_dist Minimum tolerated distance between the points
#' @param n_tot Minimum tolerated distance between the points
#' @details iteratively samples points and rejects them if min thist
#' to all existing points i the sample are not respected
#' @return SpatialPoints or SpatialPointsDataFrame with at least
#' lim_dist between all points
#' @author Philipp
#' @keyword Internal
#' @noRd
#'
upsample_strategic<-function(spdf,lim_dist,n_tot,warnig=TRUE){

  proje=grepl("longlat",spdf@proj4string)

  smps=1:nrow(spdf@coords)

  strfy=spdf[sample(smps,100),]
  strfy=thin_them(strfy,lim_n=10)
  fails=0

  strtcor=paste(strfy@coords[,1],strfy@coords[,2])
  fullcor=paste(spdf@coords[,1],spdf@coords[,2])
  smps=smps[-which(strtcor%in%fullcor)]

  while(fails<50 & length(smps)>0){

    smpi=sample(smps,min(30,length(smps)))
    cand=spdf[smpi,]

    sm.cand=apply(cand@coords,1,function(x,y,z,ld){
      dst=spDistsN1(y,x,longlat = z)
      ssd=any(dst<ld)
      return(ifelse(ssd,NA,min(dst)))
    },y=strfy,z=proje,ld=lim_dist)

    if(any(!is.na(sm.cand))){
      strfy=rbind(strfy,cand[which.min(sm.cand),])
      fails=0

      # identify samples to remove
      selected=which.min(sm.cand)
      too.close=which(is.na(sm.cand))

      smps=smps[-which(smps%in%smpi[c(selected,too.close)])]
    } else {
      smps=smps[-which(smps%in%smpi)]
      fails=fails+1
    }

  }

  if(nrow(strfy@coords)<n_tot){
    if(warnig){
      cat(paste("Only",nrow(strfy@coords),"found with enforced min distance constraint of",
                round(lim_dist,digits=2),"...\n"))
    }

  } else {
    strfy=strfy[sample(1:nrow(strfy@coords),n_tot),]
  }

  return(strfy)
}



