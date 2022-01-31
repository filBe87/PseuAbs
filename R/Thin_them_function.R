### =========================================================================
### thin them function
### =========================================================================
#' Thin a spatial points object by number of points or minimum distance
#'
#' Thin a spatial points object by number of points or minimum distance
#'
#' @param spdf SpatialPoints or SpatialPointsDataFrame object
#' @param lim_dist Minimum tolerated distance between the points
#' @param lim_n Number of points to be kept
#' @details iteratively the closest point pairs are compared and the
#' point with the lower overall distinctiveness is removed
#' @return thinned SpatialPoints or SpatialPointsDataFrame object
#' @author Philipp
#' @keyword Internal
#' @noRd
#'
thin_them<-function(spdf,lim_dist=NA,lim_n=NA){

  if(is.na(lim_dist) & is.na(lim_n)){
    stop("Please specify either minimum distance or number of points!")
  }

  ste=1:nrow(spdf@coords)

  # Create distance matrix from points
  proje=grepl("longlat",spdf@proj4string)
  dsts<-spDists(spdf,longlat = proje)
  colnames(dsts)=rownames(dsts)=ste


  # Define criteria according to specifications
  if(is.na(lim_n)){
    lim_n=0
  } else {
    lim_dist=max(dsts)
  }

  wi<-which(dsts<lim_dist)
  wi=wi[order(dsts[wi])]

  k <- arrayInd(wi, dim(dsts))
  k <- k[-which(k[,1]>=k[,2]),,drop=F]

  if(nrow(k)>0){

    remo=vector()

    while(nrow(k)>0 & length(remo)!=nrow(spdf@coords)-lim_n){

      # Decide which one to kill
      wk1=which(apply(k,1,function(x,y){y%in%x},y=k[1,1]))
      wk2=which(apply(k,1,function(x,y){y%in%x},y=k[1,2]))

      if(length(wk1)==1 & length(wk2)==1){
        chc=sample(1:2,1)
      } else if(length(wk1)==1 & length(wk2)>1){
        chc=1
      } else if(length(wk1)>1 & length(wk2)==1){
        chc=2
      } else {
        chc=ifelse(mean(wk1)>mean(wk2),1,2)
      }

      remo<-append(remo,colnames(dsts)[k[1,chc]])
      k<-k[-which(k[,1]==k[1,chc] | k[,2]==k[1,chc]),,drop=F]

    }

    out=spdf[-which(ste%in%remo),]

  } else {

    out=spdf

  }

  return(out)
}



