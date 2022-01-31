### =========================================================================
### Create environmental strata
### =========================================================================
#' Create environmental strata for environmentally stratified pseudoabsence
#' sampling.
#' @param env.stk raster stack with environmental layers of interest
#' @param rAll should raster be read entirely into memory? This is faster
#' but may blow the RAMs of your machine.
#' @param save_it  Should a stratified sample of c. 400000 points be saved
#' @param strat_dir Directory to save the stratified sample in,

#' Not to be called directly by the user
#'
#' @author Philipp
#' @keyword Internal
#' @noRd
#'


create_envstrat=function(env.stk,
                         rAll=TRUE,
                         save_it=TRUE,
                         strat_dir=NA,
                         poolsiz=5*10^6,
                         sampsiz=400000,
                         type){

  recli=list()
  # loop over environmental layers
  for(i in 1:nlayers(env.stk)){

    # read layer
    if(rAll){
      rsti=readAll(env.stk[[i]])
    } else {
      rsti=raster(env.stk[[i]])
    }

    # take a sample to check if raster has less than 10 levels
    smp=sampleRandom(rsti,1000)

    if(length(unique(smp))<10){

      warning(paste("Less than 10 unique values found for layer",names(env.stk)[i],"no reclassification applied.."))
      values(rsti)=as.numeric(as.factor(values(rsti)))
      recli[[i]]=rsti

    } else {

      # Define five equidistant bins
      strt_i=seq(from=rsti@data@min,
                 to=rsti@data@max,
                 length.out=6)

      # Reclassify raster according to these bins
      rclm=matrix(c(strt_i[1:5],strt_i[2:6],1:5*10^(i-1)),ncol=3)
      recli[[i]]=reclassify(rsti,rclm,include.lowest=TRUE)

    }


    cat(paste0("reclassified layer",i,"...\n"))
  }

  # create one raster layer representing combined bins
  recl=sum(stack(recli))

  # Create ca 400000 stratfied point samples
  if(rAll){
    sptz=as(recl,"SpatialPixelsDataFrame")
    poolsiz=nrow(sptz)
  } else {
    sptz=sampleRandom(recl,poolsiz,sp=TRUE)
  }
  sptz@proj4string=env.stk@crs

  strpts=stratify(sptz,type,sampsiz)

  if(save_it){
    lyn=paste(names(env.stk),collapse="_")
    save(strpts,file=paste0(strat_dir,type,"_",lyn,".RData"))
  }

  return(strpts)
}

