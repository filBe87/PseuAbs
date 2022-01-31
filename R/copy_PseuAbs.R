### =========================================================================
### define copy_PseuAbs
### =========================================================================
#' Copy pseudoabsences from an existing wsl.pseudoabsences object
#'
#' Copy pseudoabsences from existing wsl.pseudoabsences object to a wsl.pseudoabsences
#' object generated for new species and its presence observations.
#'
#' @param env.stack RasterStack/RasterBrick with environmental layers for sampling and
#' extraction. Needs to contain the same layers as in the wsl.pseudoabsences object.
#' @param pres SpatialPoints object with presence observations.
#' @param taxon Character; name of taxon of interest to store in meta
#' information.
#' @param x A wsl.pseudoabsences object
#' @details if the desired pseudoabsence sampling strategy is not species-specific
#' it may be more efficient to copy on the sampled pseudoabsences from another
#' object. Note that this function does not allow for any thinning of presence observations.
#' @return an object of class 'wsl.pseudoabsences'
#' @author Philipp Brun
#' @export
#' @examples
#'### =========================================================================
#'### Data preparation
#'### =========================================================================
#'
#'# Predictors
#'bio=raster::getData('worldclim',var='bio',lon=16, lat=48,res=.5)
#'bio=bio[[c(1,4,12)]]
#'
#'# install.packages("rgbif")
#'library(rgbif)
#'# extract species
#'spn1='Boletus aestivalis'
#'xt=as.vector(extent(bio))
#'baest <- occ_search(scientificName=spn1,
#'                   hasCoordinate=TRUE,
#'                   decimalLongitude=paste0(xt[1],",",xt[3]),
#'                   decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'pbaest=baest$data[,c('decimalLongitude','decimalLatitude')]
#'baest_spp=SpatialPoints(pbaest,proj4string = crs(bio))
#'
#'spn2='Boletus edulis'
#'bedu <- occ_search(scientificName=spn2,
#'                   hasCoordinate=TRUE,
#'                   decimalLongitude=paste0(xt[1],",",xt[3]),
#'                   decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'pbedu=bedu$data[,c('decimalLongitude','decimalLatitude')]
#'bedu_spp=SpatialPoints(pbedu,proj4string = crs(bio))
#'
#'### =========================================================================
#'### Sample pseudoabsences
#'### =========================================================================
#'
#'# Geo.start
#'pseu.abs1=wsl.samplePseuAbs(type="geo.strat",
#'                            n=5000,
#'                            env.stack=bio,
#'                            pres=baest_spp,
#'                            geores_fact=3,
#'                            add.strat=0,
#'                            taxon=spn1)
#'
#'plot(pseu.abs1)
#'
#'pseu.abs2=copy_PseuAbs(env.stack=bio,
#'                       pres=bedu_spp,
#'                       taxon=spn2,
#'                       x=pseu.abs1)
#'
#'plot(pseu.abs2)
#'
#'
copy_PseuAbs<-function(env.stack,
                       pres=SpatialPoints(),
                       taxon=character(),
                       x){

  ### ------------------------
  ### generate wsl.pseudoabsences object and add meta info
  ### ------------------------

  out<-preva.meta(type="pseudoabsence")

  out@meta$type=x@meta$type
  out@meta$taxon=taxon
  out@meta$template_file=x@meta$template_file

  call=match.call()
  out@call<-call

  ### ------------------------
  ### Extract and refine presences
  ### ------------------------

  if(length(pres)>0){

    xt_pres=extract(env.stack,pres)
    sna=apply(xt_pres,1,function(x){
      any(is.na(x))
    })

    if(length(which(sna))>0){
      pres=pres[-which(sna)]
      xt_pres=xt_pres[-which(sna),]
      cat(paste0(length(which(sna))," non-matching presences removed..\n"))
    }
  }

  ### ------------------------
  ### copy pseudoabsences
  ### ------------------------

  wipseu=which(x@pa==0)
  out@pa=c(rep(1,length(pres)),
           rep(0,length(wipseu)))

  out@env_vars=as.data.frame(rbind(xt_pres,x@env_vars[wipseu,]))
  out@xy=rbind(coordinates(pres),x@xy[wipseu,])

  # return
  return(out)

}
