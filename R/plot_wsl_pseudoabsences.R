### =========================================================================
### define summary function for wsl.evaluation objects
### =========================================================================
#' A simple plotting function for objects of class 'wsl.pseudoabsences'
#'
#' This function allos quickly visualizing the different thinning and sampling
#' strategies determine the spatial patterns of presences and pseudoabsences.
#'
#' @param object An object of class wsl.pseudoabsences
#' @return Presences (red crosses) and pseudoabsences (grey transparent points) are plotted in geographic space. If
#' the wsl.pseudoabsences object contains information on a template raster file (see documentation of the
#' wsl.samplePseuAbs function), this template will be plotted in the background.
#' @export
#'
plot.wsl.pseudoabsences=function(object){

  xy_pres=object@xy[which(object@pa==1),]

  xy_abs=object@xy[which(object@pa==0),]
  if(nrow(xy_abs)>10000){
    xy_abs=xy_abs[sample(1:nrow(xy_abs),10000),]
  }

  if(!is.na(object@meta$template_file)){
    rst=raster(object@meta$template_file)
    plot(rst,col=c("#f0f0f0","#99d8c9"),
         main=object@meta$type,legend=F)

  } else {
    plot(xy_abs,
         pch=16,
         xlab="",
         ylab="",
         main=object@meta$type,
         type="n")
  }

  points(xy_abs,pch=16,cex=.2,col="#00000020")
  points(xy_pres,pch=3,cex=.5,col="darkred")

  legend("bottomleft",
         pch=c(16,3),
         col=c("#00000020","darkred"),
         c(paste0("absences (n=",length(which(object@pa==0)),")"),
           paste0("presences (n=",length(which(object@pa==1)),")")))

}
