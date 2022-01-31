### =========================================================================
### meta.info function
### =========================================================================
#' Generate meta information for prediction and evaluation
#'
#' Not to be called directly by the user
#' @author Philipp
#' @keyword Internal
#' @noRd
preva.meta=function(env=parent.frame(),type=character()){

  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  m.i=list()
  m.i$author=Sys.info()[["user"]]
  m.i$date=Sys.time()

  # Generate pevaluate object
  if(type=="evaluation"){
    out<-wsl.evaluation()
    m.i$wsl.fit=env$x@meta
    m.i$cutoff=env$crit
  } else if(type=="pseudoabsence") {
    out<-wsl.pseudoabsences()
  } else {
    out<-wsl.prediction()
    m.i$wsl.fit=env$x@meta
  }

  #add Meta info
  out@meta<-m.i

  return(out)

}
