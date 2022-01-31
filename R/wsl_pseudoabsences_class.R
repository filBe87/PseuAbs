### =========================================================================
### Set Class
### =========================================================================
#' An S4 class to store pseudoabsence data
#'
#' Information on coordinates, presence/pseudoabsence, and environmental
#' conditions is organized in seperate slots. Moreover, the slot meta
#' keeps track of meta information.
#'
#' @slot meta a list with meta information
#' @slot pa a vector with 1 encoding presences and 0 encoding pseudoabsences
#' @slot env_vars a data.frame containing extractions of the env.stack at the points xy
#' @slot xy two-column matrix with coordinates of the presence points
#' and the sampled pseudoabsence points. The projection is the same as the one of the env.stack supplied
#' @slot call the call made to the wsl.samplePseuAbs function
#' @author Philipp Brun
#' @export
wsl.pseudoabsences<-setClass("wsl.pseudoabsences",slots=c(meta="list", # Meta information
                                                          pa="numeric", # store presence/pseudoabsence information
                                                          env_vars="data.frame", # store extracted env variables
                                                          xy="matrix", # store coordiantes
                                                          call="call")) # conserve function call
