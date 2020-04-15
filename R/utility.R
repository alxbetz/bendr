#' Calculate a vector of x-values for plotting
#'
#' @param m.formula
#' model formula
#' @param fdata
#' dataframe with concentration values
#'
#' @return
#' @export
#'
#' @examples
calc_xrange = function(m.formula,fdata) {
  concName = intersect(all.vars(m.formula),colnames(fdata))[2]
  concV = dplyr::pull(fdata,concName)
  loc = min(concV)
  hic = max(concV)
  if (loc>0){
    x.start = 0
  } else {
    x.start = 2 * loc
  }
  x.values = seq(x.start,hic*1.4,length.out=1E3)
}


#' calculate ECx
#'
#' @param x desired effective concentration level (ECx)
#' @param slope slope of the hill curve
#' @param ec50 Effective concentration at 50 percent
#'
#' @return ECx value
#' @export
#'
#' @examples
#' x = 10
#' slope = 2
#' ec50 = 12
#' calc_ecx(x,slope,ec50)
calc_ecx = function(x,slope,ec50) {
  ecx = ((x/(100-x))^(1/slope)) * ec50
  names(ecx) = paste0("ec",as.character(x))
  ecx
}










