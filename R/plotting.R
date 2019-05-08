require(tidyr)
require(ggplot2)
require(dplyr)


#' Plot a single dose response curve
#'
#' @param fitObject
#' drfit result object
#'
#' @return
#' ggplot object
#' @export
#'
#' @examples
plot_single_drc = function(fitObject) {
  p = ggplot(fitObject$plot.data, aes(log.concentration, curve.predict)) +
    geom_point(data=fitObject$data,aes(x=logconc,y=effect)) +
    geom_line(aes(y = curve.predict), colour = "red") +
    geom_line(aes(y = curve.predict + ci.values), colour = "red", linetype = 2) +
    geom_line(aes(y = curve.predict - ci.values), colour = "red", linetype = 2)
  p
}

#' Plot multiple dose response curves
#'
#' @param fit_df
#' data frame containing the model fits in the column m.fit
#' @param plot.layout
#' determines whether all DRC are plotted in one single plot or a plot grid.
#' one of c('single','multi')
#'
#' @return
#' ggplot object
#' @export
#'
#' @examples
plot_multi_drc = function(fit_df,plot.layout = "single") {
  plotdata = fit_df %>%
    dplyr::pull(m.fit) %>%
    map(function(x) x$plot.data)
  names(plotdata) = fit_df %>% dplyr::pull(sampleID)
  plotLong = bind_rows(plotdata,.id = "sampleID")

  rawdata = fit_df %>% dplyr::select(sampleID,data) %>% unnest()

  if(plot.layout == 'single') {
    p = ggplot(plotLong,aes(x=log.concentration,color=sampleID)) +
      geom_line(aes(y=curve.predict)) +
      geom_line(aes(y=curve.predict + ci.values), linetype = 2) +
      geom_line(aes(y=curve.predict - ci.values), linetype = 2) +
      geom_point(data=rawdata,aes(x=logconc,y=effect,color=sampleID))
  } else if(plot.layout == 'multi') {
    p = ggplot(plotLong,aes(x=log.concentration)) +
      geom_line(aes(y=curve.predict)) +
      geom_line(aes(y=curve.predict + ci.values), linetype = 2) +
      geom_line(aes(y=curve.predict - ci.values), linetype = 2) +
      geom_point(data=rawdata,aes(x=logconc,y=effect)) +
      facet_wrap(~ sampleID)
  } else {
    print("invalid plot.layout parameter, please use 'single' or 'multi'")
    return(NA)
  }
  p


}
