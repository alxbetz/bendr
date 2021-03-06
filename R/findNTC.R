#' Find the non-toxic concentration.
#'
#' If you use the result of this calculation in your research, please cite Stadnicka-Michalak et al.
#'
#' @param bendrObj
#' A bendr result object
#'
#' @return
#' The absolute (non-logged) non-toxic concentration value
#' @export
#'
#' @examples
#' require(tibble)
#'drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
#'rdata = as_tibble(
#'data.frame(
#'    conc=c(0.1,1,10,100,1000,5000),
#'    replicateA=c(1,12,27,65,88,100),
#'    replicateB=c(2,16,21,62,81,100),
#'    replicateC=c(4,15,21,61,85,96)
#'  ))
#'
#'fitObjectRep = fitdr_replicates(drc.formula,rdata,2:4,conc)
#'ntc = findNTC(fitObjectRep)
#'
findNTC = function(bendrObj) {
  if(!any(bendrObj$plot.data$ci.values.lower > 90)) {
    plot.data = extend_ci(bendrObj)
  } else {
    plot.data = bendrObj$plot.data
  }

  upperCI = plot.data$ci.values.upper
  lowerCI = plot.data$ci.values.lower
  curve.predict = plot.data$curve.predict
  log.concentration = plot.data$log.concentration
  logEC50 = coef(bendrObj$curve.fit)['logEC50']
  slope = coef(bendrObj$curve.fit)['slope']

  # In case of very narrow confidence intervals, caused e.g. by the non-optimal concentration range:
  #   If Effect_NtClower CI> 10 * Effect_NtCupperCI -> Effect_NtC = Effect_NtClowerCI / 10 where: Effect_NtCi (%) is the effect caused by the NtC
  # NtC can be calculated then from Effect_NtC (and sigmoidal dose-response curve).
  #


  ## NtC upper
  # find concentration at which the upper confidence interval falls below
  # 100
  uci100 = upperCI > 100
  uci100idx = which(uci100)
  if(length(uci100idx) > 0) {
    inds = max(uci100idx)
    Effect_NtC_upper = 100 - curve.predict[inds]
    if (Effect_NtC_upper > 0){
      NtC_upper = 10^(log.concentration[inds])
    } else {
      NtC_upper = 0
    }
  } else {
    Effect_NtC_upper = 0
    NtC_upper = 0
  }



  ## NtC lower
  # find concentration at which lower confidence interval is 90
  inds = max(which(lowerCI > 90))
  Effect_NtC_lower = 100 - curve.predict[inds]
  if (Effect_NtC_lower > 0){
    NtC_lower = 10^(log.concentration[inds])
  } else {
    NtC_lower = 0
  }




  ## Choose NtC based on upper and lower CI criteria
  if (NtC_upper > NtC_lower) {
    NtC = NtC_lower
    Effect_NtC = Effect_NtC_lower
  } else {


    if (Effect_NtC_lower/Effect_NtC_upper > 10) {
      Effect_NtC = Effect_NtC_lower/10
      NtC = 10^(logEC50-log10(100/(100-Effect_NtC)-1)/slope)

    } else {
      NtC = NtC_upper
      Effect_NtC = Effect_NtC_upper
    }
  }

  ## Include measured data if the effect of any replicate (at concentration below
  # the concentration, which has been thus far chosen as the NtC) is > 10%
  # 1. which measured concentrations are below NtC
  # 2. does any of the effects measured at these concentrations have an effect > 10%
  # 3. if yes, select the lowest concentration that is below the NtC


  belowNtC = bendrObj$data %>%
    dplyr::filter(conc < NtC,effect <= 90)
  if(nrow(belowNtC) > 0) {
    NtC = belowNtC %>% pull(conc) %>% min()
  }

  NtC
}


