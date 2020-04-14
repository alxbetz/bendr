#' two stage nonlinear least squares fit
#'
#' @param m.formula
#' model formula
#' @param start
#' vector of start parameters
#' @param fdata
#' data frame of data to fit
#'
#' @return
#' nls2 fit object
#' @export
#'
#' @examples
nlsfit = function(m.formula,start,fdata) {
  fit1 = nlmrt::nlxb(m.formula,
              start = start,
              trace = FALSE,
              data = fdata)
  curvefit <- nls2::nls2(m.formula, data = fdata, start = fit1$coefficients,
                   algorithm = "brute-force")
  curvefit
}


#' Calculate the confidence intervals
#'
#'confidence intervals are calculated beyond the concentrations, at which measurements were taken
#'
#'
#' @param m.formula
#' model formula
#' @param x.values
#' concentration values at which the curve will be plotted
#' @param curvefit
#' model fit object from nonlinear least squares
#' @param curve.predict
#' predicted effect values at the concentrations specified in x.values
#' @param dof
#' degrees of freedom for the inverse student's t distribution ; default is (# measurements)-2
#' @param level
#' Level of confidence interval to use (0.95 by default). [0,1]
#'
#' @return
#' Confidence interval values for all concentrations in x.values
#' @export
#'
#' @examples
calc_ci = function(m.formula,x.values,curvefit,curve.predict,dof,level = 0.95,debug=FALSE) {

  logEC50.value = coef(curvefit)[1]
  slope.value = coef(curvefit)[2]

  # get the covariance matrix for the parameters
  covmatrix= vcov(curvefit)

  # get the jacobian for all x.values
  residual.jacfun = nlmrt::model2jacfun(m.formula,coef(curvefit))
  nameY = all.vars(m.formula)[1]
  nameX = all.vars(m.formula)[3]
  paramY = list(curve.predict)
  paramX = list(x.values)
  names(paramY) = nameY
  names(paramX) = nameX
  residual.jacobian = residual.jacfun(coef(curvefit), paramX, paramY)
  if(debug) {
    jac_plot = ggplot(as_tibble(residual.jacobian)) + geom_point(aes(x=logEC50,y=slope))
    ggsave("./residual_jacobian.png")
  }

  # get the confidence interval for all x_values
  error = residual.jacobian %*% covmatrix %*% t(residual.jacobian)
  error = diag(error)
  error = sqrt(error)

  # t for the inverse studen't distribution depends on the number of measurements
  p = 1- ((1-level)/2)
  t_student = qt(p, df =dof)
  ci.values = error*t_student

  ci.values
}



#' Dose Response curve fit
#'
#' @param m.formula
#' Model formula.
#' @param fdata
#' Data to fit the model to.
#' @param level
#' Level of confidence interval to use (0.95 by default). [0,1]
#' @param start
#' initial parameter values for the dose response fitting. c(logEC50,slope)
#' @param verbose
#' logical, prints intermediate results if true
#'
#' @return
#' @importFrom nlstools confint2
#' @export
#'
#' @examples
fitdr = function(m.formula,fdata,level=0.95,start=vector(),verbose=FALSE,debug=FALSE) {
  responseName = all.vars(m.formula)[1]
  concName = all.vars(m.formula)[3]
  if(length(start) == 0) {
    #start = c(logEC50 = median(fdata$logconc), slope = diff(range(fdata$logconc))/ diff(range(fdata$effect)))
    responseV = fdata %>% dplyr::pull(responseName)
    dependentV = fdata %>% dplyr::pull(concName)
    #set slopeGuess to +1 if the response values are ascending and -1 if they are descending
    #slopeGuess = sign(mean(diff(responseV)))
    slopeGuess = as.double(sign(responseV[length(responseV)]-responseV[1]))
    #if the concentration values are descending in the dataframe, negate the slope guess
    #if(mean(diff(dependentV)) < 0)
    if(dependentV[length(dependentV)] - dependentV[1] < 0 ) {
      slopeGuess = - slopeGuess
    }
    start = c(logEC50 = median(dplyr::pull(fdata,concName)), slope = slopeGuess)
    if(verbose) {
      print("Initial values are")
      print(start)
    }

  }
  curvefit = nlsfit(m.formula,start=start,fdata = fdata)
  x.values = calc_xrange(m.formula,fdata)
  xdf = data.frame(x.values)
  names(xdf) = concName
  curve.predict = predict(curvefit, newdata = xdf)

  ci.values = calc_ci(m.formula,x.values,curvefit,curve.predict,dof=nrow(fdata)-2,level=level,debug=debug)
  plot.data = data.frame(log.concentration = x.values,curve.predict = curve.predict,ci.values = ci.values)
  logec50 = coefficients(curvefit)[1]
  ec50 = 10^logec50
  slope = coefficients(curvefit)[2]
  ec10 = calc_ecx(10,slope,ec50)
  logec10 = log10(ec10)
  ci.par = confint2(curvefit)
  ec50.ci = 10^ci.par['logEC50',]
  slope.ci = ci.par['slope',]
  ec10.ci = c(
    calc_ecx(10,ci.par['slope',1],ec50.ci[1]),
    calc_ecx(10,ci.par['slope',2],ec50.ci[2])
  )


  names(ec50) = c('')
  names(ec10) = c('')
  list(plot.data = plot.data ,
       curve.fit=curvefit,
       confidence.level=level,
       data=fdata,
       ec50=ec50,
       ec50.ci = ec50.ci,
       ec10 = ec10,
       ec10.ci = ec10.ci,
       slope = slope,
       slope.ci = slope.ci,
       aic = AIC(curvefit),
       xname = concName,
       yname = responseName
  )
}


#' Dose response curve fit with replicates
#'
#' @param m.formula
#' @param fdata
#' @param effectColumns
#' @param concColumn
#' @param level
#' @param start
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
fitdr_replicates= function(m.formula,fdata,effectColumns,concColumn,level=0.95,start=vector(),verbose=F,debug=F) {
  quo_concColumn = enquo(concColumn)
  data.logged = fdata %>% dplyr::mutate(logconc := log10(!! quo_concColumn))
  data.long = data.logged %>% tidyr::gather(replicateID,effect,effectColumns) %>%
    arrange(desc(logconc))
  #data.summ = data.long %>% dplyr::group_by(logconc) %>% dplyr::summarise(mean = mean(effect,na.rm=T),n=n(),sd = sd(effect,na.rm=T))

  #m.formula = as.formula(stringr::str_replace(deparse(m.formula),'effect','mean'))

  data.fit = fitdr(m.formula,data.long,level=level,start=start,verbose=verbose,debug=debug)
  data.fit$data.long = data.long
  data.fit$nreplicates = length(effectColumns)
  data.fit
}



#' Fit multiple dose response curves
#'
#' @param m.formula
#' model formula
#' @param fdata
#' data frame with one column corresponding to the (non-logged) concentration values
#' and an arbitrary number of columns corresponding to effects
#' @param effectColumns
#' vector of effect columns; can be either numeric indices or column names
#' @param concColumn
#' name of the concentration column; this parameter is lazily evaluated, so do not use quotes
#'
#' @return
#' @import rlang
#' @import dplyr
#' @export
#'
#' @examples
fitdr_multi= function(m.formula,fdata,effectColumns,concColumn) {
  quo_concColumn = enquo(concColumn)
  data.logged = fdata %>% dplyr::mutate(logconc := log10(!! quo_concColumn))
  data.long = data.logged %>% tidyr::gather(sampleID,effect,effectColumns)
  data.fit = data.long %>% dplyr::group_by(sampleID) %>% tidyr::nest() %>%
    dplyr::mutate(m.fit = purrr::map(data,fitdr,m.formula = m.formula)) %>%
    dplyr::mutate(
      ec50=purrr::map_dbl(m.fit,function(mod) 10^(coefficients(mod$curve.fit)[1])),
      slope=purrr::map_dbl(m.fit,function(mod) coefficients(mod$curve.fit)[2]),
      aic=purrr::map_dbl(m.fit,function(mod) mod$aic))
  data.fit
}
