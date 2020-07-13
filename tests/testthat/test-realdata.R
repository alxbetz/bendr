require(tidyverse)
require(bendr)

rawfiles = list.files(path = "~/polybox2/dr_database/rawdata/newnames",full.names = T)
rawdata = lapply(rawfiles,readxl::read_excel)
rawdata_clean = lapply(rawdata,cleanAndPivot)
drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))

fitObject = bendr::fitdr(drc.formula,rawdata_clean[[1]],ci_method = "bootstrap")
fitHelper = function(x) {
  if(!is.na(x)) {
    bendr::fitdr(drc.formula,x,ci_method = "delta")
  } else {
    NA
  }

}
outdir = "~/polybox2/dr_database/rfit_plots"
fits <- vector(mode = "list", length = length(rawdata))
for(i in 1:length(rawdata_clean)) {


  print(i)
  bname = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(rawfiles[[i]]))
  outfile = file.path(outdir,paste0(bname,".png"))
  if(any(is.na(rawdata_clean[[i]]))) {
    next
  }
  qd = qdiff(rawdata_clean[[i]][['effect']],lo = 0.05,hi = 0.95)
  minConcEffect = rawdata_clean[[i]] %>% slice(which.min(.$conc)) %>% pull(effect)
  maxConcEffect = rawdata_clean[[i]] %>% slice(which.max(.$conc)) %>% pull(effect)
  rawdata_clean[[i]] %>% slice(which.min(.$conc)) %>% pull(effect)


  if(
      (
        qd<50 &&
        (
          max(rawdata_clean[[i]][['effect']]) < 53 ||
          min(rawdata_clean[[i]][['effect']]) > 47
        )
      ) ||
     min(rawdata_clean[[i]][['effect']]) > 60 ||
     max(rawdata_clean[[i]][['effect']]) < 40
    ) {
    print(qd)

    p = ggplot(rawdata_clean[[i]]) + geom_point(aes(x=logconc,y=effect)) +
      ggtitle(paste0(bname,"|| qdiff: ",qd))

  } else if(minConcEffect < maxConcEffect){
    #print("Effect at minimum concentration is lower than at maximum concentration. Is the data correctly preprocessed ?")

    p = ggplot(rawdata_clean[[i]]) + geom_point(aes(x=logconc,y=effect)) +
      ggtitle(bname)
    ggsave(p,filename = outfile)
    } else{
    fits[[i]] = fitHelper(rawdata_clean[[i]])
    p = bendr::plot_single_drc(fits[[i]])

    }
  p = p +
    scale_x_continuous('log10 Concentration [mg/L]') +
    scale_y_continuous('% Cell Viability',breaks = seq(0,120,by=20), labels = as.character(seq(0,120,by=20)),limits = c(-10,120)) +
    theme_bw()
  ggsave(p,filename = outfile)
}

#toreview
#"/Users/alx/polybox2/dr_database/rawdata/newnames//LIV_kuehen_48_0024h_AB_m1_FBS00_0.5mL_nn_in_no_W_525-66-6.xlsx"
#1109
#"/Users/alx/polybox2/dr_database/rawdata/newnames//LIV_kuehen_48_0024h_NR_m1_FBS00_0.5mL_nn_di_no_E_15307-86-5.xlsx"
#1125
#"/Users/alx/polybox2/dr_database/rawdata/newnames//LIV_kuehen_48_0024h_NR_m2_FBS50_0.5mL_nn_di_no_E_15307-86-5.xlsx"
#1201
fo = fitHelper(rawdata_clean[[1109]])
bendr::plot_single_drc(fo)

j=59
qdiff(rawdata_clean[[j]]['effect'],lo = 0.05,hi = 0.95)
ggplot(rawdata_clean[[67]]) + geom_point(aes(x=logconc,y=effect))
hf = function(x,lo,hi){
  if(is.na(x)) {
    return(NA)
  }
  qdiff(x[["effect"]],lo,hi)
}
rawdata_qd = lapply(rawdata_clean,hf,lo=0.05,hi=0.95)

fitted = lapply(rawdata_clean,fitHelper)

bendr::plot_single_drc(fitObject)
fitHelper(rawdata_clean[[128]])
test_that("real data bootstrapping works", {

  drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
  sdata = as_tibble(
    data.frame(
      conc=c(0.0001,0.001,0.01,0.1,1),
      effect=c(12,27,69,86,98)
    ))
  data.logged = sdata %>% mutate(logconc = log10(conc))

  fitObject = fitdr(drc.formula,data.logged)
  summary(fitObject$curve.fit)
  pSingle = plot_single_drc(fitObject) + xlab("log(concentration)") + ylab("effect")
  if(class(pSingle)[2] == 'ggplot') {
    succeed()
  }
})

###ntc
ntcoutdir = file.path(outdir,"ntc")
ntc_v = rep(NA,length(rawdata))
for (i in 1:293) {
  print(i)
  if(length(fits[[i]]) > 0) {
    ntc_v[i] = bendr::findNTC(fits[[i]])
    bname = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(rawfiles[[i]]))
    p = bendr::plot_single_drc(fits[[i]]) + geom_vline(xintercept = log10(bendr::findNTC(fits[[i]])))
    ggsave(filename = file.path(ntcoutdir,paste0(bname,".png")),plot = p)
    }


  }

drc.inverse.function = function(effect,logEC50,slope) {
  return(logEC50 - (log10((100/effect)-1)/slope))
}

drc.inverse.function(90,ci.lower.coef[1],ci.lower.coef[2])

i=247

lox = drc.inverse.function(99.9,log10(fits[[i]]$ec50),fits[[i]]$slope)

xtend_x = seq(from = lox,
              to = max(fits[[i]]$plot.data$log.concentration),
              length.out = 1E3)

curve.predict.xtend = predict(fits[[i]]$curve.fit,newdata = data.frame(logconc=xtend_x))
ci.values.xtend = bendr::prediction_ci(drc.formula,
                     x.values = xtend_x,
                     curvefit = fits[[i]]$curve.fit,
                     curve.predict = curve.predict.xtend,
                     dof = nrow(fits[[i]]$data)-2)
plot(curve.predict.xtend-ci.values.xtend)


for (i in 1:length(rawdata)) {

  if(length(fits[[i]]) > 0) {
    if(!any(fits[[i]]$plot.data$ci.values.upper > 100)) {
      print(i)
      print(bendr::findNTC(fits[[i]]))
    }

#    ntc_v[i] = bendr::findNTC(fits[[i]])

  }


}
