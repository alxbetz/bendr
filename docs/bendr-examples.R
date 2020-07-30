## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,message=FALSE,warning=FALSE-----------------------------------
library(bendr)
library(tibble)
library(dplyr)
library(ggplot2)

## ----datai,eval=FALSE----------------------------------------------------
#  fname = "/path/to/file"
#  require(readr)
#  fdata = readr::read_csv(fname)
#  fdata = readr::read_tsv(fname)
#  require(readxl)
#  fdata = readxl::read_excel(fname)
#  names(fdata)
#  View(fdata)

## ----dummydata-----------------------------------------------------------
sdata = as_tibble(
  data.frame(
    conc=c(1,10,100,1000,5000),
    effect=c(12,27,69,86,98)
  ))

## ----formula100----------------------------------------------------------
drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))

## ----formula.dyn---------------------------------------------------------
loe = min(sdata$effect)
hie = max(sdata$effect)
drc.formula.dyn = substitute(
  effect ~ bottom + ((top-bottom)/(1+10^((logEC50-logconc)*slope))),
  list(top = hie, bottom=loe))

## ----singlefit-----------------------------------------------------------

data.logged = sdata %>% mutate(logconc = log10(conc))

fitObject = bendr::fitdr(drc.formula,data.logged)
summary(fitObject$curve.fit)
pSingle = plot_single_drc(fitObject) + xlab("log(concentration)") + ylab("effect")
pSingle 
#ggsave("/path/to/image.png",pSingle,width="8",height="6",units="in")

## ----multifit------------------------------------------------------------

mdata = as_tibble(
  data.frame(
    conc=c(1,10,100,1000,5000),
    effectA=c(12,27,69,86,98),
    effectB=c(15,23,64,82,91)
  ))

fo = fitdr_multi(drc.formula,mdata,2:3,conc)
pMulti = plot_multi_drc(fo,plot.layout='multi')
pMulti

