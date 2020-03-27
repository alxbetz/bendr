
test_that("basic dose response curve fitting works", {

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

test_that("dose response curve fitting works with replicates", {

  drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
  rdata = as_tibble(
    data.frame(
      conc=c(0.1,1,10,100,1000,5000),
      replicateA=c(1,12,27,65,88,100),
      replicateB=c(2,16,21,62,81,100),
      replicateC=c(4,15,21,61,85,96)
    ))

  fitObjectRep = fitdr_replicates(drc.formula,rdata,2:4,conc)
  pRep = plot_replicate_drc(fitObjectRep)

  if(class(pRep)[2] == 'ggplot') {
    succeed()
  }
})


test_that("multiple response curve fitting works", {
  drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
  mdata = as_tibble(
    data.frame(
      conc=c(1,10,100,1000,5000),
      effectA=c(12,27,69,86,98),
      effectB=c(15,23,64,82,91)
    ))

  fo = fitdr_multi(drc.formula,mdata,2:3,conc)
  pMulti = plot_multi_drc(fo,plot.layout='multi')

  if(class(pMulti)[2] == 'ggplot') {
    succeed()
  }
})

test_that("non-toxic concentration finding works", {
  drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
  rdata = as_tibble(
    data.frame(
      conc=c(0.1,1,10,100,1000,5000),
      replicateA=c(1,12,27,65,88,100),
      replicateB=c(2,16,21,62,81,100),
      replicateC=c(4,15,21,61,85,96)
    ))

  fitObjectRep = fitdr_replicates(drc.formula,rdata,2:4,conc)
  ntc = findNTC(fitObjectRep)

  if(class(ntc) == 'numeric' && ntc < 5000 && ntc > 0.1) {
    succeed()
  }
})
