args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  stop("No file name provided\n")
}

require(readxl)
require(dplyr)
require(readr)
require(bendr)
require(ggplot2)


if(endsWith(fname,'xls') | endsWith(fname,'xls')) {
  edata = readxl::read_excel(fnames)
} else if (endsWith(fname,'csv')){
  edata = read_csv2("fname")
}

drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
nrep = ncol(edata) - 1
names(edata) <- c('conc',paste0('rep', as.character(seq(nrep))))
edataf = edata %>% dplyr::filter(conc > 0)
if(nrep > 1) {
  fo = bendr::fitdr_replicates(drc.formula,edataf,2:nrep,conc,debug=F)
} else {
  edatafLog = edataf %>% dplyr::mutate(logconc = log10(conc))
  fo = bendr::fitdr(drc.formula,edatafLog)
}

fo$ntc = findNTC(fo)
sum_idx = c(
  'ec50',
  'ec50.ci',
  'ec10',
  'ec10.ci',
  'slope',
  'slope.ci',
  'ntc',
  'confidence.level'
)
vals = unlist(fo[sum_idx])
est_vals = tibble(params = names(vals),values=vals)


bname = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(fname))
po = bendr::plot_single_drc(fo) +
  xlab('log10 Concentration [mg/L]') +
  ylab('Effect [%]') +
  theme_light
ggsave(filename = paste0(bname,'.png'),po)
saveRDS(object = fo,file = paste0(bname,'.RDS'))
write_csv(fo$plot.data,paste0(bname,"_plotdata.csv"))
write_csv(est_vals,paste0(bname,"_predicted.csv"))

