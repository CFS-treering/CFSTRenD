# author: Xiao Jing Guo (xiaojing.guo@nrcan-rncan.gc.ca)
# purpose: compare gamm models using ML method
# input: clim.test: the climate variable to test
# m.form.base: all the predictors except the one being tested
# t00: no species effect on both intercept and clim.test
# t0: no species effect on clim.test
# t1: species effect on both intercept and clim.test
# m.canidates: the list of models to be compared.

# output:
  # 1 the winning model with smallest AIC
  # 2 the fitting statistics of all involved model in m.candidates

# note: the code can be used for prediction by changing mtd to "REML" and choosing the right m.candidates


#' species-site level of detrending BAI using gamm
#' @description
#' detrending BAI using gamm model with family = Gamma(link='log') and corAR1."REML" and "ML" methods are available.It can be used to compare models with aic using "ML" method.
#'
#' @param data : datatable
#' @param mtd : method ML or REML
#' @param m.candidates : the list of RHS of formulas
#' @param out.mod : output gamm object (TRUE or FALSE)
#' @param out.csv : output csv files (TRUE or FALSE)
#' @param out.dir : directory of output csv files
#'
#'
#'
#' @import mgcv
#' @import nlme
#' @import stats
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @export detrend_site
#'

detrend_site <- function(data, mtd, m.candidates, out.mod = FALSE, out.csv = FALSE, out.dir = NULL){
  # for comparing and selecting model on AIC
  if (mtd == "ML" & length(m.candidates) > 1){
for (i in 1:length(m.candidates)){

  m.tmp <- gamm(as.formula(m.candidates[i]),
                   random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                   method = mtd,family = Gamma(link='log'), data = data)


  aic.tmp <- data.table(i = i, form = m.candidates[i], aic =  AIC(m.tmp$lme), R2 = summary(m.tmp$gam)$r.sq, methd = mtd)


  if (i == 1) {
    aic.mn <- aic.tmp$aic
    m.sel <- m.tmp
    i.sel <- i
    aic.all <- aic.tmp
  } else{
    aic.all <- rbind(aic.all, aic.tmp)
    if (aic.mn > aic.tmp$aic){
      aic.mn <- aic.tmp$aic
      m.sel <- m.tmp
      i.sel <- i
    }
  }
# if (length(m.candidates) > 1) saveRDS(m.tmp, compress = TRUE, file =   paste0(clim.test,"/","m", i, ".",mtd,".rds"))
 rm(aic.tmp, m.tmp)
}
# aic.all$clim.test <- clim.test
aic.all[aic == min(aic), selected := "*"]
form.sel <- aic.all[aic == min(aic)]$form
# saveRDS(m.sel, compress = TRUE, file =   paste0(clim.test,"/", "sel.mod", " ",mtd,".rds"))
rm(m.sel)
# saveRDS(aic.all, compress = TRUE, file =   paste0(clim.test,"/","fitting ", mtd ,".rds"))
if (out.csv == TRUE){
  if (is.null(out.dir)) message("please provide the directiory to output csv files")  else{
    if (!(dir.exists(out.dir))) dir.create(out.dir, recursive = TRUE)
write.csv(aic.all, file =  file.path(out.dir, paste0("fitting ", mtd ,".csv"), row.names = FALSE, na = ""))
  }
}
  }

  if (mtd == "ML" & length(m.candidates) == 1) {
    print(paste0("only 1 equation, NO USE for ML method, please check"))
    return()
    }
  if(mtd == "REML" & length(m.candidates) >1) {
    print(paste0("multiple equations for REML method, please check"))
    return()}
if(mtd == "REML" & length(m.candidates) == 1) {
  # for 1 equation only
  form.sel <- m.candidates
}
  # fitting REML for prediction
  m.sel <- gamm(as.formula(form.sel),
                random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                method = "REML",family = Gamma(link='log'), data = data)

  if (out.mod) saveRDS(m.sel, compress = TRUE, file =   file.path(out.dir, paste0("sel.mod_REML",".rds")))

 pred.terms  <-as.data.frame( predict(m.sel$gam, type="terms",se.fit=TRUE))
  setnames(pred.terms, c("fit.s.ageC.", "se.fit.s.ageC."),c("fit.s.ageC", "se.fit.s.ageC"))
  fit.bai <- as.data.frame(predict(m.sel$gam, type = "response", se.fit = TRUE))
  names(fit.bai) <- c("fit.bai", "se.fit.bai")
  tmp.bai <- data.table(data, pred.terms, fit.bai)
  tmp.bai$res.bai <- residuals(m.sel$gam, type = "response")

  tmp.bai[, c("lci_bai", "uci_bai") := .(
    fit.bai - qnorm(0.975) * se.fit.bai,
    fit.bai + qnorm(0.975) * se.fit.bai)]

  # names(tmp.bai)
  # fit.s.age <- "fit.s.ageC..species.fac"



    # if (any(unique(str_detect(names(tmp.bai), fit.s.age)) == TRUE)){
  #     tmp.bai[, fit.s.ageC:=eval(parse(text= paste0(fit.s.age, "LARILAR"))) + eval(parse(text= paste0(fit.s.age, "PICEMAR")))]
  #     tmp.bai[, se.fit.s.ageC:=eval(parse(text= paste0("se.", fit.s.age, "LARILAR"))) + eval(parse(text= paste0("se.",fit.s.age, "PICEMAR")))]
  # # CI
      tmp.bai[, c("lci_ageC", "uci_ageC") := .(
        fit.s.ageC - qnorm(0.975) * se.fit.s.ageC,
        fit.s.ageC + qnorm(0.975) * se.fit.s.ageC)]
    # }
  # for (iclim.test in c("VPDsummer","VPDprevsummer" )){
  # fit.s.clim <- paste0("fit.s.", iclim.test, "..species.fac")
  # if (any(unique(str_detect(names(tmp.bai), fit.s.clim)) == TRUE)){
  #   tmp.bai[, clim.fit:=eval(parse(text= paste0(fit.s.clim, "LARILAR"))) + eval(parse(text= paste0(fit.s.clim, "PICEMAR")))]
  #   tmp.bai[, se.clim.fit:=eval(parse(text= paste0("se.", fit.s.clim, "LARILAR"))) + eval(parse(text= paste0("se.", fit.s.clim, "PICEMAR")))]
  #   # CI
  #   tmp.bai[, c("lci_clim", "uci_clim") := .(
  #     clim.fit - qnorm(0.975) * se.clim.fit,
  #     clim.fit + qnorm(0.975) * se.clim.fit)]
  #
  # setnames(tmp.bai, c("clim.fit", "se.clim.fit", "lci_clim", "uci_clim"),
  #          c(paste0("fit.s.", iclim.test),paste0("se.fit.s.", iclim.test),paste0("lci_", iclim.test),paste0("uci_", iclim.test)))
  # }
  # }

  ptable <- data.table(Parameter = row.names(summary(m.sel$gam)$p.table), summary(m.sel$gam)$p.table )
  stable <- data.table(Parameter = row.names(summary(m.sel$gam)$s.table), summary(m.sel$gam)$s.table)
  aic.reml <- data.table(form = form.sel, aic =  AIC(m.sel$lme), R2 = summary(m.sel$gam)$r.sq, methd = "REML")

  if (out.csv == TRUE){
  if (is.null(out.dir)) message("please provide the directiory to output csv files")  else{
    if (!(dir.exists(out.dir))) dir.create(out.dir, recursive = TRUE)
  write.csv(ptable, file =  file.path(out.dir, paste0("ptable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(stable, file =  file.path(out.dir, paste0("stable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(tmp.bai, file = file.path(out.dir, paste0("prediction ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(aic.reml, file = file.path(out.dir, paste0("fitting ", "REML" ,".csv")), row.names = FALSE, na = "")
  }

    }
# if (make.plot) {
#   p.bai <- ggplot(tmp.bai, aes(x = year, y = fit.bai, color = species.fac)) + geom_line(size =1.2) +
#     geom_ribbon(aes(ymin=exp(lci_bai), ymax=exp(uci_bai), fill = species.fac), alpha=0.2) +
#     labs( y = paste0("bai"), x = "year") +
#     ggtitle( paste0("bai prediction") )
#
#   p.age <- ggplot(tmp.bai, aes(x = ageC, y = exp(fit.s.ageC), color = species.fac)) + geom_line(size =1.2) +
#     geom_ribbon(aes(ymin=exp(lci_ageC), ymax=exp(uci_ageC), fill = species.fac), alpha=0.2) +
#     labs( y = paste0("fit.s(", "ageC", ")"), x = "ageC") +
#     ggtitle( paste0("term response: ", "ageC") )
#
#
#
# }

# }

# return(i.sel)
return(list(model = m.sel, fitting = aic.reml, ptable = ptable, stable = stable, predtable = tmp.bai))
}

# chron <- function(data){
# dt.long <-
#
#   chron <- data[, .(N = .N, bai = mean(bai), pred.bai = mean(fit.bai), res.bai = mean(res.bai)), by = .(year)]
#
# }
