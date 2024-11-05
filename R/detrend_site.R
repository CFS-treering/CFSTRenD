# purpose: compare gamm models using ML method


#' species-site level of detrending BAI using gamm
#' @description
#' detrending BAI using gamm model with family = Gamma(link='log') and corAR1."REML" and "ML" methods are available.It can be used to compare models with aic using "ML" method.
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.
#' @param out.csv directory of output csv files. default is NULL for no csv output.
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
#' @details
#' If users specify multiple candidate models as input to m.candiates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model. Once the optimal model is identified,
#' it will be refitted using the restricted maximum likelihood (REML) method to obtain unbiased estimates of model parameters.
#'
#' If users specify only 1 candidate model as input to m.candiates argument, the model will apply REML model directly.

#' @export detrend_site
#'

detrend_site <- function(data, resp_scale = "log", m.candidates,  out.csv = NULL){
  if (resp_scale == "log") {

    famil = gaussian("identity")
  }else {
    famil = Gamma("log")
  }
  # for comparing and selecting model on AIC
  if (length(m.candidates) > 1){
for (i in 1:length(m.candidates)){
  formul <- as.formula(m.candidates[i])
  if (resp_scale == "log") formul <- update(formul, log(.) ~ .)

  m.tmp <- gamm(formul,
                   random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                   method = "ML",family = famil, data = data)


  aic.tmp <- data.table(i = i, form = gsub("\\\\", "", paste(deparse(formul), collapse = " ")), aic =  AIC(m.tmp$lme), R2 = summary(m.tmp$gam)$r.sq, methd = "ML")


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
if (!is.null(out.csv)){

  if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
write.csv(aic.all, file =  file.path(out.csv, paste0("fitting ML.csv")), row.names = FALSE, na = "")
  }


  }

  # if (mtd == "ML" & length(m.candidates) == 1) {
  #   print(paste0("only 1 equation, NO USE for ML method, please check"))
  #   return()
  #   }
  # if(mtd == "REML" & length(m.candidates) >1) {
  #   print(paste0("multiple equations for REML method, please check"))
  #   return()}
if(length(m.candidates) == 1) {
  # for 1 equation only
  form.sel <- as.formula(m.candidates)
  if (resp_scale == "log") form.sel <- update(form.sel, log(.) ~ .)
}
  # fitting REML for prediction
  m.sel <- gamm(form.sel,
                random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                method = "REML",family = famil, data = data)


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
      if (resp_scale == "log") setnames(tmp.bai, c("fit.bai", "se.fit.bai", "res.bai", "lci_bai", "uci_bai"), c("fit.logbai", "se.fit.logbai", "res.logbai", "lci_logbai", "uci_logbai"))

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
  aic.reml <- data.table(form = gsub("\\\\", "", paste(deparse(form.sel), collapse = " ")), aic =  AIC(m.sel$lme), R2 = summary(m.sel$gam)$r.sq, methd = "REML")

  if (!is.null(out.csv)){

    if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
  write.csv(ptable, file =  file.path(out.csv, paste0("ptable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(stable, file =  file.path(out.csv, paste0("stable ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(tmp.bai, file = file.path(out.csv, paste0("prediction ", "REML" ,".csv")), row.names = FALSE, na = "")
  write.csv(aic.reml, file = file.path(out.csv, paste0("fitting ", "REML" ,".csv")), row.names = FALSE, na = "")
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
return(list(model = m.sel, fitting = aic.reml, ptable = ptable, stable = stable, pred = tmp.bai))
}

# chron <- function(data){
# dt.long <-
#
#   chron <- data[, .(N = .N, bai = mean(bai), pred.bai = mean(fit.bai), res.bai = mean(res.bai)), by = .(year)]
#
# }
