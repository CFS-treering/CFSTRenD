#' a general gam model

#' @description
#' gam model benefiting the model selection procedure using ML method
#'
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is original scale with family = Gamma("log"); "log" for predicting log-scale with family = gaussian("identity") ; "hybrid" for predicting original scale with family = gaussian("identity")
#' @param m.candidates the list of formulas.
#'
#'
#' @import mgcv
#' @import nlme
#' @importFrom MuMIn AICc
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @details
#' This function models a gam model without considering random effects or autocorrelation.

#' If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model. Once the optimal model is identified,
#' it will be refitted using the restricted maximum likelihood (REML) method and output the results.
#'
#' If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.

#' @export gam_mod

gam_mod <- function(data, resp_scale = "", m.candidates){

  gamm_main(data , resp_scale , m.option = 0, m.candidates)

}



#' detrending model on tree-ring width series

#' @description
#' models the biological growth trends in individual tree-ring width series using mgcv::gamm
#'
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.
#'
#'
#' @import mgcv
#' @import nlme
#' @import stats
#' @importFrom MuMIn AICc
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @details
#' This function models the biological growth trends in individual tree-ring width series using mcgv::gamm.
#' By integrating a first-order autoregressive (AR1) component, it accounts for temporal autocorrelation.
#' This method can provide  “normalized” residuals, which are adjusted to reflect deviations after considering the AR1 correlation structure.
#' 'Normalized' residuals are valuable for further analyses, such as investigating relationships with climatic variables.

#' If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model. Once the optimal model is identified,
#' it will be refitted using the restricted maximum likelihood (REML) method and output the results.
#'
#' If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.

#' @export gamm_radius

gamm_radius <- function(data, resp_scale = "resp_gamma", m.candidates){

  gamm_main(data , resp_scale , m.option = 1, m.candidates)

}




#' growth model at site-level
#' @description
#' models the growth trend or climate-growth relationship at site-level
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.
#'
#'
#'
#' @import mgcv
#' @import nlme
#' @importFrom MuMIn AICc
#' @import stats
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @details
#' This function accounts for within-site variability and temporal autocorrelation by including series identity as random effects
#' and a first-order autoregressive (AR1) correlation structures, respectively. “Normalized” residuals are also provided for future analysis.
#'
#' If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model.
#' Once the optimal model is identified, it will be refitted using the restricted maximum likelihood (REML) method and output the results.
#'
#' If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.


#' @export gamm_site
#'
#'
gamm_site <- function(data, resp_scale = "resp_gamma", m.candidates){

  data$uid_tree.fac <- as.factor(as.character(data$uid_tree))
  gamm_main(data , resp_scale , m.option = 2, m.candidates )
}



#' spatial growth model at regional-level (multiple sites)
#' @description
#' models the growth trend or climate-growth relationship at regional-level with multiple sites
#'
#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.

#'
#'
#'
#' @import mgcv
#' @import nlme
#' @import stats
#' @importFrom MuMIn AICc
#' @import utils
#' @import data.table
#'
#'
#' @return list including model, fitting statistics, ptable, stable, prediction table and spatial effect(moranI)
#' @details
#' This function accounts for within-site variability and temporal autocorrelation by including series identity as random effects
#' and a first-order autoregressive (AR1) correlation structures, respectively. Among-site variability and spatial effects are captured by incorporating site identity as random effects.
#' The model is refitted automatically by introducing a smooth term for latitude and longitude using the Spatial Over-Smooth ("sos") basis if significant spatial autocorrelation persists.
#' “Normalized” residuals are provided for future analysis.
#'
#' If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model.
#' Once the optimal model is identified, it will be refitted using the restricted maximum likelihood (REML) method and output the results.
#'
#' If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.

#' @export gamm_spatial

gamm_spatial <- function(data, resp_scale = "resp_gamma", m.candidates){
  data$uid_tree.fac <- as.factor(as.character(data$uid_tree))
  data$uid_site.fac <- as.factor(as.character(data$uid_site))
  m.candidates <- paste0(m.candidates, " + s(uid_site.fac, bs = 're')")
  gamm_main(data , resp_scale , m.option = 3, m.candidates)

}


#' spatial growth model for large dataset or vast geographical coverage
#' @description
#' To address the computational limitations of GAMMs for large datasets, this function offers a hybrid solution combining the efficiency of the mgcv::bam() function
#' with the itsadug package's post hoc methods for managing temporal autocorrelation.
#'

#' @param data data containing all necessary columns to run the model
#' @param resp_scale the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is
#' @param m.candidates the list of formulas.
#'
#'
#'
#' @import furrr
#' @import parallel

#' @import mgcv
#' @import itsadug
#' @import nlme
#' @import stats
#' @importFrom MuMIn AICc
#' @import utils
#' @import data.table
#' @import sp
#' @import spdep
#'
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @details
#' This function accounts for within-site variability and temporal autocorrelation by including series identity as random effects
#' and a first-order autoregressive (AR1) correlation structures, respectively. Among-site variability and spatial effects are captured by incorporating site identity as random effects.
#' The model is refitted automatically by introducing a smooth term for latitude and longitude using the thin plate ("tp") basis if significant spatial autocorrelation persists.
#' “Normalized” residuals are provided for future analysis.
#'
#' This function supports parallel computation for the large-scale, geographically distributed datasets.
#'
#' If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
#' The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model.
#' Once the optimal model is identified, it will be refitted using the restricted maximum likelihood (REML) method and output the results.
#'
#' If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.

#' @export bam_spatial

bam_spatial <- function(data, resp_scale = "resp_gamma", m.candidates){
  data$uid_tree.fac <- as.factor(as.character(data$uid_tree))
  data$uid_site.fac <- as.factor(as.character(data$uid_site))
  m.candidates <- paste0(m.candidates, " + s(uid_site.fac, bs = 're')")
  gamm_main(data , resp_scale , m.option = 4, m.candidates)

}



# main function
gamm_main <- function(data, resp_scale = "resp_gamma", m.option, m.candidates){
  if (length(m.candidates) == 0) stop("must assign the equation(s) to m.candidates")
  if ( !(resp_scale %in% c("resp_log", "resp_gamma", "resp_gaussian"))) stop(paste0( "please check resp_scale, it allows 3 options: resp_log, on log-scale of resp; resp_gaussian, on response scale assuming gaussian distribution; resp_gamma, on response scale with family Gamma(log)"))
#   if (resp_scale == "resp_log") {
# # if in log-scale, use gaussian distribution
#     famil = gaussian("identity")
#   }else {
#     # in orginal scale, use Gamma with log link function
#     famil = Gamma("log")
#     if (resp_scale == "resp_gaussian") {
#
#       famil = gaussian("identity")
#     }
#   }

  if (resp_scale == "resp_log") famil = gaussian("identity")
  if (resp_scale == "resp_gaussian") famil = gaussian("identity")
  if (resp_scale == "resp_gamma") famil = Gamma("log")

  # for comparing and selecting model on AIC
  if (length(m.candidates) > 1){
    for (i in 1:length(m.candidates)){
      formul <- as.formula(m.candidates[i])
      # in log-scale, log-transfrom response variable
      if (resp_scale == "resp_log") formul <- update(formul, log(.) ~ .)

      if (m.option == 0){

        m.tmp <- gam(formul,
                      method = "ML",family = famil, data = data)

      }
      if (m.option == 1){

        m.tmp <- gamm(formul,
                      correlation = corCAR1(value = 0.5),
                      method = "ML",family = famil, data = data)

      }
      if (m.option %in% c(2,3)){
        m.tmp <- gamm(formul,
                      random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                      method = "ML",family = famil, data = data)
      }

      if (m.option == 4)  {
               # Detect available cores for parallel processing
        available_cores <- detectCores(logical = FALSE) - 1  # Adjusted cores based on system

        # Decide if parallel processing is supported
        if (available_cores > 1) {
          plan(multisession, workers = available_cores)
        } else {
          plan(sequential)
        }

        m0.tmp <- bam(formul,

                      method = "ML",data = data)
        setorder(data, uid_site, uid_tree, ageC)
        data[, start.event := c(TRUE, rep(FALSE, .N - 1)), by = .(uid_site, uid_tree)]
        # r1 <- start_value_rho(m0, plot=TRUE)
        # print (paste0(Sys.time(), "          ar1"))
        m.tmp <- bam(formul, data=data, rho=start_value_rho(m0.tmp, plot=FALSE), AR.start=data$start.event, method = "ML")

        rm(m0.tmp)
        # Reset to sequential
        plan(sequential)

      }
      # if (m.option < 4){

      # aic.tmp <- data.table(form = gsub("\\\\", "", paste(deparse(formul), collapse = " ")), R2 = summary(m.tmp)$r.sq, methd = "ML")
      # if (m.option == 0) {
      #   aic.tmp$aic <-  AIC(m.tmp)
      #   aic.tmp$aicc <-  AICc(m.tmp)
      #   aic.tmp$bic <-  BIC(m.tmp)}else {
      #     aic.tmp$aic <-  AIC(m.tmp$lme)
      #     aic.tmp$aicc <-  AICc(m.tmp$lme)
      #     aic.tmp$bic <-  BIC(m.tmp$lme)
      #
      #   }
      # }else{

        aic.tmp <- data.table(i = i, form = gsub("\\\\", "", paste(deparse(formul), collapse = " ")), aic =  AIC(m.tmp), aicc = AICc(m.tmp), bic =  BIC(m.tmp), R2 = summary(m.tmp)$r.sq, methd = "ML")

    # }

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
      rm(aic.tmp, m.tmp)
    }
    aic.all[aicc == min(aicc), selected := "*"]
    form.sel <- as.formula(aic.all[aicc == min(aicc)]$form)
    rm(m.sel)
    # if (!is.null(out.csv)){
    #
    #   if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
    #   write.csv(aic.all, file =  file.path(out.csv, paste0("fitting ML.csv")), row.names = FALSE, na = "")
    # }


  }else{
    # for 1 equation only
    form.sel <- as.formula(m.candidates)
    if (resp_scale == "resp_log") form.sel <- update(form.sel, log(.) ~ .)
  }
  # fitting REML for prediction

  if (m.option == 0) {
    m.sel <- gam(form.sel,
                method = "REML",family = famil, data = data)

  }

  if (m.option == 1) {
    m.sel <- gamm(form.sel,
                  correlation = corCAR1(value = 0.5),
                  method = "REML",family = famil, data = data)

  }
   if (m.option %in% c(2,3)) {
     m.sel <- gamm(form.sel,
                  random = list(uid_tree.fac=~1), correlation = corCAR1(value = 0.5),
                  method = "REML",family = famil, data = data)
    data[, res.normalized:=residuals(m.sel$lme, type = "normalized")]
    # Extract the substring inside parentheses in case in log-scale
    y.char <- sub(".*\\((.*?)\\).*", "\\1", all.vars(m.sel$gam$formula)[1])

    }
  if (m.option == 4){


      # Detect available cores for parallel processing
      available_cores <- detectCores(logical = FALSE) - 1  # Adjusted cores based on system

      # Decide if parallel processing is supported
      if (available_cores > 1) {
        plan(multisession, workers = available_cores)
      } else {
        plan(sequential)
      }

      m0.sel <- bam(form.sel,

                    method = "fREML",data = data)
      setorder(data, uid_site, uid_tree, ageC)
      data[, start.event := c(TRUE, rep(FALSE, .N - 1)), by = .(uid_site, uid_tree)]

      # Step 3: Add rho and AR1 structure
      rho.start <- start_value_rho(m0.sel, plot = FALSE)
      m.sel <- bam(
        form.sel,
        data = data,
        rho = rho.start,
        AR.start = data$start.event,
        method = "fREML"
      )

      data[start.event==FALSE, res.normalized:=resid_gam(m.sel)]
      # Extract the substring inside parentheses in case in log-scale
      y.char <- sub(".*\\((.*?)\\).*", "\\1", all.vars(m.sel$formula)[1])

      # Reset to sequential
      plan(sequential)


  }

    # further for spatial effect
    if (m.option >= 3 & all(c("latitude", "longitude") %in% names(data))){

      data[, c("lon_use", "lat_use") := list(round(longitude,1), round(latitude,1) )]

      y.site <- data[!is.na(res.normalized), .(obs.med = median(get(y.char)), res.med = median(res.normalized)), by = .(lon_use, lat_use )]
      if (nrow(y.site) > 5) {
      coordinates(y.site) <- ~ lon_use+lat_use
      knea <- knearneigh(coordinates(y.site), longlat = TRUE)
      nb <- knn2nb(knea)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      moran.I.o <- moran.test(y.site$obs.med,lw)
      moran.I.o <- data.frame(statistic = moran.I.o$estimate[1], expected = moran.I.o$estimate[2], Variance = moran.I.o$estimate[3],
                              p.value = moran.I.o$p.value)
      names(moran.I.o) <- paste0(names(moran.I.o), ".obs")

      moran.I.r <- moran.test(y.site$res.med,lw)
      moran.I.r <- data.frame(statistic = moran.I.r$estimate[1], expected = moran.I.r$estimate[2], Variance = moran.I.r$estimate[3],
                              p.value = moran.I.r$p.value)
      names(moran.I.r) <- paste0(names(moran.I.r), ".res")
      moranI <- cbind(moran.I.o, moran.I.r)
      rm(nb, knea, moran.I.o, moran.I.r)
      if (moranI$p.value.obs < 0.1 & moranI$p.value.res < 0.1 ) {
      if (m.option == 3)  {
        form.sel <- update(m.sel$gam$formula, . ~ . + s(latitude, longitude, bs = "sos"))
        m.sel <- gamm(form.sel,  data = data,
                      random = list(uid_tree.fac= ~1),correlation =  corCAR1(value = 0.5),
                      family = famil)
        data[, res.normalized.LL:=residuals(m.sel$lme, type = "normalized")]
      }
        if (m.option == 4)  {
          form.sel <- update(m.sel$formula, . ~ . + s(latitude, longitude, bs = "tp", k = 5))
          # Detect available cores for parallel processing
          available_cores <- detectCores(logical = FALSE) - 1  # Adjusted cores based on system

          # Decide if parallel processing is supported
          if (available_cores > 1) {
            plan(multisession, workers = available_cores)
          } else {
            plan(sequential)
          }

          m0.sel <- bam(form.sel,

                        method = "fREML",data = data)
          setorder(data, uid_site, uid_tree, ageC)
          # data[, idrow:=seq_len(.N), by = .(uid_site, uid_tree)][,start.event := (idrow== 1)][, idrow:= NULL]
          # r1 <- start_value_rho(m0, plot=TRUE)
          # print (paste0(Sys.time(), "          ar1"))
          m.sel <- bam(form.sel, data=data, rho=start_value_rho(m0.sel, plot=FALSE), AR.start=data$start.event)

          data[start.event==FALSE, res.normalized.LL:=resid_gam(m.sel)]

          # Reset to sequential
          plan(sequential)

        }

        y.site <- data[!is.na(res.normalized.LL), .(res.med.LL = median(res.normalized.LL)), by = .(lon_use, lat_use )]

        moran.I.r_LL <- moran.test(y.site$res.med.LL,lw)
        moran.I.r_LL <- data.frame(statistic = moran.I.r_LL$estimate[1], expected = moran.I.r_LL$estimate[2], Variance = moran.I.r_LL$estimate[3],
                                p.value = moran.I.r_LL$p.value)
        names(moran.I.r_LL) <- paste0(names(moran.I.r_LL), ".res_LL")
        moranI <- cbind(moranI, moran.I.r_LL)

      }
    }
    } # m.option >= 3, for testing spatial effect, and adding s(lat, lon) if necessary
  # }



  if (m.option < 4){
    if (m.option == 0) msel.gam <- m.sel else msel.gam <- m.sel$gam
  pred.terms  <-as.data.frame( predict(msel.gam, type="terms",se.fit=TRUE))
  rhs_terms <- attr(terms(msel.gam$formula), "term.labels")
  # to set column names for only 1 term,
  if (length(rhs_terms) == 1){
  if (length(names(pred.terms)) == 2) # only for 1 term with fit and se.fit
  names(pred.terms) <-paste0(c("fit.", "se.fit."), sub("s\\(([^,\\)]+).*", "\\1", rhs_terms))
  }

  pred.terms <- format_byterm(msel.gam, pred.terms)

  fit.y <- as.data.frame(predict(msel.gam, type = "response", se.fit = TRUE))
  names(fit.y) <- c("fit.resp", "se.fit.resp")
  tmp.y <- data.table(data, pred.terms, fit.y)
  tmp.y$res.resp <- residuals(msel.gam, type = "response")
  if (m.option > 0) tmp.y$res.resp_normalized <- residuals(m.sel$lme, type = "normalized")
  if ("res.normalized" %in% names(tmp.y)) tmp.y$res.normalized <- NULL
  # if (resp_scale == "resp_log") setnames(tmp.y, c("fit.resp", "se.fit.resp", "res.resp", "res.resp_normalized"),
  #                                   c("fit.log_resp", "se.fit.log_resp", "res.log_resp", "res.log_resp_normalized"))


  ptable <- data.table(Parameter = row.names(summary(msel.gam)$p.table), summary(msel.gam)$p.table )
  stable <- data.table(Parameter = row.names(summary(msel.gam)$s.table), summary(msel.gam)$s.table)
  aic.reml <- data.table(form = gsub("\\\\", "", paste(deparse(form.sel), collapse = " ")), aic =  AIC(m.sel), aicc = AICc(m.sel), bic =  BIC(m.sel), R2 = summary(msel.gam)$r.sq,  methd = "REML")
  # if (m.option == 0) {
  #   aic.reml$aic <-  AIC(msel.gam)
  #   aic.reml$aicc <- AICc(msel.gam)
  #   aic.reml$bic <-  BIC(msel.gam)
  #
  #   }else {
  #     aic.reml$aic <-  AIC(m.sel$lme)
  #     aic.reml$aicc <- AICc(m.sel$lme)
  #     aic.reml$bic <-  BIC(m.sel$lme)
  #
  # }

  }else{
    pred.terms  <-as.data.frame( predict(m.sel, type="terms",se.fit=TRUE))
    fit.y <- as.data.frame(predict(m.sel, type = "response", se.fit = TRUE))
    names(fit.y) <- c("fit.resp", "se.fit.resp")
    tmp.y <- data.table(data, pred.terms, fit.y)
    tmp.y$res.resp <- residuals(m.sel, type = "response")
    tmp.y[start.event==FALSE,res.resp_normalized := resid_gam(m.sel)]




    ptable <- data.table(Parameter = row.names(summary(m.sel)$p.table), summary(m.sel)$p.table )
    stable <- data.table(Parameter = row.names(summary(m.sel)$s.table), summary(m.sel)$s.table)
    aic.reml <- data.table(form = gsub("\\\\", "", paste(deparse(form.sel), collapse = " ")), aic =  AIC(m.sel), aicc = AICc(m.sel), bic = BIC(m.sel), R2 = summary(m.sel)$r.sq, methd = "fREML")

  }
  # rename log-scale
  if (resp_scale == "resp_log") {
    setnames(tmp.y, c("fit.resp", "se.fit.resp", "res.resp"),
                                    c("fit.log_resp", "se.fit.log_resp", "res.log_resp"))
  if ("res.resp_normalized" %in% names(tmp.y)) setnames(tmp.y, "res.resp_normalized", "res.log_resp_normalized" )
    }
  # if (!is.null(out.csv)){
  #
  #   if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)
  #   write.csv(ptable, file =  file.path(out.csv, paste0("ptable ", "REML" ,".csv")), row.names = FALSE, na = "")
  #   write.csv(stable, file =  file.path(out.csv, paste0("stable ", "REML" ,".csv")), row.names = FALSE, na = "")
  #   write.csv(tmp.y, file = file.path(out.csv, paste0("prediction ", "REML" ,".csv")), row.names = FALSE, na = "")
  #   write.csv(aic.reml, file = file.path(out.csv, paste0("fitting ", "REML" ,".csv")), row.names = FALSE, na = "")
  # }

  return.lst <- list(model = m.sel, fitting = aic.reml, ptable = ptable, stable = stable, pred = tmp.y)
  if (m.option >= 3 & exists("moranI")) return.lst$moranI <- moranI
  if (length(m.candidates) > 1)  return.lst$fitting_ML <- aic.all

  return(return.lst)
}



#' calculate bai

#' @description
#' calculate basal area (cm2) and basal area increment (cm2)
#'
#'
#' @param dt.long data in long format containing at least 3 columns: id, year, rw
#' @param id column name of series id
#' @param rw column name of ring width measurement which is in mm
#'
#'
#' @import data.table
#'
#'
#' @return add 3 columns to the original input data, ageC for cambial age, ba_cm2_t_1 for basal area of the previous year in cm2, and bai_cm2 for annual basal area increment in cm2
#'
#' @export cal.bai


cal.bai <- function(dt.long, id , rw){
  setDT(dt.long)
  med.rw <- median(as.numeric(dt.long[[rw]]), na.rm = TRUE)
  # cat("please assure the unit of ", rw, " is mm")

  if (!all(c(id, "year", rw) %in% names(dt.long))) stop (paste0("at least one of the variables ", id, "year", rw, "not exists, please verify..."))

  if (nrow(dt.long[, .N, by = eval(c(id, "year"))][N>1]) > 0) stop (paste0(id, "-year is not a unique key, please verify..."))

  if (med.rw > 10) print(paste0("median of rw ", med.rw, " seems too big, assure it's in mm"))
  if (med.rw < 1) print(paste0("median of rw ", med.rw, " seems too small, assure it's in mm"))

  setorderv(dt.long, c(id, "year"))

  dt.long[, `:=`(ageC = seq_len(.N),
                 radius = cumsum(eval(parse(text = rw))),  # Cumulative radius (assumes RW is added each year)
                 radius_prev = shift(cumsum(eval(parse(text = rw))), fill = 0)), by = eval(id)]  # Previous radius (shifted)

  # Compute previous BA in cm2
  dt.long[, ba_cm2_t_1 := pi * (radius_prev^2)/100]

  dt.long[, bai_cm2 := pi * (radius^2 - radius_prev^2)/100]
  dt.long[bai_cm2 < 0]

  # Drop the radius columns if not needed
  dt.long[, c("radius", "radius_prev") := NULL]
  return(dt.long)
}


#' format wide to long

#' @description
#' convert the format of term prediction with by from wide to long in model$pred
#'
#'
#' @param model a gam model
#'
#' @import data.table
#'
#'
#' @return dt.pred, the term prediction by gam were formed as 1 term, for both fits and stand error
#'
#' @export format_byterm


format_byterm <- function(model, dt.pred){
  if (is.null(model) | is.null(dt.pred)) stop("please refit the model with package CFSTRenD...")

  setDT(dt.pred)
  rhs_terms <-attr(terms(model$formula), "term.labels")

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Keep only terms with 's(' and 'by ='
  smooth_with_by <- rhs_terms[str_detect(rhs_terms, "s\\([^\\)]+by\\s*=")]

  # Extract the variable and the 'by' category
  matched <- str_match(smooth_with_by, "s\\(([^,]+),[^)]*?by\\s*=\\s*([^,\\)\\s]+)")
  # matched <-setDT (as.data.frame(matched))
  # # Lists:
  # x.lst <- matched[, 2]  # first argument to s()
  # by.lst <- matched[, 3]    # value after 'by ='
  dt.pred$byterm <- NA_character_

  # not necessary as only by term appeared in matched
  # ix <- 1
  # while (ix <= nrow(matched)) {
  #
  #
  #   fit.s.clim <- paste0("fit.s.", matched[, 2][ix], "..", matched[, 3][ix])
  #   if (any(str_detect(names(dt.pred), fit.s.clim)) == TRUE) ix <- ix + 1 else matched <- matched[-ix,]
  #
  # }


  # se extraction from gratia is on the term only, so slightly smaller than mgcv "term" use mgcv to keep consistence.
  if (nrow(matched) > 0){

    for (ix in 1: nrow(matched)){

      fit.s.clim <- paste0("fit.s.", matched[, 2][ix], "..", matched[, 3][ix])

      for (se in c("", "se.")){

        cols_to_sum <- grep(paste0("^", se, fit.s.clim), names(dt.pred), value = TRUE)

        dt.pred$new_column <- rowSums(dt.pred[, cols_to_sum, with = FALSE])

        setnames(dt.pred, "new_column", paste0( se, strsplit(fit.s.clim, paste0(".", matched[, 3][ix]) )))
        dt.pred <- dt.pred[, !names(dt.pred) %in% cols_to_sum, with = FALSE]

      }
      dt.pred$byterm <- paste(paste0(matched[, 2], "-", matched[, 3]), collapse = ", ")

    }
    }

  return(dt.pred)
}
