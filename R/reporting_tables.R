



#' Generate table for species report
#'
#' This function generates the data table for data reporting at project-species level .
#' @param tr_6i data table of meta data for a specific species
#' @param dt.spc_radii data table of tree ring data for a specific species

#'
#' @return a data table.
#' @export

#'
table_spc <- function(tr_6i, dt.spc_radii){

  dt.radii.spc <- merge(tr_6i, dt.spc_radii, by = "uid_radius")
  # dt.radii.spc <- dt.radii[species == spc.lst[i.spc]]
  dt.summary.0 <- tr_6i[, .(lat = paste0(round(min(latitude),1), " ", round(max(latitude),1)),
                            lon = paste0(round(min(longitude),1), " ", round(max(longitude),1))), by = .(species)]

  setnames(dt.summary.0, c("lat", "lon"), c("Latitude range", "Longitude range"))
  dt.summary.0a <- melt(dt.summary.0, measure.vars = names(dt.summary.0), variable.name = "Category", value.name = "Value")





  dt.summary.1 <- dt.radii.spc[, .(Yspan = paste0(min(ymin), " ", max(ymax)), Nsites = length(unique(uid_site)), Ntrees = length(unique(uid_tree)),
                                   Nsamples = length(unique(uid_sample)), N = .N,
                                   summ_rw = paste0(round(mean(rw.mean), 1), " \u00B1 ", round(sd(rw.mean), 1), " ( ", round(min(rw.min), 1), ", ", round(max(rw.max), 1), " )"),
                                   summ_len = paste0(round(mean(N), 1), " \u00B1 ", round(sd(N), 1), " ( ", round(min(N), 1), ", ", round(max(N), 1), " )"),
                                   summ_ar1 = paste0(round(mean(ar1_rw), 1), " \u00B1 ", round(sd(ar1_rw), 1), " ( ", round(min(ar1_rw), 1), ", ", round(max(ar1_rw), 1), " )"))]

  dt.summary.1[, (names(dt.summary.1)) := lapply(.SD, as.character)]
  setnames(dt.summary.1, c("Yspan","Nsites", "Ntrees", "Nsamples", "N", "summ_rw", "summ_len", "summ_ar1"), c("Year span", "Number of sites", "Number of trees", "Number of samples", "Number of series", "summary rw(mm)**", "summary len.series**", "summary AR1**"  ))

  dt.summary.1a <- melt(dt.summary.1, measure.vars = names(dt.summary.1), variable.name = "Category", value.name = "Value")

  dt.summary.2a <- rbind(dt.summary.0a, dt.summary.1a)


  if ("qa_code" %in% names(dt.radii.spc))  {
    dt.Npass<- data.table(Category = "Number of series (pass)*", Value = nrow(dt.radii.spc[qa_code == "pass"]), ord = 8.5 )} else{
      dt.Npass<- data.table(Category = "Number of series (pass)*", Value = "Not Applicable", ord = 8.5 )
    }
  dt.summary.2a[, ord:=.I]
  summary_spc <- rbind(dt.summary.2a, dt.Npass)
  setorder(summary_spc, ord)
  summary_spc[,ord:=NULL]
  # dt.summary.3a$species <- unique(tr_6ispc$species)
  return(summary_spc )
}


# #' Generate table for radii report
# #'
# #' This function generates the data table for data reporting at project-species-site-radii level .
# #' @param tr_6i data table of meta data for a specific species
# #' @param tr_7i data table of tree ring data for a specific species
# #' @param rep_radii logic, for user to choose if the radii table is presented in the data report.
# #' @param ... additional arguments
# #' @param ... Additional arguments passed to \code{CFS_qa()}.
# #' @return a data table.
# #' @keywords internal
# #' @noRd
# #'
# table_spc_site_radii <- function(tr_6i, tr_7i, rep_radii = TRUE, ...){
#   args <- list(...)
#   qa.label_data <- args$qa.label_data
#   qa.min_nseries <- args$qa.min_nseries
# # if (rep_radii == TRUE & nrow(tr_6i) >= qa.min_nseries) {
#   # note thant nrow(tr_6i) checking in template_data_report.rmd 2025-08-16
#   if (rep_radii == TRUE) {
#   tr_7ispc<- copy(tr_7i)
#   # setnames(tr_7ispc, c("uid_radius", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))
#   tr_7ispc[, c("SampleID", "Year" ,"RawRing"):= .(uid_radius, year, rw_mm)]
#   # check duplication
#  if (nrow( tr_7ispc[, .N, by = .(SampleID, Year)][N>1]) > 0) stop("SampleID-Year is not unique key, please check...")
#   # include RW_trt, this is the series that qa process works on,
#   # tr_7ispc[, RW_trt:= RawRing - shift(RawRing), by = SampleID]; label_trt <- "differentiated"
#   if (!("RW_trt" %in% names(tr_7ispc))) stop("please assign the column RW_trt for running the qa process...")
#   acf.trt <- tr_7ispc[!is.na(RW_trt), .( ar1_trt = round(acf(RW_trt, plot = FALSE)$acf[2],2) ), by = .(SampleID)]
#
#   # run quality assessment procedure
#   qa.trt.ccf <-CFS_qa(dt.input = tr_7ispc, qa.label_data = qa.label_data, qa.min_nseries = qa.min_nseries)
#   dt.radii.spc <- copy(qa.trt.ccf$dt.stats)
#
#   dt.radii.spc[, uid_radius:= as.integer(SampleID)]
#   dt.radii.spc <- merge(tr_6i[, c("uid_radius", "site_id", "radius_id")], dt.radii.spc, by = "uid_radius")
#
# }else{
#   # qa.trt.ccf <- list()
#   dt.radii.spc <- tr_7i[, .(N = .N, rw.mean = mean(rw_mm), rw.sd = sd(rw_mm), rw.min = min(rw_mm), rw.max = max(rw_mm), ymin = min(year), ymax = max(year), ar1_rw = acf(rw_mm, plot = FALSE)$acf[2] ), by = .(site_id, radius_id, uid_radius)]
#   # if more than 10 series, running the correlation with mean
#   if (nrow(tr_6i) > 10){
#     # mean.rw <- tr_7i[, mean_rw:= mean(rw_mm), by = uid_radius]
#     tmp <- merge(tr_7i, tr_7i[, .(mean_rw= mean(rw_mm)), by = year], by = "year")
#     chk <- tmp[, .(ncorr_mean_rw = .N, result = cor.test(mean_rw, rw_mm)), by = .(uid_radius)]
#
#     chk.corr <- chk[, .SD[4], by = uid_radius][, corr_mean_rw:= as.numeric(unlist(result))]
#
#     chk.pvalue <- chk[, .SD[3], by = uid_radius][, pcorr_mean_rw:= as.numeric(unlist(result))]
#
#     dt.radii.spc <- dt.radii.spc[chk.corr[, c("uid_radius", "corr_mean_rw", "ncorr_mean_rw")], on = "uid_radius"][chk.pvalue[, c("uid_radius", "pcorr_mean_rw")], on = "uid_radius"]
#   }
# }
#   return(dt.radii.spc)
# }


#' Generate table for site report including the kNN result
#'
#' This function generates the data table for data reporting at project-species-site level .
#' @param dt.radii.spc result table from table_spc_site_radii
#' @param rw_ref data table of reference for kNN, in wide format
#' @param scale.N_nbs number of neighbor to search

#'
#' @return a data table.
#' @export


table_spc_site <- function(dt.radii.spc, rw_ref,  scale.N_nbs){
  dt.site_radii <-merge(rw_ref[,c("uid_site", "site_id", "species", "latitude", "longitude", "uid_tree", "uid_sample", "uid_radius")],
                        dt.radii.spc[, c("uid_radius", "N", "rw.mean", "rw.min", "rw.max")], by = "uid_radius")
  # dt.site.stats <- dt.site[, .(Ntrees = length(unique(uid_tree)), Nsamples = length(unique(uid_sample)), Nradius = .N), by = .(uid_site, site_id, longitude, latitude, species)]


  dt.site.stats <- dt.site_radii[ ,.(
    Ntrees = length(unique(uid_tree)), Nsamples = length(unique(uid_sample)), Nradius = .N,
    summ_len = paste0(round(mean(N), 1), " \u00B1 ", round(sd(N), 1), " ( ", round(min(N), 1), ", ", round(max(N), 1), " )"),
    summ_rw = paste0(round(mean(rw.mean), 1), " \u00B1 ", round(sd(rw.mean), 1), " ( ", round(min(rw.min), 1), ", ", round(max(rw.max), 1), " )")),
    by = .(uid_site, site_id, longitude, latitude, species)]
  # setnames(dt.site.stats, c("summ_rw", "summ_len"), c("summary rw(mm)**", "summary len.series**"))

  # head(dt.site.stats)
  if (length(unique(dt.site.stats$uid_site)) < scale.N_nbs + 2) {
    message("not sufficiant reference sites")
    dt.site.out <- dt.site.stats
  } else{


    list.scale <- CFS_scale(target_site = dt.site.stats[ ,c("species", "site_id")], ref_sites = rw_ref, scale.N_nbs = scale.N_nbs)

    dt.scale.all <- data.table::rbindlist(
      lapply(list.scale, `[[`, "ratio.median"),
      use.names = TRUE, fill = TRUE
    )
    # summary table of spc

    dt.site.out <- merge(dt.site.stats, dt.scale.all[, c("uid_site", "rw.median", "rw.median.nbs", "ratio_median")], all.x = TRUE,  by = "uid_site")
  }
  return(dt.site.out)
}

