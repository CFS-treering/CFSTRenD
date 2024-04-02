
#' convert data to CFS-TRenD format
#'
#' @param data : input data in wide format
#' @param out.csv : output csv file (TRUE or FALSE)
#' @param out.dir : directory of output csv files

#'
#' @import data.table
#' @import stringr
#'
#' @return 7 tables starting with tr_
#' @export CFS_format
#'

#'
CFS_format <- function (data, out.csv = FALSE, out.dir = NULL) {

  # if (is.data.frame(data)| is.data.table(data)){
  #
  #   if (is.null(N.cols_meta)) {
  #     stop("need to identify the columns of meta data")
  #     return()
  #   }else{
  #
  #     setDT(data)
  #     names.data <- names(data)
  #
  #     cols_rw <- names.data[-N.cols_meta]
  #
  #     data[,uid_radius.tmp:= .I]
  #   # separate the data into 2
  #
  #     # meta data
  #     dt.meta <- data[, c("uid_radius.tmp",-cols_rw), with = FALSE]
  #
  #   # tree ring
  #   # summary(data$uid_radius.tmp)
  #   dt.tr <- data[, c("uid_radius.tmp",cols_rw), with = FALSE]
  #
  #   dt.rwl <- melt(dt.tr, id.vars = c("uid_radius.tmp"))[!is.na(value)]
  #   setnames(dt.rwl, "value", "rw_mm")
  #   dt.rwl[, year:= as.numeric(str_extract(variable,"\\(?[0-9,.]+\\)?"))]
  #
  # }
  # }else
  # {
    if (length(data) != 2) {
      print(length(data))
      stop("pls verify data only with 2 tables in the form list(tab1, tab2)")
      }
    if (nrow(data[[1]]) > nrow(data[[2]])) {
      dt.meta <- data[[2]]
      dt.rwl <- data[[1]]
    }else{
      dt.meta <- data[[1]]
      dt.rwl <- data[[2]]
    }
    if (!any(str_detect(names(dt.meta), "id_radius.tmp")) ) stop("please add column id_radius.tmp in metadata")
    if (!any(str_detect(names(dt.rwl), "id_radius.tmp")) ) stop("please add column id_radius.tmp in rw measurement")
  # }

  meta.all0 <- variables.cfstrend()
  meta.all <- meta.all0[!str_detect(Variable, "uid_")]
  setdiff(meta.all0$Variable, meta.all$Variable)
  cols.musthave <-as.character( meta.all[Required == 1 & table != "tr_7_ring_widths"]$Variable)
  add.Vars <- setdiff(cols.musthave, names(dt.meta))

if (length(add.Vars) > 0) {
  print("step 1: checking required columns: pls add the following columns: ")
  message( paste(add.Vars, collapse = ", "))
  return()
}else{
  for (i in seq_along(cols.musthave)){
    if (nrow(dt.meta[is.na(get(cols.musthave[i]))]) > 0) message (paste0(cols.musthave[i]), " not complete, please check")
  }
  print("all mandatory variables is in table, go ahead to fill other variables")
}

  fill.meta <- meta.all[!(Variable %in% names(dt.meta))]
  vars.fill <- as.character(fill.meta$Variable)
  FormatV.fill <- as.character(fill.meta$FormatV)

  lapply( seq_along(vars.fill), function(i) dt.meta[, c(vars.fill[i]) := eval(parse(text = FormatV.fill[i]))])

  # trend variables

  dt.new <- dt.meta[, c("uid_radius.tmp", meta.all$Variable), with = FALSE]
  setDF(dt.new)
  dt.new[dt.new == ''] <- NA

  setDT(dt.new)
  dt.new[, uid_project := .GRP, by = .(project_name)]

  dt.new[, uid_site := .GRP, by = .(site_id, latitude, longitude)]


  dt.new[, uid_tree := .GRP, by = .(uid_project, uid_site, tree_id)]
  dt.new[, uid_meas := .GRP, by = .(uid_tree, meas_no)]
  dt.new[, uid_sample := .GRP, by = .(uid_meas, sample_id)]

  dt.new[, uid_radius := .GRP, by = .(uid_sample, radius_id)]

  ys <- dt.rwl[, .(ymin = min(year), ymax = max(year)), by = c("uid_radius.tmp")]
  dt.new <- dt.new[ys, on = .(uid_radius.tmp)]
  dt.new[, c("ymin.proj", "ymax.proj"):= .(min(ymin), max(ymax) ), by = .(uid_project)]
  dt.new[, yr.meas:= max(year), by = .(uid_meas)]
  dt.new[, year_range:= paste0(ymin.proj, " ; ", ymax.proj)]
  dt.new[is.na(meas_date), meas_date:= paste0(yr.meas, "-00-00")]
  #

  # in 7 tables
  fn7 <- sort(unique(meta.all0$table))

  for (i.tbl in 1:6){
    ord <- unlist(meta.all0[table == fn7[i.tbl]]$Variable)
    tmp <- dt.new[, .N, by = eval(ord)][, N:=NULL]
     assign(fn7[i.tbl], tmp)
    # check completeness
     if (i.tbl == 4) tmp[meas_no == "-999" | meas_no == -999, meas_no:=NA]
     if (i.tbl == 6) tmp[bark_thickness_mm == 999, bark_thickness_mm:= NA_real_]
    setDF(tmp)
     # tmp<-tmp[,!grepl( "uid_" , names( tmp ) )]
    tmp[tmp == ''] <- NA

    v_count <-sapply(tmp, function(y) sum(length(which(!is.na(y)))))
    v_count <- round(v_count / dim(tmp)[1] * 100, 1)

    v_count <- data.table(a = matrix(v_count, ncol = 1))
    names(v_count) <- "pct"
    v_count$var <- paste0("tr", i.tbl, "_", ord)
    # remove uids
    v_count <- v_count[!str_detect(var, "uid_") ]
    setcolorder(v_count, c("var", "pct"))
      if (i.tbl == 1)
        complete_vars <- v_count else
          complete_vars <- rbind(complete_vars, v_count)


  }



  tr_7_ring_widths <- merge(dt.new[, c("uid_radius.tmp", "uid_radius")], dt.rwl, by = "uid_radius.tmp")[, c("uid_radius","year" , "rw_mm")]
  setorder(tr_7_ring_widths, uid_radius,year)
  tr_7_ring_widths <- tr_7_ring_widths[, as.character(meta.all0[table == fn7[7]]$Variable), with = FALSE]
# wide format
  dt.new <- dt.new[, unique(meta.all0[table != "tr_7_ring_widths"]$Variable), with = FALSE]
  tr_all_wide <- dcast(tr_7_ring_widths, uid_radius ~ year, value.var = "rw_mm")
  tr_all_wide <- merge(dt.new, tr_all_wide, by = "uid_radius")
  if (out.csv == TRUE){
    if (is.null(out.dir)) message("please provide the directiory to output csv files")  else{
      if (!(dir.exists(out.dir))) dir.create(out.dir, recursive = TRUE)
      for (i.tbl in 1:7){
        write.csv(eval(parse(text = fn7[i.tbl])), file =file.path(out.dir, paste0(fn7[i.tbl], ".csv" )), na = "",  row.names = FALSE)

      }
      write.csv(complete_vars, file =file.path(out.dir, paste0("complete_vars", ".csv" )), na = "",  row.names = FALSE)
      save(tr_all_wide, file =paste0("tr_all_wide", ".Rdata" ))
      }
    }

  return(list(tr_1_projects = tr_1_projects, tr_2_sites= tr_2_sites, tr_3_trees = tr_3_trees, tr_4_meas = tr_4_meas, tr_5_samples = tr_5_samples,
               tr_6_radiuses = tr_6_radiuses,  tr_7_ring_widths = tr_7_ring_widths, tr_all_wide = tr_all_wide, complete_vars = complete_vars)
         )

}

#' seperate 1 table into 2: dt.meta and dt.rwl with the link id_radius.tmp
#'
#' @param data : dataframe or datatable in wide format
#' @param N.cols_meta : all column index except ring measurement
#' @import data.table
#' @import stringr
#' @return 2 tables: dt.meta and dt.rwl
#' @export dt.sep2
#'

dt.sep2 <- function(data, N.cols_meta=NULL){

  if (is.null(N.cols_meta)) stop("pls specify the column index of meta data")
    setDT(data)
  names.data <- names(data)
  cols_rw <- names.data[-N.cols_meta]
  if (any(names.data[N.cols_meta] == "uid_radius")) uidr.tmp <- "uid_radius" else {uidr.tmp <-"uid_radius.tmp"
   data[,uid_radius.tmp:= .I]}
  # separate the data into 2

    # meta data
  if (any(names.data[N.cols_meta] == "uid_radius")) dt.meta <- data[, c(names.data[N.cols_meta]), with = FALSE] else
    dt.meta <- data[, c(uidr.tmp,names.data[N.cols_meta]), with = FALSE]

  # tree ring
  # summary(data$uid_radius.tmp)
  dt.tr <- data[, c(uidr.tmp,cols_rw), with = FALSE]

  dt.rwl <- melt(dt.tr, id.vars = c(uidr.tmp))[!is.na(value)]
  setnames(dt.rwl, "value", "rw_mm")
  dt.rwl[, year:= as.numeric(str_extract(variable,"\\(?[0-9,.]+\\)?"))]
  dt.rwl <- dt.rwl[ ,c(uidr.tmp, "year", "rw_mm"), with = FALSE]
  return(list(dt.meta = dt.meta, dt.rwl = dt.rwl))
  }
