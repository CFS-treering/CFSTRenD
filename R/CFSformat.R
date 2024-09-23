
#' convert data to CFS-TRenD format
#'
#' @param data : a list, first is input data in wide format; second is a flat sequence referring to the column indices of meta variables
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
  if (length(data) != 2) {
    print(length(data))
    stop("pls verify data is a list with 2 items, first is input data in wide format; second is a flat sequence referring to the column indices of meta variables")
  }
  if (!(is.data.frame(data[[1]])| is.data.table(data[[1]]))) stop("the first item of data must be the complete data in wide format")
  if (is.list(data[[2]] | is.null(data[[2]]) | class(data[[2]]) != "integer" )) stop("the second item of data must be a flat sequence referring to the column indices of meta variables")

    # {

    # if (is.null(N.cols_meta)) {
    #   stop("need to identify the columns of meta data")
    #   return()
    # }else{

      # setDT(data)
      names.data.o <- names(data[[1]])
      if (any(str_detect(names.data.o, "uid_radius.tmp"))) stop("please rename the column 'uid_radius.tmp'")

      cols_rw <- names.data.o[-data[[2]]]
      # check if all numbers
      # # is_numeric <- grepl("^[0-9]+$", x)
      # if ( all(grepl("^[0-9]+$", cols_rw))) {
      #   names(data[[1]])[-data[[2]]]<- paste0("Y", names(data[[1]])[-data[[2]]])
      #   } else{
      #
      #   # check if year with character prefix,
      #   if ( !(all(grepl("^[A-Za-z]+$", str_sub(cols_rw,1,1))) & all(grepl("^[0-9]+$", str_sub(cols_rw,2,-1))))) stop("colname of rw accepts 2 formats: year in numeric(NNNN) or year with 1-letter prefix (CNNNN)")
      #   }
      #

      # pattern <- "^[A-Za-z]+[0-9]+$"
     if (!(all(grepl("^[0-9]+$", cols_rw)) | all(grepl("^[A-Za-z]+[0-9]+$", cols_rw)))) stop("colname of rw accepts 2 formats: year in numeric(NNNN) or prefix-year (CNNNN)")
      idx.y5 <- as.numeric(str_extract(cols_rw,"\\(?[0-9,.]+\\)?")) > 9999
      if (any(idx.y5)) stop(paste0("please check the colname(s): ", paste(cols_rw[idx.y5], collapse = ", ")))

       # separate the data into 2

    # meta data
    dt.meta <- data.table(uid_radius.tmp = 1:nrow(data[[1]]),data[[1]][,data[[2]]])
    # tree ring
    dt.tr <- data.table(uid_radius.tmp = 1:nrow(data[[1]]),data[[1]][,-data[[2]]])

    dt.rwl <- melt(dt.tr, id.vars = c("uid_radius.tmp"))[!is.na(value)]
    setnames(dt.rwl, "value", "rw_mm")
    dt.rwl[, year:= as.numeric(str_extract(variable,"\\(?[0-9,.]+\\)?"))]

  # }
  # }
  # else
  # {
  # if (length(data) != 2) {
  #   print(length(data))
  #   stop("pls verify data only with 2 tables in the form list(tab1, tab2)")
  #   }
  #   if (nrow(data[[1]]) > nrow(data[[2]])) {
  #     dt.meta <- data[[2]]
  #     dt.rwl <- data[[1]]
  #   }else{
  #     dt.meta <- data[[1]]
  #     dt.rwl <- data[[2]]
  #   }
  #   if (!any(str_detect(names(dt.meta), "id_radius.tmp")) ) stop("please add column id_radius.tmp in metadata")
  #   if (!any(str_detect(names(dt.rwl), "id_radius.tmp")) ) stop("please add column id_radius.tmp in rw measurement")
  # # }

  meta.all0 <- variables.cfstrend()
  meta.all <- meta.all0[!str_detect(Variable, "uid_")]
  setdiff(meta.all0$Variable, meta.all$Variable)
  # must-have meta variables for tr_1 to tr_6
  cols.musthave <-as.character( meta.all[Required == 1 & !(table %in% c("tr_7_ring_widths", "tr_8_uids_removed"))]$Variable)
  add.Vars <- setdiff(cols.musthave, names(dt.meta))

if (length(add.Vars) > 0) {
  print("step 1: checking mandatory columns: please verify the the following columns.... ")
  message( paste(add.Vars, collapse = ", "))
  return()
}else{
  for (i in seq_along(cols.musthave)){
    if (nrow(dt.meta[is.na(get(cols.musthave[i]))]) > 0) message (paste0(cols.musthave[i]), " not complete, please check")
  }
  print("you have filled all the mandatory information")
  for (i in 1:6){
    i.table <- meta.all[str_detect(table, paste0("tr_", i)) & Required != 1 ]
    if (nrow(i.table) > 0) print(paste0(unique(i.table$table), ": ", paste0(i.table$Variable, collapse = ", ")))
    }
  # print("if you have filled all the information, go ahead to fill other variables")
}
  user_input <- readline(prompt = "do you have info of the above unmandatory columns? Y to return to complete info, or N to allow the system to fill with NA  (Y/N) : ")
  if (toupper(user_input) == "Y") return() else{
  fill.meta <- meta.all[!(Variable %in% names(dt.meta))]
  vars.fill <- as.character(fill.meta$Variable)
  FormatV.fill <- as.character(fill.meta$FormatV)

  lapply( seq_along(vars.fill), function(i) dt.meta[, c(vars.fill[i]) := eval(parse(text = FormatV.fill[i]))])
}
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

  ys <- dt.rwl[, .(rw_ystart = min(year), rw_yend = max(year)), by = c("uid_radius.tmp")]
  dt.new <- dt.new[ys, on = .(uid_radius.tmp)]
  dt.new[, c("ymin.proj", "ymax.proj"):= .(min(rw_ystart), max(rw_yend) ), by = .(uid_project)]
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

  # ylast <- tr_7_ring_widths[, .(ystart = min(year), ylast= max(year)), by = .(uid_radius)]
# wide format
  dt.new <- dt.new[, unique(meta.all0[table != "tr_7_ring_widths"]$Variable), with = FALSE]
  # dt.new <- merge(dt.new, ylast, by = "uid_radius")
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

  # return(list(tr_1_projects = tr_1_projects, tr_2_sites= tr_2_sites, tr_3_trees = tr_3_trees, tr_4_meas = tr_4_meas, tr_5_samples = tr_5_samples,
  #              tr_6_radiuses = tr_6_radiuses,  tr_7_ring_widths = tr_7_ring_widths, tr_all_wide = tr_all_wide, tr_meta = dt.new, complete_vars = complete_vars)
  #        )

  return(list(  tr_rw_long = tr_7_ring_widths, tr_all_wide = tr_all_wide, complete_vars = complete_vars)
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



#' plot the median of tree ring measurement of 1 site and of its nearest neighbors
#'
#' @param site2chk : data.table with columns uid_project, uid_site, site_id, species, longitude and latitude
#' @param tr_all_wide : reference sites including including rw measurement in wide format
#' @param ver.tr : release version of CFS-TRenD
#' @param N.nb : number of nearest-neighbors
#'
#' @import data.table
#' @import stringr
#' @import RANN
#' @import patchwork

#' @export CFS_scale
#'

#'

CFS_scale <- function(site2chk, tr_all_wide, ver.tr = "V1.2",  N.nb){
  # scale
  if (nrow(site2chk) > 1) stop("we can only process 1 species in 1 site each time ...")

  site.all.spc <- tr_all_wide[species == site2chk$species,c("uid_project", "uid_site", "site_id", "species", "latitude", "longitude")][, .N, by = .(uid_project, uid_site, site_id, species, latitude, longitude)]
  rw.all.spc <- merge(tr_all_wide, site.all.spc[, c("uid_project", "uid_site", "species")], by = c("uid_project", "uid_site", "species"))

  closest <- RANN::nn2(site.all.spc[,c("longitude", "latitude")], site2chk[ ,c("longitude", "latitude")], k = N.nb + 1)
  plot(1:length(closest$nn.dists), closest$nn.dists)
  uid_site.closest <- data.table(ord = 0:N.nb, uid_project = site.all.spc$uid_project[closest$nn.idx], uid_site = site.all.spc$uid_site[closest$nn.idx], nn.dists = as.numeric(closest$nn.dists))
  uid_site.closest <- merge(uid_site.closest, site.all.spc[, c("uid_project", "uid_site", "longitude", "latitude")], by = c("uid_project", "uid_site"))

  uid_radius.chk <- merge(rw.all.spc, uid_site.closest[, c("uid_project", "uid_site")], by = c("uid_project", "uid_site"))
  if (ncol(tr_all_wide) != ncol(uid_radius.chk)) stop("check the colname of uid_radius.chk")
  uid_radius.chk <- cbind(uid_radius.chk[, c("uid_project", "uid_site", "species", "uid_radius")], uid_radius.chk[, 42: ncol(uid_radius.chk)])
  uid_radius.chk.long <- melt(uid_radius.chk, id.vars = c( "species", "uid_project", "uid_site", "uid_radius"))[!is.na(value)]
  uid_radius.chk.long[, year:= as.numeric(as.character(variable))]
  med.radius <- uid_radius.chk.long[, .(N = .N, rw.median = median(value), yr.mn = min(year), yr.max = max(year)), by = .(species, uid_project, uid_site, uid_radius)]
  med.site <- med.radius[, .(Ncores = .N, rw.median = median(rw.median), yr.mn = min(yr.mn), yr.max = max(yr.max)), by = .(species, uid_project, uid_site)]

  med.radius.yr <- uid_radius.chk.long[, .(N = .N, rw.median = median(value) ), by = .(species, uid_project, uid_site, uid_radius, year)]
  med.site.yr <- med.radius.yr[, .(Ncores = .N, rw.median = median(rw.median)), by = .(species, uid_project, uid_site, year)]
  med.site.yr <- merge(med.site.yr, uid_site.closest, by = c("uid_project", "uid_site"))

  chk.site <- merge(uid_site.closest, med.site, by = c("uid_project", "uid_site"))
  setorder(med.site.yr, ord, uid_project, uid_site, year)



  g2 <- ggplot(med.site.yr, aes(x = year, y = rw.median, group = uid_site)) +
    geom_line(aes(color = factor(ifelse(ord == 0, "red", "darkblue"))), alpha = 0.6) +
    scale_color_manual(values = c("red" = "red", "darkblue" = "darkblue"),
                       labels = c("red" = "Site to Check", "darkblue" = "Neighbours"),
                       name = "") +
    theme(legend.position = "right",
          plot.caption = element_text(hjust = 0, face = "italic")) +  # Position the legend on the right

    labs(title = paste0( site2chk$site_id, " (",  site2chk$species, ")" ),

         y = "rw.median(mm)",
         caption = paste0("Data source: CFS-TRenD ", ver.tr)  # Add the data source caption
    )


  g1 <- ggplot(chk.site, aes(x = longitude, y = latitude, size = rw.median)) +
    geom_point(aes(color = ifelse(ord == 0, "red", "darkblue")), alpha = 0.6) +
    scale_color_identity() + #  Use the actual colalpha = 0.6, color = "darkblue") +
    # scale_shape_manual(values = c("circle" = 20)) + # Map "circle" to shape 16 and "cross" to shape 4
    scale_size_continuous(range = c(1, 20)) + #+  Adjust size range as needed
    # scale_x_continuous(breaks = sort(unique(data.tmp$lon)))


    theme_minimal() +
    theme(strip.text = element_text(size = 16),# Increase the size of facet labels
          panel.grid.minor = element_blank() , # Remove minor grid lines
          plot.title = element_text(size = 18), # Set title size
          plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt"), # Adjust plot margins
          plot.caption = element_text(hjust = 0, face = "italic")) + # Customize caption appearance

    labs(
      x = "Longitude",
      y = "Latitude",
      # caption = paste0("Data source: CFS-TRenD ", ver.tr) , # Add the data source caption
      size = paste0("rw median(mm)"))
  # labs()
  # +
  # geom_text(aes(label = ifelse(ord %in% c(0, 1), as.character(round(rw.median, 1)), "")),
  # hjust = 0.5, vjust = -1, check_overlap = TRUE)

  # gridExtra::grid.arrange(g2, g1, ncol = 2)
  g2 + g1 + plot_layout(widths = c(2, 1))
}
