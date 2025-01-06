
#' Convert tree-ring data into CFS-TRenD format
#'
#' @param data  a list, first is input data in wide format; second is a flat sequence referring to the column indices of ring measurement variables
#' @param out.csv  output csv file (default is NULL, otherwise to specify the directory to output)


#' @description
#' converts tree-ring data from various formats into a format compatible with hierarchical structure of the CFS-TRenD data collection (Girardin et al., 2021).
#'
#' @references Girardin, M.P., Guo, X.J., Metsaranta, J., Gervais, D., Campbell, E., Arsenault, A., Isaac-Rentone, M., Harvey, J.E., Bhatti, J., Hogg, E.A. 2021. A national tree-ring repository for Canadian forests (CFS-TRenD): structure, synthesis and applications. Environmental Reviews, 29 (999), 1-17. https://doi.org/10.1139/er-2020-0099

#'
#' @import data.table
#' @import stringr
#'
#' @return A list of 3 elements:
#' 1) A list containing seven tables compatible with CFS-TRenD data structure;
#' 2) A data table containing all the meta-data and ring width measurement in wide format;
#' 3) A data table for the percentage of completeness of each variable.

#'
#' @export CFS_format
#'

#'
CFS_format <- function (data, out.csv = NULL) {
  if (length(data) != 2) {
    print(length(data))
    stop("pls verify data is a list with 2 items, first is input data in wide format; second is a flat sequence referring to the column indices of ring measurement variables")
  }
  if (!(is.data.frame(data[[1]])| is.data.table(data[[1]]))) stop("the first item of data must be the complete data in wide format")
  if (is.list(data[[2]] | is.null(data[[2]]) | class(data[[2]]) != "integer" )) stop("the second item of data must be a flat sequence referring to the column indices of ring measurement variables")


      names.data.o <- names(data[[1]])
      # if (any(str_detect(names.data.o, "uid_radius.tmp"))) stop("please rename the column 'uid_radius.tmp'")

      # Check the columns starting with uid_, these are reserved column names

      uid_matches <- names.data.o[str_detect(names.data.o, str_c( c("uid_project", "uid_site", "uid_tree", "uid_meas", "uid_sample", "uid_radius", "uid_radius.tmp"), collapse = "|"))]

      # Print the detected uid_matches

      if (length(uid_matches) > 0) stop(paste0("please rename the column(s): ", paste(uid_matches, collapse = ", ")))

      cols_rw <- names.data.o[data[[2]]]

    # print(cols_rw)
      # pattern <- "^[A-Za-z]+[0-9]+$"
     if (!(all(grepl("^[0-9]+$", cols_rw)) | all(grepl("^[A-Za-z]+[0-9]+$", cols_rw)))) stop("colname of rw accepts 2 formats: year in numeric(NNNN) or prefix-year (CNNNN)")
      idx.y5 <- as.numeric(str_extract(cols_rw,"\\(?[0-9,.]+\\)?")) > 9999
      if (any(idx.y5)) stop(paste0("please check the colname(s): ", paste(cols_rw[idx.y5], collapse = ", ")))

       # separate the data into 2

    # meta data
    dt.meta <- data.table(uid_radius.tmp = 1:nrow(data[[1]]),data[[1]][,-data[[2]], with = FALSE])
    # tree ring
    dt.tr <- data.table(uid_radius.tmp = 1:nrow(data[[1]]),radius_id = data[[1]]$radius_id, data[[1]][,data[[2]], with = FALSE])

    dt.rwl <- melt(dt.tr, id.vars = c("uid_radius.tmp", "radius_id"))[!is.na(value)]
    setnames(dt.rwl, "value", "rw_mm")
    dt.rwl[, year:= as.numeric(str_extract(variable,"\\(?[0-9,.]+\\)?"))]


  meta.all0 <- variables.cfstrend()
  meta.all <- meta.all0[!str_detect(Variable, "uid_")]
  setdiff(meta.all0$Variable, meta.all$Variable)
  # must-have meta variables for tr_1 to tr_6
  cols.musthave <-as.character( meta.all[Required == 1 & !(table %in% c("tr_7_ring_widths", "tr_8_uids_updated"))]$Variable)
  add.Vars <- setdiff(cols.musthave, names(dt.meta))

if (length(add.Vars) > 0) {
  print("step 1: checking mandatory columns: please verify the the following columns.... ")
  message( paste(add.Vars, collapse = ", "))
  return()
}
  # checking variable completeness
  for (i in seq_along(cols.musthave)){
    if (nrow(dt.meta[is.na(get(cols.musthave[i]))]) > 0) message (paste0(cols.musthave[i]), " not complete, please check")
  }

  print("you have filled all the mandatory information")

  # for un-mandatory columns

  if (length(setdiff(meta.all[Required != 1]$Variable, names.data.o)) == 0) print("you have filled all the un-mandatory information") else{
  cat("\n the following are unmandatory information and not found in your data\n")

  # not required, and not reported, prompt to fill with NA

  i.table <- meta.all[Required != 1 & !(Variable %in% names.data.o)]
  i.table <- i.table[, .(var = paste(Variable, collapse = ",")), by = table]
  print (i.table)
  # for (i in 1:6){
  #   # ith table
  #   i.table <- meta.all[str_detect(table, paste0("tr_", i)) & Required != 1 & !(Variable %in% names.data.o) ]
  #   if (nrow(i.table) > 0) {
  #     print(paste0(unique(i.table$table), ": ", paste(i.table$Variable, collapse = ", ")))
  # }
  #   }
  # print("if you have filled all the information, go ahead to fill other variables")


    user_input <- readline(prompt = "do you have info of the above variables? Y to return to complete info, or N to allow the system to fill with NA  (Y/N) : ")
    if (toupper(user_input) == "Y") return() else{
    fill.meta <- meta.all[!(Variable %in% names(dt.meta))]
    vars.fill <- as.character(fill.meta$Variable)
    FormatV.fill <- as.character(fill.meta$FormatV)

    lapply( seq_along(vars.fill), function(i) dt.meta[, c(vars.fill[i]) := eval(parse(text = FormatV.fill[i]))])
  }
  }

  # trend variables

  dt.new <- copy(dt.meta)
  setDF(dt.new)
  dt.new[dt.new == ''] <- NA

  setDT(dt.new)
  dt.new[, uid_project := .GRP, by = .(project_name)]
  # chk site
  site_LL <- dt.new[, .N, by = .(site_id, latitude, longitude)]
  dup.LL <- site_LL[, .N, by = .(site_id)][N>1]
  dup.site <- site_LL[, .N, by = .(latitude, longitude)][N>1]
  if (nrow(dup.LL) > 0) stop(paste0("site_id: ", paste(dup.LL$site_id, collapse = ", "), " associated with multiple lat-lon. please verify..."))
  if (nrow(dup.site) > 0) {
    dup.site[, coord:= paste0("(", longitude, ", ", latitude, ")")]
    stop(paste0("lon-lat: ", paste(dup.site$coord, collapse = ", "), " associated with multiple site_id. please verify..."))
  }
# Canada's boundaries are roughly:
#
# Latitude range: 41.7°N to 83.1°N
# Longitude range: 52.6°W to 141.0°W
# source:  https://www12.statcan.gc.ca/census-recensement/2011/ref/dict/geo016-eng.cfm

  if (min(dt.new$latitude) < 41.7 | max(dt.new$latitude) > 83.1) stop(paste0("latitude not in range: (41.7 , 83.1)"))
  if (min(dt.new$longitude) < -141 | max(dt.new$longitude) > -52.6) stop(paste0("longitude not in range: (-141 , -52.6)"))

  dt.new[, uid_site := .GRP, by = .(site_id, latitude, longitude)]

  chk.spc<-dt.new[, .N, by = .(site_id, tree_id, species)][, .N, by = .(site_id, tree_id)][N > 1]

  if (nrow(chk.spc) > 0) {
    chk.spc[, site.tree:= paste0(site_id, "$", tree_id)]
    stop(paste0("site$tree: ", paste(chk.spc$site.tree, collapse = ", ")), " associated with multiple species, please verify...")
    }
  # chk if species is included in species_nficode
  # species_nficode is the species list of tree source, and was saved as internal data
  chk.spc2 <- setdiff(unique(dt.new$species), species_nficode$nfi_species_code)
  if (length(chk.spc2)  > 0) stop(paste0("species: ", paste(chk.spc2, collapse = ", "), " not recognized, please verify..."))

  dt.new[, uid_tree := .GRP, by = .(uid_project, uid_site, tree_id)]


  dt.new[, uid_meas := .GRP, by = .(uid_tree, meas_no)]
  dt.new[, uid_sample := .GRP, by = .(uid_meas, sample_id)]

  dt.new[, uid_radius := .GRP, by = .(uid_sample, radius_id)]

  ys <- dt.rwl[, .(rw_ystart = min(year), rw_yend = max(year)), by = c("uid_radius.tmp")]
  setorder(dt.rwl, uid_radius.tmp, year)
  dt.rwl[, ydif:= year - shift(year), by = "uid_radius.tmp"]
  dt.rwl[, rwinc:= rw_mm - shift(rw_mm), by = "uid_radius.tmp"]
  chk.ydif <- dt.rwl[, .SD[-1], by = uid_radius.tmp][is.na(ydif)][, .N, by = .(radius_id)]
  if (nrow(chk.ydif) > 0) stop(paste0("radius_id: ", paste(chk.ydif$radius_id, collapse = ", "), " with missing year, please verify..." ))
  # str(dt.rwl)
  chk.rw <- dt.rwl[ rw_mm < 0][, .N, by = .(radius_id)]
  if (nrow(chk.rw) > 0) stop(paste0("radius_id: ", paste(chk.ydif$radius_id, collapse = ", "), " with negative rw measurement, please verify..." ))

  dt.new <- dt.new[ys, on = .(uid_radius.tmp)]
  dt.new[, c("ymin.proj", "ymax.proj"):= .(min(rw_ystart), max(rw_yend) ), by = .(uid_project)]
  dt.new[, yr.meas:= max(ys), by = .(uid_meas)]

  dt.new[is.na(meas_date), meas_date:= paste0(yr.meas, "-00-00")]
  dt.new[, year_range:= paste0(ymin.proj, " ; ", ymax.proj)]
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
  dt.new <- dt.new[, unique(meta.all0[!(table %in% c("tr_7_ring_widths", "tr_8_uids_updated"))]$Variable), with = FALSE]
  # dt.new <- merge(dt.new, ylast, by = "uid_radius")
  tr_all_wide <- dcast(tr_7_ring_widths, uid_radius ~ year, value.var = "rw_mm")
  tr_all_wide <- merge(dt.new, tr_all_wide, by = "uid_radius")
  if (!is.null(out.csv)){

      if (!(dir.exists(out.csv))) dir.create(out.csv, recursive = TRUE)

      for (i.tbl in 1:7){
        write.csv(eval(parse(text = fn7[i.tbl])), file =file.path(out.csv, paste0(fn7[i.tbl], ".csv" )), na = "",  row.names = FALSE)

      }
      write.csv(tr_all_wide, file =file.path(out.csv, paste0("tr_all_wide", ".csv" )), na = "",  row.names = FALSE)
      write.csv(complete_vars, file =file.path(out.csv, paste0("complete_vars", ".csv" )), na = "",  row.names = FALSE)


      }


  return(list(  tr_all_long = list(tr_1_projects = tr_1_projects, tr_2_sites= tr_2_sites, tr_3_trees = tr_3_trees, tr_4_meas = tr_4_meas, tr_5_samples = tr_5_samples,
                                                tr_6_radiuses = tr_6_radiuses,  tr_7_ring_widths = tr_7_ring_widths), tr_all_wide = tr_all_wide, complete_vars = complete_vars)
  )

}


#' compare the median of tree ring measurement of a specific site to those of nearby sites
#' @description
#' Apply a k-nearest neighbors (k-NN) method based on geographic location for the same species.
#' It compares the median tree-ring measurements of a specific site to those of nearby sites
#'
#' @param site2chk  data.table with columns uid_project, uid_site, site_id, species, longitude and latitude
#' @param ref_sites  reference sites including ring width measurements in wide format
#' @param N.nbs  number of nearest-neighbors
#' @param make.plot  a list , first item for plot figures or not, second item for the caption of data source
#'
#' @import data.table
#' @import stringr
#' @import RANN
#' @importFrom scales pretty_breaks
#' @import patchwork
#'
#' @return A data table containing the median ring-width measurements of the involved sites, along with the distances from the specific site

#' @export CFS_scale
#'

#'

CFS_scale <- function(site2chk, ref_sites, N.nbs = 10, make.plot = c(FALSE, "")){
  # scale
  if (nrow(site2chk) > 1) stop("we can only process 1 species in 1 site each time ...")

  site.all.spc <- ref_sites[species == site2chk$species,c("uid_project", "uid_site", "site_id", "species", "latitude", "longitude")][, .N, by = .(uid_project, uid_site, site_id, species, latitude, longitude)]
  rw.all.spc <- merge(ref_sites, site.all.spc[, c("uid_project", "uid_site", "species")], by = c("uid_project", "uid_site", "species"))

  closest <- RANN::nn2(site.all.spc[,c("longitude", "latitude")], site2chk[ ,c("longitude", "latitude")], k = N.nbs + 1)
  # plot(1:length(closest$nn.dists), closest$nn.dists)
  uid_site.closest <- data.table(ord = 0:N.nbs, uid_project = site.all.spc$uid_project[closest$nn.idx], uid_site = site.all.spc$uid_site[closest$nn.idx], nn.dists = as.numeric(closest$nn.dists))
  uid_site.closest <- merge(uid_site.closest, site.all.spc[, c("uid_project", "uid_site", "longitude", "latitude")], by = c("uid_project", "uid_site"))

  uid_radius.chk <- merge(rw.all.spc, uid_site.closest[, c("uid_project", "uid_site")], by = c("uid_project", "uid_site"))
  if (ncol(ref_sites) != ncol(uid_radius.chk)) stop("check the colname of uid_radius.chk")
  uid_radius.chk <- cbind(uid_radius.chk[, c("uid_project", "uid_site", "species", "uid_radius")], uid_radius.chk[, 45: ncol(uid_radius.chk)])
  # only consider the value > 0
  uid_radius.chk.long <- melt(uid_radius.chk, id.vars = c( "species", "uid_project", "uid_site", "uid_radius"))[!is.na(value)][value > 0]
  uid_radius.chk.long[, year:= as.numeric(as.character(variable))]
  med.radius <- uid_radius.chk.long[, .(N = .N, rw.median = median(value), yr.mn = min(year), yr.max = max(year)), by = .(species, uid_project, uid_site, uid_radius)]
  med.site <- med.radius[, .(Ncores = .N, rw.median = median(rw.median), yr.mn = min(yr.mn), yr.max = max(yr.max)), by = .(species, uid_project, uid_site)]

  med.radius.yr <- uid_radius.chk.long[, .(N = .N, rw.median = median(value) ), by = .(species, uid_project, uid_site, uid_radius, year)]
  med.site.yr <- med.radius.yr[, .(Ncores = .N, rw.median = median(rw.median)), by = .(species, uid_project, uid_site, year)]
  med.site.yr <- merge(med.site.yr, uid_site.closest, by = c("uid_project", "uid_site"))

  chk.site <- merge(uid_site.closest, med.site, by = c("uid_project", "uid_site"))

  # ratio defines as site2chk/neighbor , when it's smaller means the site2chk's magnitude is smaller,
  chk.site$rw.ratio <- chk.site[ord == 0]$rw.median/chk.site$rw.median

  chk.site$size_class <- cut(
    chk.site$rw.ratio,
    breaks = c(-Inf, 0.1, 0.2, 0.5, 2, 5, 10, Inf), # Define the breaks
    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.5", "0.5 - 2", "2 - 5", "5 - 10", "> 10"), # Define labels
    include.lowest = TRUE
  )
  if (any(is.na(chk.site$size_class))) {
    warning("Some values in rw.ratio are outside the defined breaks and will be excluded.")
  }


  setorder(med.site.yr, ord, uid_project, uid_site, year)

  rw.median <- data.table(chk.site[ord == 0][,c("uid_project", "uid_site","species", "rw.median")],
                      chk.site[ord > 0][, .(N.nbs = .N, rw.median.nbs = median(rw.median))])[,ratio_median := round(rw.median/rw.median.nbs,1)]


  if (make.plot[1] == TRUE){
  g2 <- ggplot(med.site.yr, aes(x = year, y = rw.median, group = uid_site)) +
    geom_line(aes(color = factor(ifelse(ord == 0, "red", "darkblue"))), alpha = 0.6) +
    scale_color_manual(values = c("red" = "red", "darkblue" = "darkblue"),
                       labels = c("red" = "Site", "darkblue" = "Nbs"),
                       name = "") +
    theme(legend.position = "right",
          plot.caption = element_text(hjust = 0, face = "italic")) +  # Position the legend on the right

    labs(
            y = "rw.median (mm)"

    )


  g1 <- ggplot(chk.site, aes(x = longitude, y = latitude, size = size_class)) +
    geom_point(aes(color = ifelse(ord == 0, "red", "darkblue")), alpha = 1) +
    scale_color_identity() + # Use the actual colors
    scale_size_manual(
      values = c("< 0.1" = 1, "0.1 - 0.2" = 2, "0.2 - 0.5" = 3, "0.5 - 2" = 5, "2 - 5" = 7, "5 - 10" = 9, "> 10" = 11), # Map size classes to point sizes
      name = "ratio(./nb)"
    ) +
    scale_x_continuous(breaks = pretty_breaks(n = 3), expand = expansion(mult = 0.3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5), expand = expansion(mult = 0.2)) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 16), # Increase the size of facet labels
      panel.grid.minor = element_blank(), # Remove minor grid lines
      plot.title = element_text(size = 18) # Set title size
    ) +
    labs(
      x = "Longitude",
      y = "Latitude"
    )


  print(g2 + g1 + plot_layout(widths = c(1.2, 1)) +
          plot_annotation(
            title = paste0( "uid_site ", site2chk$uid_site, " and its ", N.nbs, " neighbors ( species: ", site2chk$species, ")"  ),
            caption = paste0("Data source: ", make.plot[2]),
            tag_levels = 'a',
            tag_suffix = ")") &
          theme(
            plot.title = element_text(face = "bold"),
            plot.tag = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0.2, vjust = -1, face = "italic" )

            # plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt") # Adjust plot margins
            ))
  }
  setnames(chk.site, "uid_site", "uid_site.nb")
  chk.site<- data.table(rw.median[, "uid_site"], chk.site)
  return(list(chk.site = chk.site, ratio.median = rw.median))
}


#' frequency distributions by geo-location per species
#'
#' @param tr_meta   meta table from function CFS_format()
#' @param uid_level which uid level to count(uid_project, uid_site, uid_tree, uid_meas, uid_sample, uid_radius)
#' @param cutoff_yr cut-off year for a subset which series were recorded on or after
#' @param geo_resolution resolution of longitude and latitude in degree, default: c(5,5)


#' @import data.table
#' @import stringr
#' @import ggplot2
#'
#' @return a data table of counts of uid by latitude-longitude per species
#' @export CFS_freq
#'

CFS_freq <- function(tr_meta,  uid_level, cutoff_yr = -999,geo_resolution = c(5,5)){
  if (!(uid_level %in% c("uid_project", "uid_site", "uid_tree", "uid_meas", "uid_sample", "uid_radius"))) stop("uid_level should be in c('uid_project', 'uid_site', 'uid_tree', 'uid_meas', 'uid_sample', 'uid_radius')")


  dt.meta.sel <- tr_meta[rw_yend >= cutoff_yr, c(unique(c("uid_radius", uid_level)), "longitude", "latitude",  "species"), with = FALSE]
  dt.meta.sel[,lon:=round(longitude/geo_resolution[1], 0)*geo_resolution[1]]
  dt.meta.sel[, lat:= round(latitude/geo_resolution[2], 0)*geo_resolution[2]]

  uids.sll <- dt.meta.sel[, .N, by = c("lat", "lon", "species", uid_level)][,N:= NULL]
  uids.spc <- uids.sll[,.N, by = .(species)]
  setorder(uids.spc, -N)
  uids.spc[, pct.species := N/sum(N)]
  uids.spc[, pct.species := round(pct.species * 100,0)]
  uids.spc[, ord:= .I]
  uids.sll <- uids.sll[uids.spc[,c("species", "pct.species", "N", "ord")], on = .(species)]

  dist_uids <- uids.sll[, .(nuids = .N), by = .(ord, species, pct.species, N, lat, lon)]


  dist_uids <- dcast(dist_uids, ord + species + N + pct.species + lat ~ lon, value.var = "nuids")
  if (cutoff_yr > 0) dist_uids <- data.table(uid_label = paste0(uid_level, "_yr", cutoff_yr), dist_uids) else
    dist_uids <- data.table(uid_label = uid_level, dist_uids)
  setorder(dist_uids, -pct.species, species, -lat)


  return(dist_uids)
}


