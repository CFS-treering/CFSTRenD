#' data spatial distribution by species measured later than certain year
#'
#' @param tr_meta :  meta table from CFS_format
#' @param crit_yr : year criteria
#' @param resol_lat :latitude resolution
#' @param resol_lon :longitude resolution
#' @param uid : counts on which uid
#' @param spc.lst species list
#' @param N.species : number of top species to be output (999 for all)
#' @param make.plot plot the distribution (TRUE or FALSE)
#' @param out.csv : output csv file (TRUE or FALSE)
#' @param out.dir : directory of output csv files
#'
#' @import data.table
#' @import stringr
#' @import ggplot2
#'
#' @return 1 table species-lat * lon, values are number of uids
#' @export CFS_freq
#'

CFS_freq <- function(tr_meta, crit_yr,resol_lat, resol_lon,  uid, spc.lst, N.species = 999, make.plot = TRUE, out.csv = FALSE, out.dir = NULL){
  if (!(uid %in% c("uid_project", "uid_site", "uid_tree", "uid_meas", "uid_sample", "uid_radius"))) stop("uid should be in c('uid_project', 'uid_site', 'uid_tree', 'uid_meas', 'uid_sample', 'uid_radius')")
  # if (length(data) == 1){
  #
  #   # separate the data into 2
  #
  #   # tree ring
  #   # summary(data$uid_radius.tmp)
  #   dt.tr <- data[, c("uid_radius",names(data)[42:ncol(data)]), with = FALSE]
  #
  #   dt.rwl <- melt(dt.tr, id.vars = c("uid_radius"))[!is.na(value)]
  #   setnames(dt.rwl, "value", "rw_mm")
  #   dt.rwl[, year:= as.numeric(str_extract(variable,"\\(?[0-9,.]+\\)?"))]
  #
  #   # meta data
  #   dt.meta <- data[, c(names(data)[1:41]), with = FALSE]
  #
  # }
  # if (length(data)  == 2) {
  # if (length(data) != 2) {
  #   print(length(data))
  #   stop("pls verify data only with 2 tables in the form list(tab1, tab2)")
  # }
  # if (nrow(data[[1]]) > nrow(data[[2]])) {
  #   dt.meta <- data[[2]]
  #   dt.rwl <- data[[1]]
  # }else{
  #   dt.meta <- data[[1]]
  #   dt.rwl <- data[[2]]
  # }
  # }

# ylast <- dt.rwl[, .(ylast= max(year)), by = .(uid_radius)][ylast >= crit_yr]

dt.meta.sel <- tr_meta[rw_yend >= crit_yr, c(unique(c("uid_radius", uid)), "longitude", "latitude",  "species"), with = FALSE]
dt.meta.sel[, lat:= round(latitude/resol_lat, 0)*resol_lat]
dt.meta.sel[,lon:=round(longitude/resol_lon, 0)*resol_lon]

uids.sll <- dt.meta.sel[, .N, by = c("lat", "lon", "species", uid)][,N:= NULL]
uids.spc <- uids.sll[,.N, by = .(species)]
setorder(uids.spc, -N)
uids.spc[, pct.species := N/sum(N)]
uids.spc[, pct.species := round(pct.species * 100,0)]
uids.spc[, ord:= .I]
if (!is.null(spc.lst)) {
  N.species <- NULL
  uids.sel <- uids.sll[uids.spc[species %in% spc.lst][,c("species", "pct.species", "N", "ord")], on = .(species)]
}
if (!is.null(N.species)) uids.sel <- uids.sll[uids.spc[ord <= N.species][,c("species", "pct.species", "N", "ord")], on = .(species)]
dist_uids <- uids.sel[, .(nuids = .N), by = .(ord, species, pct.species, N, lat, lon)]
if (make.plot){
  data.tmp <- dist_uids[ord <10]
  data.tmp[, spc.pct := paste0(species, " N", N, " ", pct.species, "%")]
  data.tmp$spc.pct <- factor(data.tmp$spc.pct, levels = unique(data.tmp$spc.pct[order(data.tmp$ord)]))

p1 <- ggplot(data.tmp, aes(x = lon, y = lat, size = nuids)) + facet_wrap(~spc.pct) +
  geom_point(alpha = 0.6, color = "darkblue") +
  scale_size_continuous(range = c(1, 10)) + # Adjust size range as needed
  scale_x_continuous(breaks = sort(unique(data.tmp$lon))) + # Set x-axis ticks to unique values of 'lat'

  theme_minimal() +
  theme(strip.text = element_text(size = 16),# Increase the size of facet labels
        panel.grid.minor = element_blank() , # Remove minor grid lines
        plot.title = element_text(size = 20), # Set title size
        plot.margin = margin(t = 10, r = 10, b = 30, l = 30, unit = "pt"), # Adjust plot margins
        plot.caption = element_text(hjust = 0, face = "italic")) + # Customize caption appearance

    labs(title = paste0( str_sub(uid, 5), " distribution yr",crit_yr),
      x = "Longitude",
       y = "Latitude",
       caption = "Data source: CFS-TRenD mar-2024" , # Add the data source caption
       size = paste0("n.", str_sub(uid, 5), "s"))

print(p1)
}

dist_uids <- dcast(dist_uids, species + pct.species + lat ~ lon, value.var = "nuids")
setorder(dist_uids, -pct.species, species, -lat)





if (out.csv == TRUE){
  if (is.null(out.dir)) message("please provide the directiory to output csv files")  else{
    if (!(dir.exists(out.dir))) dir.create(out.dir, recursive = TRUE)
    write.csv(dist_uids, file =  file.path(out.dir, paste0("count_", str_sub(uid, 5) ,".csv")), row.names = FALSE, na = "")
  }

}
return(dist_uids)
}
#
# dist_samples <- CFS_freq(list(dt.rwl, dt.meta), 1990,5, 10,  "uid_sample", 3)
# dist_sites <- CFS_freq(list(dt.rwl.all, dt.meta.all), 2000,5, 5,  "uid_site", 5)

# data <- list(dt.rwl.all, dt.meta.all)

# yend <- 1990; resol_lat <- 5; resol_lon <- 10; uid<-"uid_sample"; N.species <- 3
