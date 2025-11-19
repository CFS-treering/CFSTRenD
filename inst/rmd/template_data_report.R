params <-
list(qa.min_nseries = 50, scale.N_nbs = 10)

## ----setup_data, include=FALSE------------------------------------------------
# This chunk sets global options to hide all R code from appearing in the output
if (!requireNamespace("htmltools", quietly = TRUE)) {
  # stop("Package 'htmltools' is required to run this function. Please install it.")
  install.packages("htmltools")
}
knitr::opts_chunk$set(echo = FALSE)
library(data.table)
library(growthTrendR)
library(gt)
# library(htmltools)
library(stringr)
# 
 robj <- params$robj
 reports.sel <- params$data_report.reports_sel
 label_data <- params$qa.label_data
 label_trt <- params$qa.label_trt
 qa.min_nseries <- params$qa.min_nseries

 # max.dist_km <- params$parms.qa["max.dist_km"]
 N.nbs.qa <- params$scale.N_nbs

 
 # for test
 
 
  if (!inherits(robj, "cfs_format")) stop("please check the input , make sure it's the result of CFS_format() function")

# library(ggplot2)
# library(sf)

# read ring measurement
 # samples69.o <- fread(system.file("extdata", "samples69.csv", package = "growthTrendR"))
 #  samples69.trt <- CFS_format(data = list(samples69.o, 39:140), usage = 1, out.csv = NULL)
 #
#  f <- system.file("extdata", "IDF.trt.rdata", package = "growthTrendR")
# load(f)
# robj <- IDF.trt; reports.sel <- 1:4; qa.min_nseries <- 20; rep_site <- 3; rep_species <- 2; N.nbs.qa <- 5; i.spc <- 1
# usage <- params$usage
  tr_long <- robj$tr_all_long
  
 tr_w6 <- robj$tr_all_wide[, c("uid_site", "site_id", "uid_tree", "uid_sample",  "uid_radius", "radius_id", "rw_yend", "rw_ystart", "latitude", "longitude", "species", "dbh_cm", "ht_tot_m")]
tr_7 <- merge(tr_w6, tr_long$tr_7_ring_widths, by ="uid_radius")

spc.lst <- unique(tr_long$tr_3_trees$species)

dt.freq <- CFS_freq(tr_w6) 
dist.uids <- melt(dt.freq$dist_uids, id.vars = c("uid_label", "species" , "ord", "N", "pct.species" , "lat"),
                  variable.name = "lon",
                  value.name = "nuids")[!is.na(nuids)]
dist.uids[, lon:= as.numeric(as.character(lon))]
setorder(dist.uids, ord, lon, lat)

# label_trt <- "differentiated"

# min.nseries.qa <- 20; N.nbs.qa <- 6;
# by species
# qa<- vector("list", length(spc.lst))
# for (i.spc in seq_along(spc.lst)){
# print(i.spc)

rep_proj <- 1; rep_species <- 2; rep_site <- 3; rep_radii <- 4

