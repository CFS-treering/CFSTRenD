## ----setup_general, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
   echo = TRUE,
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup_lib, include=FALSE-------------------------------------------------
if (!requireNamespace("htmltools", quietly = TRUE)) {
  stop("Package 'htmltools' is required to run this function. Please install it.")
}
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required to run this function. Please install it.")
}
library(CFSTRenD)
library(data.table)
library(ggplot2)
library(sf)



## ----data_qa, echo= FALSE, eval = FALSE---------------------------------------
# # Example: load a dataset included in the package
# 
# # ring measurement
# samples69.o <- fread(system.file("extdata", "samples69.csv", package = "CFSTRenD"))
# 
# # formatting the users' data conformed to CFSTRenD
# samples69.trt <- CFS_format(data = list(samples69.o, 39:140), usage = 1, out.csv = NULL)
# class(samples69.trt)

## ----scale--------------------------------------------------------------------
# loading processed data
samples69.trt <- readRDS(system.file("extdata", "samples69.trt.rds", package = "CFSTRenD"))
# , message = FALSE, warning = FALSE, results = 'hide'
all.sites <- samples69.trt$tr_all_wide[,.N, by = c("species", "uid_site", "site_id")][, N:=NULL]
if (nrow(all.sites[, .N, by = .(species, site_id)][N>1]) > 0) stop("species-site_id is not unique...")
# e.g. taking the target sites
target_site <- all.sites[c(1,2), -"uid_site"]

ref.sites <- merge(samples69.trt$tr_all_wide[,c("species", "uid_site", "site_id", "latitude","longitude", "uid_radius")], samples69.trt$tr_all_long$tr_7_ring_widths, by = c("uid_radius"))

dt.scale <- CFS_scale( target_site = target_site, ref_sites = ref.sites, scale.label_data_ref = "CFSTRenD V1.2-proj69", scale.max_dist_km = 20, scale.N_nbs = 5)

## ----scale_report_demo, eval= FALSE-------------------------------------------
# generate_report(robj = dt.scale, output_file = outfile_scale)

## ----scale_report, include = FALSE, message = FALSE, warning = FALSE, results = 'hide'----
outfile_scale <- tempfile(fileext = ".html")
generate_report(robj = dt.scale, output_file = outfile_scale)

## ----scale_inclu, echo = FALSE, results = 'asis'------------------------------
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9;">')
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9; 
#           max-height: 400px; overflow-y: scroll;">')


cat('
<div style="margin: 20px 0; padding: 10px; 
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            border-radius: 8px;">
')

htmltools::includeHTML(outfile_scale)  # embed in vignette
cat('</div>')

## ----qaa, message = FALSE, warning = FALSE, results = 'hide'------------------

# data processing
samples69_long <- merge(samples69.trt$tr_all_wide[, c("uid_site", "site_id", "species", "uid_tree", "uid_sample", "sample_id", "radius_id", "uid_radius")],
                        samples69.trt$tr_all_long$tr_7_ring_widths, by = "uid_radius")

names(samples69_long)
# rename to the reserved column name
setnames(samples69_long, c("sample_id", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))

# assign treated series
samples69_long[, RW_trt:= RawRing - shift(RawRing), by = SampleID]

dt.qa <-CFS_qa(dt.input = samples69_long, qa.label_data = "CFSTRenD V1.2-proj69", qa.label_trt = "difference", qa.min_nseries = 50)

## ----qa_report_demo, eval=FALSE, message = FALSE, warning = FALSE, results = 'hide'----
# 
# # e.g. series to check
# # chk.lst <- dt.qa$dt.stats[1:2,]$SampleID
# 
# generate_report(robj = dt.qa,  qa.out_series = c("X003_101_001", "X00._101_002"), output_file = NULL)

## ----qa_report, include= FALSE, message = FALSE, warning = FALSE, results = 'hide'----
chk.lst <- dt.qa$dt.stats[1:2,]$SampleID 

outfile_qa <- tempfile(fileext = ".html")
generate_report(robj = dt.qa,  qa.out_series = chk.lst, output_file = outfile_qa)

## ----qa_inclu, echo = FALSE, results = 'asis'---------------------------------
# cat('<div style="border: 5px solid #4682B4; border-radius: 8px; padding: 10px; 
#           margin: 15px 0; background-color: #f9f9f9; 
#           max-height: 400px; overflow-y: scroll;">')

cat('
<div style="margin: 20px 0; padding: 10px; 
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            border-radius: 8px;">
')

htmltools::includeHTML(outfile_qa)  # embed in vignette
cat('</div>')

