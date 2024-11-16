## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(data.table)
library(CFSTRenD)
# read ring measurement
 samples69.o <- fread(system.file("extdata", "samples69.csv", package = "CFSTRenD"))
 # read climate
clim69 <- fread(system.file("extdata", "clim69.csv", package = "CFSTRenD"))


## -----------------------------------------------------------------------------

samples69.trt <- CFS_format(data = list(samples69.o, 39:140), out.csv = NULL)
samples69 <- samples69.trt$tr_all_wide
# print(samples69[, .N, by = .(latitude, longitude, uid_site)])


## ----fig.width=7, fig.height=5, fig.align='center'----------------------------
dt.freq <- CFS_freq(samples69, uid_level = "uid_tree",  geo_resolution = c(1,1))
plots_freq(dt.freq, out.species = "all", caption = "CFS-TRenD1.2 samples69", out.pdf= NULL )

## ----fig.width=7, fig.height=5, fig.align='center'----------------------------
# CFS_scale()
# specify the site and species to check
chk.site <- samples69[uid_site == 1 & species == "PSEUMEN",c("species", "uid_site", "site_id", "latitude", "longitude")][, .SD[1]]

dt.scale <- CFS_scale(site2chk = chk.site, ref_sites = samples69, N.nbs = 10, make.plot = list(TRUE, "CFS-TRenD1.2 samples69"))


## ----fig.width=7, fig.height=5, fig.align='center'----------------------------


 samples69_long <- merge(samples69.trt$tr_all_wide[, c("uid_site", "site_id", "species", "uid_tree", "uid_sample", "uid_radius")], 
   samples69.trt$tr_all_long$tr_7_ring_widths, by = "uid_radius")
# Calculate BAI
# rename to the reserved column name
setnames(samples69_long, c("uid_sample", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))
# check duplication
samples69_long[, .N, by = .(SampleID, Year)][N>1]
# run quality assessment procedure
qa.trt.ccf <-CFS_qa(dt.input = samples69_long)

# SampleID list
# unique(qa.trt.ccf$dt.ccf$SampleID)

# plot specific samples, figure output to screen
plots_qa(qa.trt.ccf$plot.lst, out.series = c("5"), out.pdf = NULL )

# plot all series, figures output to a pdf file
# plots_qa(qa.trt.ccf$plot.lst, out.series = "all", out.pdf = "c:/tmp/qa.pdf" )


## -----------------------------------------------------------------------------

samples69_clim <- merge(samples69.trt$tr_all_wide[, c("uid_site", "site_id", "species", "uid_tree", "uid_radius")], samples69.trt$tr_all_long$tr_7_ring_widths, by = "uid_radius")
# Calculate BAI
setorder(samples69_clim, uid_radius, year)
samples69_clim[, .N, by = .(uid_radius, year)][N>1]
samples69_clim[, `:=`(ageC = seq_len(.N),
  radius = cumsum(rw_mm),  # Cumulative radius (assumes RW is added each year)
          radius_prev = shift(cumsum(rw_mm), fill = 0)), by = uid_radius]  # Previous radius (shifted)

# Compute BA in cm2
samples69_clim[, ba_cm2_t_1 := pi * (radius_prev^2)/100]

samples69_clim[, bai_cm2 := pi * (radius^2 - radius_prev^2)/100]


# Drop the radius columns if not needed
samples69_clim[, c("radius", "radius_prev") := NULL]

samples69_clim <- merge(samples69_clim, clim69, by = c("site_id", "year"))

samples69_clim[, uid_tree.fac:= as.factor(as.character(uid_tree))]
setorder(samples69_clim, uid_tree, year)
# choose 1 site, ignore the first 10 years as Girardin 2022
site.i <- samples69_clim[uid_site == 10][ageC > 10]
det.i <- detrend_site(data = site.i, 
             m.candidates = c(
  "bai_cm2 ~ log(ba_cm2_t_1) + s(ageC)",
  "bai_cm2 ~ log(ba_cm2_t_1) + s(ageC) + FFD",
  "bai_cm2 ~ log(ba_cm2_t_1) + s(ageC) + s(FFD)"
  
)
)
print(summary(det.i$model$gam))

