library(data.table)

library(CFSTRenD)
species_nficode <- fread("P:/Jing/2010-08/Martin/Treering_bank/Git/tr/references/species_treesource.csv")

# province
loc.shp <- "P:/Jing/2010-08/Pierre/kNN 30m/Data treatment/J data"

require(rgdal)
require(ggplot2)
require(sf)
require(tmap)
shp <- read_sf(dsn = file.path(loc.shp, "Province_reprojected.shp"), stringsAsFactors = F)


usethis::use_data(species_nficode, internal = TRUE, overwrite = FALSE)
usethis::use_data(species_nficode, shp, internal = TRUE, overwrite = TRUE)

library(CFSTRenD)
library(data.table)
library(CFSTRenD)
# read ring measurement
samples69.o <- fread(system.file("extdata", "samples69.csv", package = "CFSTRenD"))
# read climate
clim69 <- fread(system.file("extdata", "clim69.csv", package = "CFSTRenD"))

load("R/sysdata.rda")
ls()  # Check if species_nficode is present

library(stringr)
samples69.trt <- CFS_format(data = list(samples69.o, 39:140), out.csv = NULL)
devtools::load_all()
rm(species_nficode)
CFSTRenD:::species_nficode


devtools::document()  # Update documentation.
devtools::build()     # Build the package tarball.
devtools::install()   # Install the package.
