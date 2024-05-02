#' prepare variable names and column order for each tables
#'
#'
#' @import data.table
#' @import stringr
#'
#' @return two lists: meta.all and ord.all
#' @export
#'
#' @examples
#' meta.all <- variables.cfstrend()
variables.cfstrend <- function(){
  # data structure version, it changes with the structure changes,
  Ver.stru <- "s1.2"
meta.Projects <-data.table(Variable = c( "uid_project", "submission_id", "project_name", "description",  "year_range", "reference", "open_data",  "contact1", "contact2"),
                         Format = c("integer", "integer", "character",  "character",  "character", "character",  "logical","character","character"),
                          Description = c("unique project identification used in treeSource","submission identification", "original project identification",  "Project description",
                                          "year range of the project", "References", "open to public","Name of primary person responsible for the data","Name of secondary person responsible for the data"),
                         Required = c(0,1,1,2,2,2,1,1,2)
)

meta.Sites <-data.table(Variable = c( "uid_site",  "site_id",    "latitude" ,      "longitude",  "datasource", "investigators"),
                           Format = c("integer", "character",     "numeric",         "numeric", "character", "character"),

                           Description = c("unique site identification used in treeSource",
                                           "original site identification",   "  decimal degrees (positive for N)",
                                  "decimal degrees (negative for W)", "inventory/target", "Investigators"),
                           Required = c(0,1,1,1,2,2)
)


meta.Trees <-data.table(Variable = c( "uid_tree", "uid_site", "tree_id",  "species", "uid_project"),
                       Format = c("integer","integer","character",  "character","integer"),
                       Description = c("unique tree identification used in treeSource", "unique site identification used in treeSource","original tree identification",
                                       "use 7 character species codes from NFI with no spaces (e.g.; PICEMAR; POPUTRE)", "unique project identification used in treeSource"),
                       Required = c(0,0,1,1,0)
                      )



meta.Meas <-data.table(Variable = c( "uid_meas", "uid_tree",  "meas_no",   "meas_date", "status",    "dbh_cm",    "ht_tot_m"),
                       Format = c("integer","integer", "integer", "character",  "character", "numeric","numeric"),
                       Description = c("unique measurement identification used in treeSource", "unique tree identification used in treeSource",
                                       "original measurement identification, 0 for missing info", "measurement date, last year in tree_RW in case missing",
                                       "health status: LIVE/DEAD","Stem diameter (cm) at breast height (ca. 1.3 m)", "Total tree height (m)"),
                       Required = c(0,0,1,1,2,2,2)
)


meta.Samples <-data.table(Variable = c( "uid_sample", "uid_meas", "sample_id", "sample_type",   "sample_ht_m",     "sample_diameter_cm"),
                         Format = c( "integer",  "integer", "character", "character",    "numeric",   "numeric"),
                         Description = c("unique sample identification used in treeSource", "unique measurement identification used in treeSource", "sample ID",
                                        "sample type: core/disk",  "height at point of sample (m)", "diameter at point of sample (cm)"),
                         Required = c(0,0,1,2,2,2)
)


meta.Radius <-data.table(Variable = c( "uid_radius", "uid_sample", "radius_id", "cofecha_id", "ring_meas_method",
                                      "crossdating_visual", "crossdating_validation", "age_corrected", "bark_thickness_mm", "radius_inside_cm",
                                      "dtc_measured_mm", "dtc_estimated_mm", "comments" ),
                           Format = c( "integer",  "integer", "character",  "character", "character",
                                       "logical", "character", "numeric",  "numeric", "numeric",
                                        "numeric" , "numeric","character"),
                           Description = c("unique radius/core identification used in treeSource", "unique sample identification used in treeSource",
                                          "original radius/core identification, \"O\" for missing", "identification used in cofecha",  "ring measurement method: windendro/velmex/coorecorder",
                                          "visaully cross-dated? yes/no", "validation cross-dated by: CDENDRO/COFECHA",
                                          "Corrected age (with rings added based on fresh DTC divided by fresh pith increment)",
                                          "thickness of the bark (mm)", "Core length corrected to Fresh tree dimensions outside bark (cm)",
                                          "Measured distance from start of earliest ring to centre (mm)",
                                          "Gap-filled distance from start of earliest ring to centre (mm)","user notes on cores"),
                           Required = c(0,0,1,2,2,2,2,2,2,2,2,2,2)
                           )


# v1.2 uid_ringwidth not necessary, using uid_radius + year as primary key
meta.Ringwidths <-data.table(Variable = c("uid_radius", "year",  "rw_mm"),
                         Format = c(  "integer",  "integer", "numeric"),
                         Description = c( "unique radius/core identification used in treeSource",
                                         "Year", "ring width (mm)"),
                         Required = c(0,1,1)
)

# V1.2 uids being removed
meta.uids_removed <-data.table(Variable = c("modification_id", "uid_removed", "uid_level",  "submission_id", "uid_project",  "ver_firstremoved",  "ver_lastexist", "reason", "investigator"),
                               Format = c(  "integer",  "integer",  "character", "integer",  "integer", "character","character","character","character"),
                               Description = c("modification id", "uid removed", "the level of uid below which all uids will be removed('uid_XXX')", "submission_id which uids were removed",
                                                "uid_project in which uids were removed","the first version in which uid_removed disappear from the dataset",
                                                "the last version in which uid_removed still exist in the dataset", "reason for removing the uid_removed", "investigator"),
                               Required = c(1,0,1,1,0,1,1,1,1)
)

# BEGIN;
# CREATE TABLE tr.tr_8_uids_removed (
#   modification_id integer PRIMARY KEY,
#   uid_removed integer,
#   uid_level TEXT,
#   submission_id integer,
#   uid_project integer,
#   ver_firstremoved TEXT,
#   ver_lastexist TEXT,
#   reason TEXT,
#   investigator TEXT
# );
#
# COMMENT ON COLUMN tr.tr_8_uids_removed.modification_id IS 'modification id.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.uid_removed IS 'uid removed.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.uid_level IS 'the level of uid below which all uids will be removed(uid_XXX).';
# COMMENT ON COLUMN tr.tr_8_uids_removed.submission_id IS 'submission_id which uids were removed.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.uid_project IS 'uid_project in which uids were removed.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.ver_firstremoved IS 'the first version in which uid_removed disappear from the dataset.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.ver_lastexist IS 'the last version in which uid_removed still exist in the dataset.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.reason IS 'reason for removing the uid_removed.';
# COMMENT ON COLUMN tr.tr_8_uids_removed.investigator IS 'investigator.';
# COMMIT;


# ts.lst <- c("tr_1_projects", "tr_2_sites", "tr_3_trees", "tr_4_meas", "tr_5_samples", "tr_6_radiuses", "tr_7_ring_widths")

meta.all <- rbind(meta.Projects[,table := "tr_1_projects"],
                  meta.Sites[,table := "tr_2_sites"],
                  meta.Trees[,table := "tr_3_trees"],
                  meta.Meas[,table := "tr_4_meas"],
                  meta.Samples[,table := "tr_5_samples"],
                  meta.Radius[,table := "tr_6_radiuses"],
                  meta.Ringwidths[,table := "tr_7_ring_widths"],
                  meta.uids_removed[,table := "tr_8_uids_removed"]
                  )
meta.all$ver.structure <- Ver.stru
# print(is.data.table(meta.all))
# Separate the elements with commas

meta.all[Format == "integer", FormatV:= "NA_integer_"]
meta.all[Format == "character", FormatV:= "NA_character_"]
meta.all[Format == "numeric", FormatV:= "NA_real_"]
meta.all[Format == "logical", FormatV:= "NA"]

# meta.all <- meta.all[!str_detect(Variable, "uid_")]

# ord
# {
#   ord.project <- as.character(meta.Projects$Variable)
#   ord.site <- as.character(meta.Sites$Variable)
#   ord.tree <- as.character(meta.Trees$Variable)
#   ord.meas <- as.character(meta.Meas$Variable)
#   ord.sample <- as.character(meta.Samples$Variable)
#   ord.radius <- as.character(meta.Radius$Variable)
#   ord.width <- as.character(meta.Ringwidths$Variable)
#   ord.all <- list(ord.project, ord.site, ord.tree, ord.meas, ord.sample, ord.radius, ord.width)
# }
return(meta.all)
}


# update.YN <- 0
#
# if (update.YN == 1){
#   str.require <- data.table(Value = c(0,1,2), comments = c("treeSource ID", "must-have variables", "good-to-have variables"))
# # rm(list=ls(pattern="^meta\\."))
# # update data structure
# list.meta <- c("meta.Projects", "meta.Sites", "meta.Trees", "meta.Meas", "meta.Samples", "meta.Radius", "meta.Ringwidths")
# require(XLConnect)
# wb <- XLConnect::loadWorkbook("E:/Git/tr/documents/data structure.xlsx", create = TRUE)
#
# #creating sheets within an Excel workbook
# for (i in 1:7){
#   XLConnect::createSheet(wb, name = list.meta[i])
#   clearSheet(wb, sheet = list.meta[i])
#   XLConnect::writeWorksheet(wb, eval(parse(text = list.meta[i])), sheet = list.meta[i], startRow = 1, startCol = 1)
#   XLConnect::setColumnWidth(wb, sheet = list.meta[i], column = c(1,2,3,4), width = c(6000, 3000, 18000, 2500))
#
# }
# }
# #writing into sheets within an Excel workbook :
# XLConnect::createSheet(wb, name = "info")
# XLConnect::writeWorksheet(wb, paste0("Version: ", Ver.df), sheet = "info",  header = FALSE,startRow = 2, startCol = 1)
# XLConnect::writeWorksheet(wb, "Required value", sheet = "info", header = FALSE, startRow = 5, startCol = 1)
# XLConnect::writeWorksheet(wb, str.require, sheet = "info", header = FALSE, startRow = 6, startCol = 1)
# setColumnWidth(wb, sheet = "info", column = c(1,2), width = c(3500, 8000))
# XLConnect::saveWorkbook(wb)
#
# }




