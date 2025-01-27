# to avoid warnings in checking for data.table variables
utils::globalVariables(c( ":=", ".", "!!"))
utils::globalVariables(c( "Format", "FormatV"))

utils::globalVariables(c(
  "N", "Required", "Variable", "latitude", "longitude",
  "meas_no","project_name", "radius_id", "sample_id", "site_id", "tr_1_projects", "tr_2_sites",
  "tr_3_trees", "tr_4_meas", "tr_5_samples", "tr_6_radiuses", "tree_id", "uid_meas",
  "uid_project", "uid_radius", "uid_sample", "uid_site", "uid_tree",
  "coord", "rwinc", "site.tree", "species_nficode", "ydif",
  "aic",  "decade.yr",  "dif.decade",   "fit.s.ageC", "fn",   "i",
  "id.core",    "id2",  "lines","nc.yr",  "rest", "rw", "i.na",
  "rw.N", "rw.rest",    "rw.sc","rw.tmp",     "rw_int" ,    "rw_mm","se.fit.bai",
  "se.fit.s.ageC", "selected",   "site_ID1.01.0", "specialchar",   "v.idx",
  "lat", "lon", "meas_date", "ord", "species", "uid_radius.tmp", "value", "variable",
  "year_range", "ymax", "ymax.proj", "ymin", "ymin.proj", "yr.meas",
  "rw_yend", "rw_ystart",

  # CFS_scale
  "rw.median",  "yr.max", "yr.mn", "ratio_median", "rw.median.nbs", "size_class",
  "list.tbl", "pct.species", "spc.pct", "nuids", "bark_thickness_mm",
  "block_id","rw.N.sc",

  # CFS_qa
  "acf.trt", "ccf.ord", "max_lag", "max_ccf", "SampleID.chr", "RawRing", "RW_trt", "Year",
  "mean.rw.dif","mean.rw", "dt.trt.wide", "qa_code", "col.ord",
  "SampleID",  "colr", "id.label", "rw.treated",

  # gamm_main
  "ageC", "res.normalized", "res.normalized.LL",  "res.resp_normalized",
  "lat_use", "lon_use", "start.event",
  # cal.bai
  "ba_cm2_t_1", "bai_cm2", "radius", "radius_prev", "med.rw"
  )
)


