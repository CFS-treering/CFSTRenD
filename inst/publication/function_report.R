
tr_6ispc <- tr_6i; tr_7ispc <- tr_7i

report_spc <- function(tr_6ispc, tr_7ispc, min.nseries){
  if (nrow(tr_6ispc) > 1) {
    tr_7i<- copy(tr_7ispc)
    setnames(tr_7i, c("uid_radius", "year", "rw_mm"), c("SampleID", "Year" ,"RawRing"))
    # check duplication
    tr_7i[, .N, by = .(SampleID, Year)][N>1]
    # include RW_trt, this is the series that qa process works on,
    tr_7i[, RW_trt:= RawRing - shift(RawRing), by = SampleID]; label_trt <- "differentiated"
    acf.trt <- tr_7i[!is.na(RW_trt), .( ar1_trt = round(acf(RW_trt, plot = FALSE)$acf[2],2) ), by = .(SampleID)]
    # run quality assessment procedure
    qa.trt.ccf <-CFS_qa(dt.input = tr_7i)
    dt.radii.spc <- copy(qa.trt.ccf$dt.stats)
    dt.radii.spc[, uid_radius:= as.integer(SampleID)]
    dt.radii.spc <- merge(tr_6ispc, dt.radii.spc, by = "uid_radius")

  }
# dt.radii.spc <- dt.radii[species == spc.lst[i.spc]]
dt.summary.0 <- tr_6ispc[, .(lat = paste0(round(min(latitude),1), " ", round(max(latitude),1)),
                          lon = paste0(round(min(longitude),1), " ", round(max(longitude),1))), by = .(species)]

setnames(dt.summary.0, c("lat", "lon"), c("Latitude range", "Longitude range"))
dt.summary.0a <- melt(dt.summary.0, measure.vars = names(dt.summary.0), variable.name = "Category", value.name = "Value")





dt.summary.1 <- dt.radii.spc[, .(Yspan = paste0(min(ymin), " ", max(ymax)), Nsites = length(unique(uid_site)), Ntrees = length(unique(uid_tree)),
                                 Nsamples = length(unique(uid_sample)), N = .N,
                                 summ_rw = paste0(round(mean(rw.mean),2), " ± ", round(sd(rw.mean),2), " ( ", round(min(rw.min),2), ", ", round(max(rw.max),2), " )"),
                                 summ_len = paste0(round(mean(N),2), " ± ", round(sd(N),2), " ( ", round(min(N),2), ", ", round(max(N),2), " )"),
                                 summ_ar1 = paste0(round(mean(ar1_rw),2), " ± ", round(sd(ar1_rw),2), " ( ", round(min(ar1_rw),2), ", ", round(max(ar1_rw),2), " )"))]

dt.summary.1[, (names(dt.summary.1)) := lapply(.SD, as.character)]
setnames(dt.summary.1, c("Yspan","Nsites", "Ntrees", "Nsamples", "N", "summ_rw", "summ_len", "summ_ar1"), c("Year span", "Number of sites", "Number of trees", "Number of samples", "Number of series", "summary rw(mm)**", "summary len.series**", "summary AR1**"  ))

dt.summary.1a <- melt(dt.summary.1, measure.vars = names(dt.summary.1), variable.name = "Category", value.name = "Value")

dt.summary.2a <- rbind(dt.summary.0a, dt.summary.1a)

if (nrow(tr_6ispc) > min.nseries) {
dt.Npass<- data.table(Category = "Number of series (pass)*", Value = nrow(dt.radii.spc[qa_code == "pass"]), ord = 8.5 )} else{
  dt.Npass<- data.table(Category = "Number of series (pass)*", Value = "N.A.", ord = 8.5 )
  }
dt.summary.2a[, ord:=.I]
dt.summary.3a <- rbind(dt.summary.2a, dt.Npass)
setorder(dt.summary.3a, ord)
dt.summary.3a[,ord:=NULL]
# dt.summary.3a$species <- unique(tr_6ispc$species)
return(list(summary_spc = dt.summary.3a, qa.trt.ccf = qa.trt.ccf ))
}


tr_w6 <- rwl.all.trt$tr_all_wide[, c("uid_site", "site_id", "uid_tree", "uid_sample",  "uid_radius", "radius_id", "rw_yend", "rw_ystart", "latitude", "longitude", "species", "dbh_cm", "ht_tot_m")]
tr_7 <- merge(tr_w6,rwl.all.trt$tr_all_long$tr_7_ring_widths, by ="uid_radius")

qa<- vector("list", length(spc.lst))
for (i.spc in seq_along(spc.lst)){

  tr_6i <- tr_w6[species == spc.lst[i.spc]]
  tr_7i <- tr_7[species == spc.lst[i.spc]]

  rpt.spc <- report_spc(tr_6i, tr_7i)
  dt.summary.3a <- rpt.spc[[1]]
  # qa in a list
  rpt.spc[[2]][[1]]$species <- spc.lst[i.spc]
  rpt.spc[[2]][[2]]$species <- spc.lst[i.spc]
  qa[[i.spc]] <- rpt.spc[[2]]
}
library(data.table)

# Assume list_of_lists contains 10 lists, each with 4 data.tables

dt.stats <- rbindlist(lapply(qa, `[[`, 4), use.names = TRUE, fill = TRUE)

dt.stats[, uid_radius:= as.integer(SampleID)]
# dt.radii <- merge(robj$tr_all_wide[,c("species", "uid_radius", "uid_site", "site_id", "uid_tree", "uid_sample")], dt.radii, by = "uid_radius")

i.spc <- 1
rw_ref <- rwl.all.trt$tr_all_wide ; N.nbs <- 5

report_spc_site <- function(tr_6i, rw_ref, N.nbs){

dt.site <-tr_6i[,c("uid_site", "site_id",  "species", "latitude", "longitude", "uid_tree", "uid_sample", "uid_radius")]
dt.site.stats <- dt.site[, .(Ntrees = length(unique(uid_tree)), Nsamples = length(unique(uid_sample)), Nradius = .N), by = .(uid_site, site_id, longitude, latitude, species)]

dt.radii.spc <- dt.radii[species == spc.lst[i.spc]]

dt.site.summ <- dt.radii.spc[ ,.(
  summ_len = paste0(round(mean(N),2), " ± ", round(sd(N),2), " ( ", round(min(N),2), ", ", round(max(N),2), " )"),
  summ_rw = paste0(round(mean(rw.mean),2), " ± ", round(sd(rw.mean),2), " ( ", round(min(rw.min),2), ", ", round(max(rw.max),2), " )")), by = .(uid_site)]
# setnames(dt.site.stats, c("summ_rw", "summ_len"), c("summary rw(mm)**", "summary len.series**"))
dt.site.stats <- merge(dt.site.stats, dt.site.summ, by = "uid_site")
head(dt.site.stats)
if (length(unique(dt.site.stats$uid_site)) < N.nbs + 2) {
  message("not sufficiant reference sites")
  dt.site.out <- dt.site.stats
  } else{
dt.scale.all <- data.table()
for (uid_site.i in dt.site.stats$uid_site){

  site.chk <- dt.site.stats[uid_site == uid_site.i ,c("species", "uid_site", "site_id", "latitude", "longitude")][, .SD[1]]

  dt.scale.i <- CFS_scale(site2chk = site.chk, ref_sites = rw_ref, N.nbs = N.nbs)$ratio.median
  dt.scale.all <- rbind(dt.scale.all, dt.scale.i)
  }

dt.site.out <- merge(dt.site.stats, dt.scale.all[, c("uid_site", "rw.median", "rw.median.nbs", "ratio_median")], by = "uid_site")
}
if (nrow(dt.site.out[!(Ntrees == Nradius)]) == 0) dt.site.out[, c("Nsamples", "Nradius") := NULL] else {
  if (nrow(dt.site.out[!(Nsamples == Nradius)]) == 0) dt.site.out[, c("Nradius") := NULL] else {
    if (nrow(dt.site.out[!(Nsamples == Ntrees)]) == 0) dt.site.out[, c("Nsamples") := NULL]
  }
}
return(dt.site.out)
}

tr_6i <- tr_w6[species == spc.lst[i.spc]]
dt.site.out <- report_spc_site(tr_6i, rw_ref, N.nbs = 5)
