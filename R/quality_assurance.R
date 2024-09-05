# step 1: pair-wise ccf on all the samples. find all pairs that satisfies the condition that max_ccf @lag0
#   the results of this step serves as the initial sample list of master chronology for step 2
# step 2: run ccf between each sample and the mean of initial master chronology
#   this step is to check ccf between each sample and the the mean master chronology,
#   the results is the final master chronology satisfying the condition that max_ccf @lag0 with each sample that form the master
#   the out table contains the correlation with master at lag0 and its rank(from 1 to 21, 1 represents the maximum one, 21 represents minimum),
#   qa_code(Pass, borderline, highpeak, pm1) for each sample


# functions

#' tree ring data quality classification
#' @description
#' classify tree ring data in 4 classes ("pass","borderline", "highpeak", "pm1") using the detrended series based on 2 criteria: 1: treated
#'
#' @param dt.trt_wide : treated series in wide format for pair-wise ccf
#' @param dt.rw_long : ring width series in long format for calculating mean of master chronology
#' @param mtd.paired_ccf : processing approach for paired ccf c("S", "LB", "LP")
#' @param batch_size :  batch size for processing large dataset
#' @param max.lag : maximum lag for ccf
#' @param max.trial : maximum loop for step 2 to form the master dendrochronology

#'
#'
#'
#' @rawNamespace import(purrr, except = c(transpose))
#' @import data.table
#' @import parallel
#'
#' @return list including model, fitting statistics, ptable, stable and prediction table
#' @export tr.fullmaster
#'




tr.fullmaster <- function(dt.trt_wide , dt.rw_long , mtd.paired_ccf = "S", batch_size = 5000, max.lag = 10, max.trial = 100){
  if (!(mtd.paired_ccf %in% c("S", "LB", "LP"))) stop("please specify mtd.paired_ccf as one of c('S', 'LB', 'LP')")
  # dt.trt_wide: treated series in wide format for pair-wise ccf
  # dt.rw_long: ring width series in long format for calculating mean of master chronology

  # step 1: pair-wise ccf to find all the samples which can find at least 1 sample to reach max_ccf @ lag0

  # Generate all pairs of columns
  col_pairs <- combn(names(dt.trt_wide), 2, simplify = FALSE)

  # Use map to calculate pairwise cross-correlation for each pair
  # this is the most efficient way for this calculation so far i found, 20-06-24, 2 mins for 490 series
  # system.time(ccf.pairs <- map(col_pairs, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]])))


  # ccf.pairs <- map(col_pairs, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]]))
  # # Convert the list of ccf.pairs to a data.table using col_pairs
  # dt.ccf.pairs <- rbindlist(lapply(seq_along(col_pairs), function(i) {
  #   pair <- col_pairs[[i]]
  #   res <- ccf.pairs[[i]]
  #   data.table( ts1 = pair[1], ts2 = pair[2], max_lag = res$max_lag, max_ccf = res$max_ccf)
  # }))


  # Generate the data.table with the ccf.pairs and column pairs in one statement
  if (mtd.paired_ccf == "S"){
   dt.ccf.pairs <- rbindlist(map2(col_pairs,
                                 map(col_pairs, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]], max.lag = max.lag)),
                                 ~ {
                                   data.table(ts1 = .x[1], ts2 = .x[2], max_lag = .y$max_lag, max_ccf = .y$max_ccf)
                                 }))
  }else{

  # Number of columns
   n_cols <- length(names(dt.trt_wide))

   # Generate all pairs of columns
   col_pairs <- combn(names(dt.trt_wide), 2, simplify = FALSE)

   # Split pairs into batches
   pair_batches <- split(col_pairs, ceiling(seq_along(col_pairs) / batch_size))
  # batch process
   if (mtd.paired_ccf == "LB"){

     # Process each batch separately and combine results
     dt.ccf.pairs <- rbindlist(lapply(pair_batches, function(batch) {
       rbindlist(map2(batch,
                      map(batch, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]], max.lag = max.lag)),
                      ~ {
                        data.table(ts1 = .x[1], ts2 = .x[2], max_lag = .y$max_lag, max_ccf = .y$max_ccf)
                      }))
     }))


   }

   # parallel process for large dataset
   if (mtd.paired_ccf == "LP"){

    # Number of cores (use one less than the total number of cores)
    num_cores <- detectCores() - 1

    # Parallel processing of each batch and combining results
    dt.ccf.pairs <- rbindlist(mclapply(pair_batches, function(batch) {
      rbindlist(map2(batch,
                     map(batch, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]], max.lag = max.lag)),
                     ~ {
                       data.table(ts1 = .x[1], ts2 = .x[2], max_lag = .y$max_lag, max_ccf = .y$max_ccf)
                     }))
    }, mc.cores = num_cores))

  }

 }

  result_dt.sel <- dt.ccf.pairs[max_lag == 0 & !is.na(max_ccf)]

  ts.sel <- data.table(ts.sel = unlist(result_dt.sel$ts1, result_dt.sel$ts2))
  ts.sel <- ts.sel[, .N, by = .(ts.sel)]
  id.candi <- unique(ts.sel$ts.sel)

  # the result of step 1 is id.candi, it serves as the initial sample list of master chronology for step 2
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  # step 2: find the sample list of master chronology satisfying the condition that max_ccf @lag0 for each sample in this list with the master chronology

  # algorithm on pass

  s2.end <- FALSE;  i.iter <- 1;
  while(!s2.end & i.iter < max.trial){

    # mean of master chronology
    dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.candi][, .(.N, mean.rw = mean(RawRing)), by = .(Year)]
    setorder(dt.s2.avg, Year)
    dt.s2.avg [, mean.rw.dif:= mean.rw - shift(mean.rw)]

    dt.s2.wide <- merge(dt.s2.avg[, c("Year", "mean.rw.dif")], dt.trt.wide, by = "Year")
    # ccf of all samples with the master chronology
    dt.s2.ccf <-rbindlist(lapply(3:ncol(dt.s2.wide), ccf_avg, data = dt.s2.wide, max.lag = max.lag))
    # valid samples for master chronology
    id.pass <- unique(dt.s2.ccf[qa_code == "Pass"]$SampleID.chr)
    # id.dif <- setdiff(union(id.candi, id.pass) ,intersect (id.candi, id.pass))
    s2.end <- length(setdiff(union(id.candi, id.pass) ,intersect (id.candi, id.pass))) == 0
    print(paste0(i.iter, " N.pass: ", length(id.candi)))
    # selection list stops changing
    if ( s2.end != TRUE) {
      id.candi <- id.pass
      i.iter <- i.iter + 1

    }
  }
  setorder(dt.s2.ccf, SampleID.chr)
  dt.s2.avg[, c("s2.success", "iteration") := .(s2.end, i.iter)]
  return(list(dt.ccf = dt.s2.ccf, master = dt.s2.avg ))
  # the result of step 2 is dt.s2.ccf, samples with qa_code = "Pass" to form the master chronology

}



# Define a function to calculate max cross-correlation lag
ccf_pairs <- function(ts1, ts2, max.lag = 10) {
  ccf.chk <- ccf(ts1, ts2, lag.max = max.lag, na.action = na.pass,  plot = FALSE)
  max_ccf <- max(ccf.chk$acf)
  max_lag <- ccf.chk$lag[which.max(ccf.chk$acf)]
  return(list(max_lag = max_lag, max_ccf = max_ccf))
}

ccf_avg <- function(icol, data, blcrit = 0.1, max.lag = 10, qa_code="Fail"){

  dt.pairs <- ccf(data[,2],data[,icol, with = FALSE],lag.max=max.lag,plot=FALSE, na.action = na.pass)
  dt.ccf <- data.table(SampleID.chr = names(data)[icol], lag = as.vector(dt.pairs$lag), acf.trt = as.vector(dt.pairs$acf))
  setorder(dt.ccf, -acf.trt)
  dt.ccf[, ccf.ord:= 1:.N]
  setorder(dt.ccf, lag)
  # max acf at lag 0
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 1 ) qa_code="Pass"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 2]) == 1 & abs(dt.ccf[ccf.ord == 1]$acf.trt - dt.ccf[ccf.ord == 2]$acf.trt) <blcrit) qa_code="borderline"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 0 & dt.ccf[ccf.ord == 1]$acf.trt / dt.ccf[ccf.ord == 2]$acf.trt > 2) qa_code="highpeak"
  if (nrow(dt.ccf[lag == 1 & ccf.ord == 1]) + nrow(dt.ccf[lag == -1 & ccf.ord == 1]) == 1) qa_code="pm1"
  dt.ccf$qa_code <- qa_code
  return(dt.ccf)
}




