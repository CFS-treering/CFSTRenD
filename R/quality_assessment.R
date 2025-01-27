

# functions


# step 1: pair-wise ccf on all the samples. find all pairs that satisfies the condition that max_ccf @lag0
#   the results of this step serves as the initial sample list of master chronology for step 2
# step 2: run ccf between each sample and the mean of initial master chronology
#   this step is to check ccf between each sample and the the mean master chronology,
#   the results is the final master chronology satisfying the condition that max_ccf @lag0 with each sample that form the master
#   the out table contains the correlation with master at lag0 and its rank(from 1 to 21, 1 represents the maximum one, 21 represents minimum),
#   qa_code(Pass, borderline, highpeak, pm1) for each sample


#' tree-ring data measurement assessment
#' @description
#' Assess tree-ring measurement accuracy using a treated series based on the differences between two consecutive tree-ring measurements.
#'
#' @param dt.input tree ring data with at least 3 columns (SampleID, Year, RawRing)
#' @param batch_size  number of pairs to run in a batch, to avoid memory issues in processing large dataset
#' @param max.lag maximum lag up to which the correlation should be calculated in CCF
#' @param max.iter maximum number of iterations of step 2(see Details)

#'
#'
#'

#' @import data.table

#' @import future
#' @import furrr
#' @import parallel
#' @importFrom purrr map map2
#'
#' @return A list of 3 elements:
#' 1) dt.ccf:	A data table containing the CCF results for all samples, including the quality assessment code (qa_code).
#' 2) master:	The final master chronology, including both raw master chronology, the mean of ring measurement of the series with pass,
#' and the treated  master chronology, calculated as the difference of two consecutive raw master chronology.
#' 3) plot.lst:	plots of the raw and treated ring measurements with year, along with the CCF plots for each sample.


#'
#' @details Assess tree-ring measurement accuracy using a treated series based on the differences between two consecutive tree-ring measurements.
#' The algorithm consists of two main steps:
#'
#'  Step 1. Perform pairwise CCF on the treated series of all possible combinations of the samples.
#' raw ring master chronology is calculated as the average of the ring measurements of samples that achieve the maximum correlation at lag 0 with at least one other sample in the cross-correlation function (CCF).
#'  The difference series of the raw ring master chronology is used as the initial treated master chronology for next step.
#'
#'  Step 2. Perform CCF between the treated series of each sample and the initial treated master chronology.
#'  Samples that do not meet the criteria will be removed from the recalculation of the treated master chronology.
#'  This step is repeated until all remaining samples in the master chronology meet the criteria.
#'
#' This process results in five categories classification for all samples (pass, borderline, pm1, highpeak, fail)
#'
#'
#' To enhance efficiency and mitigate potential memory issues,
#' The function supports both parallel (multi-session) and sequential modes, and also offers a batch processing option for users.

#' @export CFS_qa
#'

CFS_qa <- function(dt.input , batch_size = 10000, max.lag = 10, max.iter = 100){

  if (length(setdiff(c("SampleID", "Year","RawRing", "RW_trt" ), names(dt.input))) > 0) stop("at least one of the mandatory columns (SampleID, Year, RawRing, RW_trt) doesn't exist, please check...")


  if (nrow(dt.input[, .N, by = .(SampleID, Year)][N > 1]) > 0) stop("SampleID-Year is not unique key, please check...")
  dt.rw_long <- dt.input[, c("SampleID", "Year","RawRing", "RW_trt")]
  setorder(dt.rw_long, SampleID,Year)
  # the series to be used for ccf should be in dt.input 2024-12-10
  # dt.rw_long[, RW_trt:= RawRing - shift(RawRing), by = SampleID]

  # SampleID.chr is the key sampleID, we use it in the functions to represent a sample.
  # starting with a character to ensure the validity as column name and value of a column

  dt.rw_long[, SampleID.chr:= paste0("d_", SampleID)]
  setorder(dt.rw_long, SampleID.chr)
  sample.lst <- sort(unique(dt.rw_long$SampleID.chr))

  dt.rw_wide <- dcast(dt.rw_long[!is.na(RW_trt)], Year ~ SampleID.chr, value.var = "RawRing")
  dt.trt_wide <- dcast(dt.rw_long[!is.na(RW_trt)], Year ~ SampleID.chr, value.var = "RW_trt")

  setcolorder(dt.trt_wide, c("Year", sample.lst))



  dt.trt_wide.o <- copy(dt.trt_wide); dt.trt_wide$Year <- NULL;
  # dt.trt_wide: treated series in wide format for pair-wise ccf
  # dt.rw_long: ring width series in long format for calculating mean of master chronology

  # step 1: pair-wise ccf to find all the samples which can find at least 1 sample to reach max_ccf @ lag0


  # Generate all pairs of columns
  col_pairs <- combn(names(dt.trt_wide), 2, simplify = FALSE)

  # Detect available cores for parallel processing
  available_cores <- detectCores(logical = FALSE) - 1  # Adjusted cores based on system

  # Decide if parallel processing is supported
  if (available_cores > 1) {
    plan(multisession, workers = available_cores)
  } else {
    plan(sequential)
  }

  # if batch size is not specified, run as 1 batch
  if (is.null(batch_size)) {
    pair_batches <- list(col_pairs)
  }else{
    pair_batches <- split(col_pairs, ceiling(seq_along(col_pairs) / batch_size))
  }

  # Run processing on batches with or without parallel
 cat("Progress pair-wise ccf...\n")
  dt.ccf.pairs <- future_map(pair_batches, process_batch, dt.wide = dt.trt_wide, max.lag = 10, .progress = TRUE) %>% rbindlist()


  cat("\n")

  result_dt.sel <- dt.ccf.pairs[max_lag == 0 & !is.na(max_ccf)]

  ts.sel <- data.table(ts.sel = unlist(result_dt.sel$ts1, result_dt.sel$ts2))
  ts.sel <- ts.sel[, .N, by = .(ts.sel)]
  id.candi <- unique(ts.sel$ts.sel)

  # the result of step 1 is id.candi, it serves as the initial sample list of master chronology for step 2
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  # step 2: find the sample list of master chronology satisfying the condition that max_ccf @lag0 for each sample in this list with the master chronology

  # algorithm on pass

  s2.end <- FALSE;  i.iter <- 1;
  while(!s2.end & i.iter <= max.iter){

    # mean of master chronology
    # dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.candi][, .(.N, mean.rw = mean(RawRing)), by = .(Year)]
    # setorder(dt.s2.avg, Year)
    # dt.s2.avg [, mean.rw.dif:= mean.rw - shift(mean.rw)]

    # mean of treated chronology as master chronology 2024-12-10
    dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.candi][, .(.N, mean.rw = mean(RawRing), mean.rw_trt = mean(RW_trt, na.rm = TRUE)), by = .(Year)]

    dt.s2.wide <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")
    # ccf of all samples with the master chronology
    dt.s2.ccf <-rbindlist(lapply(3:ncol(dt.s2.wide), ccf_avg, data = dt.s2.wide, max.lag = max.lag, qa_code = "Fail"))
    # valid samples for master chronology
    id.pass <- unique(dt.s2.ccf[qa_code == "pass"]$SampleID.chr)
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
  dt.s2.avg[, c("success", "iteration") := .(s2.end, i.iter)]

  # for generating plots
  # pre series data for both rw and treated in wide format
  dt.raw.series <- merge(dt.s2.avg[, c("Year", "mean.rw")], dt.rw_wide, by = "Year")
  dt.trt.series <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")

  setcolorder(dt.raw.series, c("Year", "mean.rw", sample.lst))
  setcolorder(dt.trt.series, c("Year", "mean.rw_trt", sample.lst))


  # pre data for bar plotting on ccf with master for raw and treated series

  dt.trt.ccf <- copy(dt.s2.ccf)
  dt.ccf.idlabel <- dt.trt.ccf[ccf.ord==1, c("SampleID.chr", "qa_code", "lag")]
  dt.ccf.idlabel[, id.label:= paste0(str_sub(SampleID.chr, 3, -1),"$", qa_code,"$", lag)]
  dt.trt.ccf<- merge(dt.trt.ccf, dt.ccf.idlabel[, c("SampleID.chr", "id.label")],by = "SampleID.chr")
  dt.trt.ccf <- dcast(dt.trt.ccf, lag ~ id.label, value.var = "acf.trt")
  names(dt.trt.ccf)
  idlabel.lst <- sort(unique(dt.ccf.idlabel$id.label))
  # test if in the same order as others
  idlabel.lst2 <- str_split_fixed(idlabel.lst, "\\$",3)[,1]
  if (!all.equal(idlabel.lst2,str_split_fixed(sample.lst, "\\_",2)[,2]  )) print("check the order of id.label in dt.ccf.idlabel")
  setcolorder(dt.trt.ccf, c("lag", idlabel.lst))

  # input data structure for ccf_avg Year, mean.value, sampleIDs...

  dt.raw.ccf <- rbindlist(lapply(3:ncol(dt.raw.series), ccf_avg, data = dt.raw.series))
  dt.raw.ccf <- dcast(dt.raw.ccf, lag ~ SampleID.chr, value.var = "acf.trt")
  names(dt.raw.ccf)
  setcolorder(dt.raw.ccf, c("lag", sample.lst))



  # ccf of all samples with the master chronology
  plot.trt.series <- lapply(3:ncol(dt.trt.series), create_plot.series, data = dt.trt.series)
  plot.raw.series <- lapply(3:ncol(dt.raw.series), create_plot.series, data = dt.raw.series)
  names(plot.trt.series) <- str_split_fixed(colnames(dt.trt.series)[3:ncol(dt.trt.series)], "\\_",2)[,2]
  names(plot.raw.series) <- str_split_fixed(colnames(dt.raw.series)[3:ncol(dt.raw.series)], "\\_",2)[,2]
  # print(plot.trt.series[[2]])
  # print(plot.raw.series[[2]])

names(plot.trt.series)

  plot.trt.ccf <- lapply(2:ncol(dt.trt.ccf), create_barplot, data = dt.trt.ccf)
  plot.raw.ccf <- lapply(2:ncol(dt.raw.ccf), create_barplot, data = dt.raw.ccf)

  names(plot.trt.ccf) <- str_split_fixed(colnames(dt.trt.ccf)[2:ncol(dt.trt.ccf)], "\\_",2)[,2]
  names(plot.raw.ccf) <- str_split_fixed(colnames(dt.raw.ccf)[2:ncol(dt.raw.ccf)], "\\_",2)[,2]

  # print(plot.trt.series[[2]])
  # print(plot.raw.ccf[[2]])
  # print(plot.trt.ccf[[2]])
  # plot.lst <- list(plot.raw.series, plot.trt.series, plot.raw.ccf, plot.trt.ccf )
  plot.lst <- list(plot.raw.series = plot.raw.series, plot.trt.series = plot.trt.series, plot.raw.ccf = plot.raw.ccf, plot.trt.ccf = plot.trt.ccf)
  # return(list(plot.raw.series = plot.raw.series, plot.trt.series = plot.trt.series, plot.raw.ccf = plot.raw.ccf, plot.trt.ccf = plot.trt.ccf))


  # for statistics per radius
  dt.s2.ccf[, SampleID := str_split_fixed(SampleID.chr, "\\_",2)[,2] ]


  # reports on radii
  dt.radii <- dt.rw_long[, .(N = .N, rw.mean = mean(RawRing), rw.sd = sd(RawRing), rw.min = min(RawRing), rw.max = max(RawRing), ymin = min(Year), ymax = max(Year), ar1_rw = acf(RawRing, plot = FALSE)$acf[2] ), by = .(SampleID.chr)]

  acf.trt <- dt.rw_long[!is.na(RW_trt), .( ar1_trt = round(acf(RW_trt, plot = FALSE)$acf[2],2) ), by = .(SampleID.chr)]

  # correlations
  dt_wide.rw <- merge(dt.s2.avg[, c("Year", "mean.rw")], dt.rw_wide, by = "Year")
  dt_wide.trt <- merge(dt.s2.avg[, c("Year", "mean.rw_trt")], dt.trt_wide.o, by = "Year")

  dt.cor <- merge(cor_avg(dt_wide.rw), cor_avg(dt_wide.trt), by = "SampleID.chr")

  stats_radii <- merge(dt.radii, acf.trt, by = "SampleID.chr")


  stats_radii <- merge(stats_radii, dt.cor, by = "SampleID.chr")


  stats_radii <- merge(dt.s2.ccf[lag == 0, c("SampleID", "SampleID.chr", "qa_code")], stats_radii, by = "SampleID.chr")


  dt.s2.ccf[, SampleID.chr := NULL]
  setcolorder(dt.s2.ccf, "SampleID")
  stats_radii[, SampleID.chr := NULL]

  # Reset to sequential
  plan(sequential)

  return(list(dt.ccf = dt.s2.ccf, master = dt.s2.avg, plot.lst = plot.lst, dt.stats = stats_radii))
  # the result of step 2 is dt.s2.ccf, samples with qa_code = "Pass" to form the master chronology

}



# Define a function to calculate max cross-correlation lag

#' pair-wise ccf
#' @description
#' pair-wise ccf
#' @param ts1 first series
#' @param ts2 second series
#' @param max.lag maximum lag for ccf


#' @import data.table

#'

#' @export ccf_pairs
ccf_pairs <- function(ts1, ts2, max.lag = 10) {
  ccf.chk <- ccf(ts1, ts2, lag.max = max.lag, na.action = na.pass,  plot = FALSE)
  max_ccf <- max(ccf.chk$acf)
  max_lag <- ccf.chk$lag[which.max(ccf.chk$acf)]
  return(list(max_lag = max_lag, max_ccf = max_ccf))
}


#' ccf with master/mean chronology
#' @description
#' ccf with master/mean chronology
#' @param icol column index
#' @param data data in wide format, first column is "Year", second column is master/mean chronology
#' @param blcrit criteria for borderline
#' @param max.lag maximum lag for ccf
#' @param qa_code quality classification code (NULL, no code)

#' @import data.table

#'

#' @export ccf_avg
ccf_avg <- function(icol, data, blcrit = 0.1, max.lag = 10, qa_code="fail"){

  dt.pairs <- ccf(data[,2],data[,icol, with = FALSE],lag.max=max.lag,plot=FALSE, na.action = na.pass)
  dt.ccf <- data.table(SampleID.chr = names(data)[icol], lag = as.vector(dt.pairs$lag), acf.trt = as.vector(dt.pairs$acf))
  setorder(dt.ccf, -acf.trt)
  dt.ccf[, ccf.ord:= 1:.N]
  setorder(dt.ccf, lag)

  if (!is.null(qa_code)){
    qa_code <- "Fail"
  # max acf at lag 0
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 1 ) qa_code="pass"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 2]) == 1 & abs(dt.ccf[ccf.ord == 1]$acf.trt - dt.ccf[ccf.ord == 2]$acf.trt) <blcrit) qa_code="borderline"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 0 & dt.ccf[ccf.ord == 1]$acf.trt / dt.ccf[ccf.ord == 2]$acf.trt > 2) qa_code="highpeak"
  if (nrow(dt.ccf[lag == 1 & ccf.ord == 1]) + nrow(dt.ccf[lag == -1 & ccf.ord == 1]) == 1) qa_code="pm1"
  dt.ccf$qa_code <- qa_code
  }
  return(dt.ccf)
}


# Define the main function to process each pair

# Use map to calculate pairwise cross-correlation for each pair
# this is the most efficient way for this calculation so far i found, 20-06-24, 2 mins for 490 series
# system.time(ccf.pairs <- map(col_pairs, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]])))

#' @keywords internal
process_batch <- function(batch,dt.wide, max.lag = 10) {
  rbindlist(
    map2(batch,
         map(batch, ~ ccf_pairs(dt.wide[[.x[1]]], dt.wide[[.x[2]]], max.lag = max.lag)),
         ~ {
           data.table(ts1 = .x[1], ts2 = .x[2], max_lag = .y$max_lag, max_ccf = .y$max_ccf)
         }
    )
  )
}

#' @keywords internal
create_plot.series <- function(icol, data) {
  if (str_detect(deparse(substitute(data)), "trt")) {

    label <- "treated "}else{

      label <- "raw"}
  sampleID<- names(data)[icol]
  p <- ggplot(data, aes(x = Year)) +
    geom_line(aes(y = get(sampleID), color = 'Tree'), na.rm = TRUE) +
    geom_point(aes(y = get(sampleID), color = 'Tree'), shape = 21, size = 2, fill = "white", na.rm = TRUE) +
    geom_line(aes(y = get(names(data)[2]), color = 'Master'), na.rm = TRUE) +
    geom_point(aes(y = get(names(data)[2]), color = 'Master'), shape = 21, size = 2, fill = "white", na.rm = TRUE) +
    labs(
      # title = paste(label," ", str_sub(sampleID, 3)),
      # title = bquote(bold(.(a_d)) ~ .(label) ~ .(str_sub(sampleID, 3))),
      #
         x = "Year", y = paste0(label, " rw (mm)"),
         color = "Series") +
    theme_minimal()
  return(p)
}


#' @keywords internal
create_barplot <- function(icol, data) {
  sampleID.chr <- names(data)[icol]
  if (str_detect(deparse(substitute(data)), "trt")) {
    parts <- str_split_fixed(names(data)[icol], "\\$",3)
    test <- parts[, 2]
    lagmax <- as.integer(parts[, 3])
    sampleID.o <- parts[, 1]

    label <- "treated"
    dt.clrs <- data.table(lag = -10:10)
    dt.clrs[, colr:= ifelse(lag == lagmax, "red", "black")]
    dt.clrs[lag == 0 & colr == 'black', colr:= "blue"]
    clrs <- setNames(as.character(dt.clrs$colr), as.character(dt.clrs$lag))
  } else {

    label <- "raw"
    test <- ""
    sampleID.o <- str_split_fixed(sampleID.chr, "\\_",2)[,2]
    dt.clrs <- data.table(lag = -10:10)
    dt.clrs[, colr:=  "black"]
    dt.clrs[lag == 0 , colr:= "blue"]
    clrs <- setNames(as.character(dt.clrs$colr), as.character(dt.clrs$lag))
  }


  # Create the bar plot
  p <- ggplot(data, aes(x = as.factor(lag), y = get(sampleID.chr), fill = as.factor(lag))) +
    geom_bar(stat = "identity",na.rm = TRUE) +
    scale_fill_manual(values = clrs) +
    # scale_fill_manual(values = clrs) +
    labs(
      # title = bquote(bold(.(a_d)) ~ "ccf_" * .(label) ~ .(sampleID.o) ~ .(test)),
         # title = bquote(bold(.(var)) ~ "ccf" * .(label))
      title = ifelse(label == "raw", "", paste( " qa_code: ", test)),
         x = "lag", y = paste0("correlation with ", label, " master"),
         color = "Series") +
    ylim(-1, 1) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}

# Helper function column-wise correlation with mean
cor_avg <- function(data){
  # Helper function using cor.test
  cor_with_stats <- function(x, y) {
    valid <- complete.cases(x, y)  # Remove NA values
    x <- x[valid]
    y <- y[valid]

    if (length(x) > 2) {
      test <- cor.test(x, y)
      list(corr = test$estimate, n = length(x), p_value = test$p.value)
    } else {
      list(corr = NA, n = length(x), p_value = NA)  # Insufficient data
    }
  }
  label <- tail(str_split(deparse(substitute(data)), "\\.")[[1]], 1)



  corr_stats <- lapply(3:ncol(data), function(i) {
    cor_with_stats(data$mean.rw, data[[i]])
  })

  # Create data.table with the results
  correlations_df <- data.table(
    SampleID.chr = names(data)[3:ncol(data)],
    corr_mean = round(sapply(corr_stats, function(x) x$corr),2),
    ncorr_mean = sapply(corr_stats, function(x) x$n),
    pcorr_mean = round(sapply(corr_stats, function(x) x$p_value),2)
  )
  names(correlations_df)[-1] <- paste0(names(correlations_df)[-1], "_", label)
  return(correlations_df)

}



