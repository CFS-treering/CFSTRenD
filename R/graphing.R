# functions

#' plotting ccf
#' @description
#' plotting ccf results for each tree
#'
#' @param plot.lst : plot list from generate_plots
#' @param out_series : 3 options(A: all series; P: specific series of out_N; R: randomly selected samples)
#' @param out_N : for out_series A, not necessary; for P: list of index of specific series; R: number of series to be randomly selected
#' @param outpdf_YN :  option to output the results in pdf ("Y", "N")
#' @param outpdf_fn : output file name including path/filename.pdf



#' @import data.table
#' @import gridExtra
#' @importFrom grDevices dev.off pdf

#' @export plot_results

plot_results <- function(plot.lst, out_series = "A", out_N = NULL, outpdf_YN ,outpdf_fn ) {

  if (!(toupper(out_series) %in% c("A", "P", "R"))) stop("out_series has 3 options: 'A', 'P' or 'R'" )
  if (!(toupper(outpdf_YN) %in% c("Y", "N"))) stop("outpdf_YN has 3 options: 'Y' or 'N'" )

  idx.all <- seq_along(plot.lst$plot.raw.series)

  if ((out_series) == "A") idx.lst <- idx.all

  if ((out_series) == "R"){
    if (is.null(out_N)) stop("out_N cannot be null to output random-selected") else
      idx.lst <- sample(idx.all, out_N)

  }
  if ((out_series) == "P"){
    if (is.null(out_N)) stop("out_N cannot be null to output random-selected") else
      idx.lst <- out_N
  }


  # plot.lst <- generate.plots(master.trt, dt.input)

  if (toupper(outpdf_YN) == "Y")  {
    if (is.null(outpdf_fn)) stop("please input the path/filename as outpdf_fn" )
    # Open a PDF device
    pdf(outpdf_fn, onefile = TRUE)
  }

    # Loop through the plot lists and arrange them on each page
    for (i in idx.lst) {
      plots <- list(plot.lst$plot.raw.series[[i]], plot.lst$plot.trt.series [[i]], plot.lst$plot.raw.ccf[[i]], plot.lst$plot.trt.ccf[[i]])
      grid.arrange(grobs = plots, ncol = 2, nrow = 2)
    }

  # Close the PDF device
  if (toupper(outpdf_YN) == "Y") dev.off()
}

# pre data for plotting series vs year
# pre.plot <- function(masterchron, dt.rw.wide, dt.trt.wide, sample.lst){

#' generating plots for all series
#' @description
#' generating plots for all series
#' @param master.trt : master chronology from tr.fullmaster function
#' @param dt.input : tables from pre.dataFormat function


#' @import data.table
#' @import ggplot2

#' @export generate_plots

generate_plots <- function(master.trt, dt.input){

  # pre series data for both rw and treated in wide format
  dt.raw.series <- merge(master.trt$master[, c("Year", "mean.rw")], dt.input$dt.rw_wide, by = "Year")
  dt.trt.series <- merge(master.trt$master[, c("Year", "mean.rw.dif")], dt.input$dt.trt_wide, by = "Year")

  setcolorder(dt.raw.series, c("Year", "mean.rw", dt.input$sample.lst))
  setcolorder(dt.trt.series, c("Year", "mean.rw.dif", dt.input$sample.lst))


  # pre data for bar plotting on ccf with master for raw and treated series

  dt.trt.ccf <- master.trt$dt.ccf
  dt.ccf.idlabel <- dt.trt.ccf[ccf.ord==1, c("SampleID.chr", "qa_code", "lag")]
  dt.ccf.idlabel[, id.label:= paste0(SampleID.chr,"$", qa_code,"$", lag)]
  dt.trt.ccf<- merge(dt.trt.ccf, dt.ccf.idlabel[, c("SampleID.chr", "id.label")],by = "SampleID.chr")
  dt.trt.ccf <- dcast(dt.trt.ccf, lag ~ id.label, value.var = "acf.trt")
  names(dt.trt.ccf)
  idlabel.lst <- sort(unique(dt.ccf.idlabel$id.label))
  # test if in the same order as others
  idlabel.lst2 <- str_split_fixed(idlabel.lst, "\\$",3)[,1]
  if (!all.equal(idlabel.lst2, dt.input$sample.lst)) print("check the order of id.label in dt.ccf.idlabel")
  setcolorder(dt.trt.ccf, c("lag", idlabel.lst))

  # input data structure for ccf_avg Year, mean.value, sampleIDs...

  dt.raw.ccf <- rbindlist(lapply(3:ncol(dt.raw.series), ccf_avg, data = dt.raw.series))
  dt.raw.ccf <- dcast(dt.raw.ccf, lag ~ SampleID.chr, value.var = "acf.trt")
  names(dt.raw.ccf)
  setcolorder(dt.raw.ccf, c("lag", dt.input$sample.lst))



  # ccf of all samples with the master chronology
  plot.trt.series <- lapply(3:ncol(dt.trt.series), create_plot.series, data = dt.trt.series)
  plot.raw.series <- lapply(3:ncol(dt.raw.series), create_plot.series, data = dt.raw.series)

  # print(plot.trt.series[[2]])
  # print(plot.raw.series[[2]])



  plot.trt.ccf <- lapply(2:ncol(dt.trt.ccf), create_barplot, data = dt.trt.ccf)
  plot.raw.ccf <- lapply(2:ncol(dt.raw.ccf), create_barplot, data = dt.raw.ccf)

  # print(plot.raw.ccf[[2]])
  # print(plot.trt.ccf[[2]])
  return(list(plot.raw.series = plot.raw.series, plot.trt.series = plot.trt.series, plot.raw.ccf = plot.raw.ccf, plot.trt.ccf = plot.trt.ccf))
}

create_plot.series <- function(icol, data) {
  if (str_detect(deparse(substitute(data)), "trt")) label <- "treated" else label <- "raw"
  sampleID<- names(data)[icol]
  p <- ggplot(data, aes(x = Year)) +
    geom_line(aes(y = get(sampleID), color = 'Tree'), na.rm = TRUE) +
    geom_point(aes(y = get(sampleID), color = 'Tree'), shape = 21, size = 2, fill = "white", na.rm = TRUE) +
    geom_line(aes(y = get(names(data)[2]), color = 'Master'), na.rm = TRUE) +
    geom_point(aes(y = get(names(data)[2]), color = 'Master'), shape = 21, size = 2, fill = "white", na.rm = TRUE) +
    labs(title = paste(label," ", str_sub(sampleID, 3)),
         x = "Year", y = paste0(label, " rw"),
         color = "Series") +
    theme_minimal()
  return(p)
}



create_barplot <- function(icol, data) {
  sampleID <- names(data)[icol]
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
    sampleID.o <- sampleID
    dt.clrs <- data.table(lag = -10:10)
    dt.clrs[, colr:=  "black"]
    dt.clrs[lag == 0 , colr:= "blue"]
    clrs <- setNames(as.character(dt.clrs$colr), as.character(dt.clrs$lag))
  }


  # Create the bar plot
  p <- ggplot(data, aes(x = as.factor(lag), y = get(sampleID), fill = as.factor(lag))) +
    geom_bar(stat = "identity",na.rm = TRUE) +
    scale_fill_manual(values = clrs) +
    # scale_fill_manual(values = clrs) +
    labs(title = paste("ccf_", label, " ", str_sub(sampleID.o, 3), " ", test),
         x = "lag", y = paste0("correlation with ", label, " master"),
         color = "Series") +
    ylim(-1, 1) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}
