# functions

#' plotting quality assessment
#' @description
#' plotting ccf results for each tree
#'
#' @param plot.lst plot list from generate_plots
#' @param out.series SampleID list to be plotted, default -999 for output all samples
#' @param out.pdf output file name including path/filename.pdf, default NULL for device screen



#' @import data.table
#' @import gridExtra
#' @importFrom grDevices dev.off pdf

#' @export plots_qa

plots_qa <- function(plot.lst, out.series = "all", out.pdf= NULL ) {



  samples.all <- names(plot.lst$plot.raw.series)

  if (all(tolower(out.series) == "all")) idx.lst <- seq_along(samples.all) else {
    idx.lst <- sapply(out.series, function(item) match(item, samples.all))
    idx.lst <- idx.lst[!is.na(idx.lst)]

  }


  if (all(is.na(idx.lst))) stop("cannot find out.series, please verify if they exist in sample.lst")
  # plot.lst <- generate.plots(master.trt, dt.input)


    if (!is.null(out.pdf))
    # Open a PDF device
    pdf(out.pdf, onefile = TRUE)


    # Loop through the plot lists and arrange them on each page
    for (i in idx.lst) {
      plots <- list(plot.lst$plot.raw.series[[i]], plot.lst$plot.trt.series [[i]], plot.lst$plot.raw.ccf[[i]], plot.lst$plot.trt.ccf[[i]])
      grid.arrange(grobs = plots, ncol = 2, nrow = 2)
    }

  # Close the PDF device
  if (!is.null(out.pdf)) dev.off()
}







#' plot frequency distribution by geo-location per species
#' @description
#' This function plots the frequency distribution by geo-location per species
#'
#' @param dt.freq a table resulting from function CFS_freq()
#' @param out.species species list, default is 'all' to output all species
#' @param out.pdf output file name including path/filename.pdf, default is NULL for screen display



#' @import data.table
#' @importFrom grDevices dev.off pdf

#' @export plots_freq

plots_freq <- function(dt.freq, out.species = "all", out.pdf= NULL ) {
  uid_label <- str_split(unique(dt.freq$uid_label), "_yr",2)[[1]]
  print(uid_label)
  dist.uids <- melt(dt.freq, id.vars = c("uid_label", "species" , "ord", "N", "pct.species" , "lat"),
                  variable.name = "lon",
                  value.name = "nuids")[!is.na(nuids)]
  dist.uids[, lon:= as.numeric(as.character(lon))]
  setorder(dist.uids, ord, lon, lat)
  if (all(tolower(out.species) == "all")) data.tmp <- dist.uids
  else{

    data.tmp <- dist.uids[species %in% out.species]
  }





  if (nrow(data.tmp) == 0) stop("cannot find out.series, please verify if they exist in sample.lst")
  # plot.lst <- generate.plots(master.trt, dt.input)



  data.tmp[, spc.pct := paste0(species, " N:", N, ", ", pct.species, "% (", ord, ")")]
  data.tmp$spc.pct <- factor(data.tmp$spc.pct, levels = unique(data.tmp$spc.pct[order(data.tmp$ord)]))
  # Get the unique spc.pct values
  spc_pcts <- unique(data.tmp$spc.pct)

  # Split into chunks of 4
  chunks <- split(spc_pcts, ceiling(seq_along(spc_pcts) / 2))


  if (!is.null(out.pdf))
    # Open a PDF device
    pdf(out.pdf, onefile = TRUE)

  for (i in seq_along(chunks)) {
    # Subset data for the current chunk
    data.i <- data.tmp[data.tmp$spc.pct %in% chunks[[i]], ]

    p1 <- ggplot(data.i, aes(x = lon, y = lat, size = nuids)) + facet_wrap(~spc.pct, ncol = 1, nrow = 2) +
      geom_point(alpha = 0.6, color = "darkblue") +
      scale_size_continuous(range = c(1, 10)) + # Adjust size range as needed
      scale_x_continuous(breaks = sort(unique(data.tmp$lon))) + # Set x-axis ticks to unique values of 'lat'

      theme_minimal() +
      theme(strip.text = element_text(size = 16),# Increase the size of facet labels
            panel.grid.minor = element_blank() , # Remove minor grid lines
            plot.title = element_text(size = 25), # Set title size
            plot.margin = margin(t = 10, r = 10, b = 30, l = 30, unit = "pt"), # Adjust plot margins
            plot.caption = element_text(hjust = 0, face = "italic")) + # Customize caption appearance

      labs(title = paste0( str_sub(uid_label[1], 5), " distribution by year ",uid_label[2]),
           x = "Longitude",
           y = "Latitude",
           caption = "Data source: CFS-TRenD V1.2 (mar-2024)" , # Add the data source caption
           size = paste0("n.", str_sub(uid_label[1], 5), "s"))

    print(p1)

  }

  # Close the PDF device
  if (!is.null(out.pdf)) dev.off()
}

