#' plotting quality assessment
#' @description
#' plotting ccf results for each tree
#'
#' @param plot.lst plot list from generate_plots
#' @param out.series SampleID list to be plotted, default -999 for output all samples
#' @param caption caption of the plots
#' @param out.pdf output file name including path/filename.pdf, default NULL for device screen



#' @import data.table
#' @import gridExtra
#' @import patchwork
#' @importFrom grDevices dev.off pdf

#' @export plots_qa

plots_qa <- function(plot.lst, out.series = "all", caption = "CFS-TRenD V1.2 samples69", out.pdf= NULL ) {

  if (nchar(caption) == 0) stop("please specify caption")

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
      # plots <- list(plot.lst$plot.raw.series[[i]], plot.lst$plot.trt.series [[i]], plot.lst$plot.raw.ccf[[i]], plot.lst$plot.trt.ccf[[i]])
      # grid.arrange(grobs = plots, ncol = 2, nrow = 2)
    p.i<-  (plot.lst$plot.raw.series[[i]] | plot.lst$plot.trt.series [[i]]) / (plot.lst$plot.raw.ccf[[i]] | plot.lst$plot.trt.ccf[[i]]) +
        plot_annotation(
          title = paste0("SampleID: ", samples.all[idx.lst[1]]) ,
          caption = paste0("Data source: ", caption),
          tag_levels = 'a',
          tag_suffix = ")") &
        theme(
          plot.title = element_text(face = "bold"),
          plot.tag = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0.2, face = "italic" )

          # plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt") # Adjust plot margins
          )
    print(p.i)
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
#' @param caption caption of data source
#' @param out.pdf output file name including path/filename.pdf, default is NULL for screen display



#' @import data.table
#' @importFrom grDevices dev.off pdf

#' @export plots_freq

plots_freq <- function(dt.freq, out.species = "all", caption = "", out.pdf= NULL ) {
  uid_label <- str_split(unique(dt.freq$uid_label), "_yr",2)[[1]]
  if (!is.na(uid_label[2])) title.tmp <- paste0( str_sub(uid_label[1], 5), " distribution by year ",uid_label[2]) else
    title.tmp <- paste0( str_sub(uid_label[1], 5), " distribution")
  # print(uid_label)
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

      labs(title = title.tmp,
           x = "Longitude",
           y = "Latitude",
           caption = paste0("Data source: ", caption) , # Add the data source caption
           size = paste0("n.", str_sub(uid_label[1], 5), "s"))

    print(p1)

  }

  # Close the PDF device
  if (!is.null(out.pdf)) dev.off()
}


#' plot data summary and location
#' @description
#' This function plots the site location and frequency distribution of series length and ring width measurement per species. It's used for generating the data report
#'
#' @param data a list of 2 tables : 1st is meta data with 1 species only ; 2nd is the ring width measurement in long-format


#' @import patchwork
#' @import ggplot2
#' @import sf

#' @export plots_ds

plots_ds <- function(data){
  dt.tr <- data[[1]]
  dt.rw <- data[[2]]
  spc <- unique(dt.tr$species)
  # if (nrow(dt.tr[!is.na(dbh_cm)]) > 0)
  #   p.dbh <- ggplot(dt.tr[!is.na(dbh_cm), .N, by = c("dbh_cm", "uid_tree", "species")], aes(x = dbh_cm)) +
  #   geom_histogram(binwidth = 5, color = "black", fill = "blue") +
  #   labs(title = "dbh (cm)", x = "dbh (cm)", y = "Frequency") +
  #   theme_classic() else p.dbh <- ggplot() + labs(title = "dbh frequency (no data)") +
  #   theme_minimal() +  # Use a minimal theme as a starting point
  #   theme(
  #     panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  #     plot.background = element_rect(fill = "transparent", color = NA)   # Transparent plot area
  #   )
  #
  # if (nrow(dt.tr[!is.na(ht_tot_m)]) > 0)
  #   p.ht <- ggplot(dt.tr[!is.na(ht_tot_m), .N, by = c("ht_tot_m", "uid_tree", "species")], aes(x = ht_tot_m)) +
  #   geom_histogram(binwidth = 5, color = "black", fill = "yellow") +
  #   labs(title = "height (m) ", x = "height (m)", y = "Frequency") +
  #   theme_classic() else p.ht <- ggplot() + labs(title = "height frequency (no data)") +
  #   theme_minimal() +  # Use a minimal theme as a starting point
  #   theme(
  #     panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  #     plot.background = element_rect(fill = "transparent", color = NA)   # Transparent plot area
  #   )

  # p.yr <- ggplot(dt.rw, aes(x = year)) +
  #   geom_bar(fill = "lightgrey", color = "lightblue") +
  #   labs(title = paste0("year frequency (max: ", nrow(tr_w6), ")"), x = "Year", y = "") +
  #   theme_minimal()
  #
  p.rw_hist <- ggplot(dt.rw, aes(x = rw_mm)) +
    geom_histogram(binwidth = 1, color = "black", fill =  "lightgreen") +
    labs(title = "ring width", x = "ring width (mm)", y = "Frequency") +
    theme_classic()
  # `map_bounds` <- st_bbox(shp)
  p.age <- ggplot(dt.rw, aes(x = rw_yend - rw_ystart + 1)) +
    geom_bar(fill = "lightblue", color = "lightblue") +
    labs(title = paste0("series length"), x = "series length", y = "Frequency") +
    theme_minimal()


    tmp <- dt.tr[, .N, by = .(longitude, latitude)]
  pts <- tmp %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84") %>%
    st_cast("POINT")

  p.loc <- ggplot() +
    geom_sf(data = shp, fill = NA, color = "lightgrey", alpha = 0) + # Polygons
    geom_sf(data = pts, color = "blue", size = 2) +          # Points
    labs(title = "site location ") +
    # annotate("text",
    #              x = map_bounds["xmin"] + 0.1 * (map_bounds["xmax"] - map_bounds["xmin"]), # Slightly inset from the left
    #   y = map_bounds["ymax"] + 0.2 * (map_bounds["ymax"] - map_bounds["ymin"]), # Slightly above the top
    #
    #          label = paste0("location"),
    #          hjust = 0, vjust = 1, color = "black", size = 7) +   # Credits
    coord_sf(expand = FALSE) +                                  # Add graticules with labels
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey", linetype = "dotted"), # Graticule styling
      panel.grid.minor = element_blank(),                                   # No minor lines
      axis.text = element_text(size = 7),                                 # Graticule label styling
      axis.ticks = element_line(size = 0.5),
      axis.title = element_blank()
    )

  p.rw <- ggplot(dt.rw, aes(x = year, y = rw_mm)) +
    geom_point(color = "lightgrey") +
    labs(title = "ring width(mm) ", y = "ring width (mm)", x = "year") +
    theme_classic()

  p.ispc<-  (p.loc | p.rw)/ (p.age | p.rw_hist) +
    plot_annotation(
      title = paste0("Species: ", spc) ,
      tag_levels = 'a',
      tag_suffix = ")") &
    theme(
      plot.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold")


      # plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt") # Adjust plot margins
    )
  return(p.ispc)
}


#' plot multiple series by wrap_facet, and put them in a list
#' @description
#' This function plots dynamic xyplot per category using wrap_facet. each facet contains nrow*ncol plots. It's used for generating the data report
#'
#' @param data a data table including all the necessary columns in varcols
#' @param varcols a list containing 3 column names for x, y and facet category.
#' @param xylabels a list containing the labels for x-axis and y-axis
#' @param nrow number of plots per row
#' @param ncol number of plots per col


#' @import ggforce
#' @import ggplot2


#' @export plots_facet

plots_facet <- function(data, varcols, xylabels, nrow, ncol) {
  # Calculate the total number of pages
  # print(paste0(nrow, " ", ncol))
  total_pages <- ceiling(length(unique(data[[varcols[[3]]]])) / (nrow*ncol))
  # print(total_pages)

  # Generate plots for each page
  plot_list <- lapply(1:total_pages, function(page) {
    ggplot(data, aes(x = .data[[varcols[[1]]]], y = .data[[varcols[[2]]]])) +
      geom_point(color = "lightgrey") +
      facet_wrap_paginate(
        ~.data[[varcols[[3]]]],
        ncol = ncol,
        nrow = nrow,
        page = page
      ) +
      scale_x_continuous(labels = scales::number_format(scale = 1)) +  # Ensure continuous year labels
      labs(x = xylabels[[1]], y = xylabels[[2]]) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  # 2 per page
  if ((floor(total_pages/2))*2 != total_pages) plot_list[[total_pages+1]] <- ggplot() +
    theme(
      panel.background = element_blank(),  # No background
      plot.background = element_blank(),   # No plot area background
      axis.ticks = element_blank(),        # No axis ticks
      axis.text = element_blank(),         # No axis text
      axis.title = element_blank()         # No axis titles
    )

  return(plot_list)
}


#' Generate a data/model Report
#'
#' This function generates an HTML data report using an R Markdown template.
#' @param robj an r object, either a list of data tables running from CFS_format()  or a list running from gamm_series, gamm_site, gamm_spatial or bam_spatial
#' @param usage for data report, 1 for data submission, hence the complete report; 2 for modelling use, simplified report. It has no effect on model report.
#' @param output_file A string specifying the name of the output HTML file.
#'
#' @return An HTML file containing the data report.
#' @export generate_report

generate_report <- function(robj, usage = 1, output_file = NULL) {
  # if (!input[[1]] %in% c("data", "model")) stop("the first item needs to be data or model")
  # robj <- input[[2]]
  if ( "data.frame" %in% class(robj$tr_all_wide) | "data.table" %in% class(robj$tr_all_wide)) type.templt <- "data" else {
    if ( "list" %in% class(robj) & !is.null(robj$model)) type.templt <- "model"
  }
  # Path to the R Markdown template
  rmd_file <- system.file("rmd", paste0(type.templt,"_report_template.Rmd"), package = "CFSTRenD")
  # rmd_file <- file.path("P:/Jing/2010-08/Martin/Treering_bank/Git/CFSTRenD/inst/rmd", paste0(type.templt,"_report_template.Rmd"))


                          # Check if the template exists
                          if (rmd_file == "") {
                            stop("Template not found! Ensure 'data_report_template.Rmd' is in inst/rmd/ directory.")
                          }

                          # Render the report
                          if (is.null(output_file) ) rstudioapi::viewer( rmarkdown::render(input = rmd_file,
                                                                                                      params = list(robj = robj, usage = usage))) else

                            rmarkdown::render(
                              input = rmd_file,
                              output_file = output_file,
                              params = list(robj = robj, usage = usage),
                              envir = new.env(parent = globalenv()) # Isolate environment
                            )
}


# generate_vignette <- function(output_file = NULL) {
#
#   # Path to the R Markdown template
#   # rmd_file <- system.file("rmd", paste0(type.templt,"_report_template.Rmd"), package = "CFSTRenD")
#   rmd_file <- file.path("P:/Jing/2010-08/Martin/Treering_bank/Git/CFSTRenD/inst", paste0("vignette_CFSTRenD.Rmd"))
#
#   # Check if the template exists
#   if (rmd_file == "") {
#     stop("Template not found! Ensure 'data_report_template.Rmd' is in inst/rmd/ directory.")
#   }
#
#   # Render the report
#   if (is.null(output_file) ) rstudioapi::viewer( rmarkdown::render(rmd_file)) else
#
#     rmarkdown::render(
#       input = rmd_file,
#       output_file = output_file,
#       envir = new.env(parent = globalenv()) # Isolate environment
#     )
# }
