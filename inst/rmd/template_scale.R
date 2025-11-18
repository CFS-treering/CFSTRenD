params <-
list()

## ----setup_scale, include=FALSE-----------------------------------------------
# This chunk sets global options to hide all R code from appearing in the output
if (!requireNamespace("htmltools", quietly = TRUE)) {
  stop("Package 'htmltools' is required to run this function. Please install it.")
}
knitr::opts_chunk$set(echo = FALSE)
robj <- params$robj

## ----results='asis', warning=FALSE, message=FALSE, fig.width=10, fig.height=5, fig.align='left'----
library(data.table)
library(CFSTRenD)
library(gt)

library(stringr)
library(patchwork)
library(ggplot2)
# robj <- qa.trt.ccf; out.series.sel <- c("X011_104_003", "X016_104_003")
# robj <- params$robj; out.series.sel = params$parms.reports["qa.out_series"]
# print(params$parms["qa.out_series"])
# plt.lst <- plot_qa(params$robj, out.series = out.series.sel)

# plt.lst <- plot_scale(robj)
# 
# 
# 
#     # Loop through the plot lists and arrange them on each page
#    # for (i in 1:length(plt.lst$plot.year)) {
#   # plots <- list(plot.lst$plot.raw.series[[i]], plot.lst$plot.trt.series [[i]], plot.lst$plot.raw.ccf[[i]], plot.lst$plot.trt.ccf[[i]])
#   # grid.arrange(grobs = plots, ncol = 2, nrow = 2)
#   # if (length(plt.lst$plot.year) > 1) p.i<-  (plt.lst$plot.year[[i]] | plt.lst$plot.ll [[i]])  else
#     p.i<-  (plt.lst$plot.year | plt.lst$plot.ll) +
#     plot_annotation(
#       title = paste0("Site: ", robj$scale.parms["scale.label_data_site"]) ,
#       caption = paste0("Data source: ", robj$scale.parms["scale.label_data_ref"]),
#       tag_levels = 'a',
#       tag_suffix = ")") &
#     theme(
#       plot.title = element_text(face = "bold"),
#       plot.tag = element_text(face = "bold"),
#       plot.caption = element_text(hjust = 0.2, face = "italic" )
#       
#       # plot.margin = margin(t = 10, r = 10, b = 30, unit = "pt") # Adjust plot margins
#     )
#   print(p.i)

# }

# robj: a list of cfs_scale objects
# Each robj[[i]] can be used with plot_scale()

for (i in seq_along(robj)) {
  
  # Generate plots for this scale object
  plt.lst <- plot_scale(robj[[i]])
  
  # Combine the year and lat-long plots side-by-side
  p.i <- (plt.lst$plot.year | plt.lst$plot.ll) +
    plot_annotation(
      title = paste0(
    "Species: ", sub("_.*", "", names(robj)[[i]]), "\n",
    "Site: ", sub(".*_", "", names(robj)[[i]])
  ),
      caption = paste0("Data source: ", robj[[i]]$scale.parms["scale.label_data_ref"]),
      tag_levels = 'a',
      tag_suffix = ")"
    ) &
    theme(
      plot.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0.2, face = "italic")
    )
  
  print(p.i)
}


