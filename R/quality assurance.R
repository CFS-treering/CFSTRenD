library(purrr)
library(data.table)
# step 1: pair-wise ccf on all the samples. find all pairs that satisfies the condition that max_ccf @lag0
#   the results of this step serves as the initial sample list of master chronology for step 2
# step 2: run ccf between each sample and the mean of initial master chronology
#   this step is to check ccf between each sample and the the mean master chronology, 
#   the results is the final master chronology satisfying the condition that max_ccf @lag0 with each sample that form the master
#   the out table contains the correlation with master at lag0 and its rank(from 1 to 21, 1 represents the maximum one, 21 represents minimum), 
#   test(Pass, borderline, highpeak, pm1) for each sample


# functions

# Define a function to calculate max cross-correlation lag
ccf_pairs <- function(ts1, ts2, lag.max = 10) {
  ccf.chk <- ccf(ts1, ts2, lag.max = lag.max, na.action = na.pass,  plot = FALSE)
  max_ccf <- max(ccf.chk$acf)
  max_lag <- ccf.chk$lag[which.max(ccf.chk$acf)]
  return(list(max_lag = max_lag, max_ccf = max_ccf))
}

ccf_avg <- function(icol, data, blcrit = 0.1, test="Fail"){
  
  dt.pairs <- ccf(data[,2],data[,icol, with = FALSE],lag.max=10,plot=FALSE, na.action = na.pass)
  dt.ccf <- data.table(SampleID.chr = names(data)[icol], lag = as.vector(dt.pairs$lag), acf.trt = as.vector(dt.pairs$acf))
  setorder(dt.ccf, -acf.trt)
  dt.ccf[, ccf.ord:= 1:.N]
  setorder(dt.ccf, lag)
  # max acf at lag 0
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 1 ) test="Pass"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 2]) == 1 & abs(dt.ccf[ccf.ord == 1]$acf.trt - dt.ccf[ccf.ord == 2]$acf.trt) <blcrit) test="borderline"
  if (nrow(dt.ccf[lag == 0 & ccf.ord == 1]) == 0 & dt.ccf[ccf.ord == 1]$acf.trt / dt.ccf[ccf.ord == 2]$acf.trt > 2) test="highpeak"
  if (nrow(dt.ccf[lag == 1 & ccf.ord == 1]) + nrow(dt.ccf[lag == -1 & ccf.ord == 1]) == 1) test="pm1"
  dt.ccf$test <- test
  return(dt.ccf)
}


find.master <- function(dt.trt_wide , dt.rw_long , max.trial = 100){
  # dt.trt_wide: treated series in wide format for pair-wise ccf
  # dt.rw_long: ring width series in long format for calculating mean of master chronology
  
  # step 1: pair-wise ccf to find all the samples which can find at least 1 sample to reach max_ccf @ lag0
  
  # Generate all pairs of columns
  col_pairs <- combn(names(dt.trt_wide), 2, simplify = FALSE)

  # Use map to calculate pairwise cross-correlation for each pair
  # this is the most efficient way for this calculation so far i found, 20-06-24, 2 mins for 490 series
  system.time(ccf.pairs <- map(col_pairs, ~ ccf_pairs(dt.trt_wide[[.x[1]]], dt.trt_wide[[.x[2]]])))
  
  # Name the ccf.pairs
  names(ccf.pairs) <- sapply(col_pairs, paste, collapse = "*")
  print(ccf.pairs)
  
  # Convert the list of ccf.pairs to a data.table
  dt.ccf.pairs <- rbindlist(lapply(names(ccf.pairs), function(pair) {
    res <- ccf.pairs[[pair]]
    data.table(pair = pair, max_lag = res$max_lag, max_ccf = res$max_ccf)
  }))
  
  dt.ccf.pairs[, c("ts1", "ts2"):=tstrsplit(pair, "*", fixed = TRUE)]
  
  result_dt.sel <- dt.ccf.pairs[max_lag == 0 & !is.na(max_ccf)]
  # summary(result_dt.sel$max_ccf)
  # quantile(result_dt.sel$max_ccf, 0.8, na.rm = TRUE)
  # result_dt.sel <- dt.ccf.pairs[max_ccf >= 0.5]
  # summary(result_dt.sel$max_ccf)
  ts.sel <- data.table(ts.sel = unlist(result_dt.sel$ts1, result_dt.sel$ts2))
  ts.sel <- ts.sel[, .N, by = .(ts.sel)]
  
  
  # the result of step 1 is id.sel, it serves as the initial sample list of master chronology for step 2
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  
  # step 2: find the sample list of master chronology satisfying the condition that max_ccf @lag0 for each sample in this list with the master chronology
  
  # algorithm on pass
  
  s2.end <- FALSE; id.sel <- unique(ts.sel$ts.sel); i.trial <- 0;
  while(!s2.end & i.trial < max.trial){
    
    # mean of master chronology
    dt.s2.avg <- dt.rw_long[SampleID.chr %in%  id.sel][, .(.N, mean.rw = mean(RawRing)), by = .(Year)]
    setorder(dt.s2.avg, Year)
    dt.s2.avg [, mean.rw.dif:= mean.rw - shift(mean.rw)]
    
    dt.s2.wide <- merge(dt.s2.avg[, c("Year", "mean.rw.dif")], dt.trt.wide, by = "Year")
    # ccf of all samples with the master chronology
    ccf.s2 <- lapply(3:ncol(dt.s2.wide), ccf_avg, data = dt.s2.wide)
    dt.ccf.s2 <-rbindlist(ccf.s2)
    # valid samples for master chronology
    id.sel2 <- unique(dt.ccf.s2[test == "Pass"]$SampleID.chr)
    str_ <- setdiff(union(id.sel, id.sel2) ,intersect (id.sel, id.sel2))
    
    # selection list stops changing
    if (length(str_) == 0) s2.end <- TRUE else{
      id.sel <- id.sel2
      i.trial <- i.trial + 1
      print(paste0(i.trial, " N.pass: ", length(id.sel)))
    } 
  }
  setorder(dt.ccf.s2, SampleID.chr)
  return(list(dt.ccf = dt.ccf.s2, master = dt.s2.avg ))
  # the result of step 2 is dt.ccf.s2, samples with test = "Pass" to form the master chronology
  
}



# step 0: read data, column names

# sourcedir="C:/Juha/TreeRingProjects/TR_QA_Example/"
dir.src <- "P:/Jing/2010-08/Martin/Treering_bank/TR_QA_Example_Juha/"

TreeData=fread(paste(dir.src,"RW_Hember.csv",sep=""))

# rename the column containing the data from "TRW" to "RawRing"
setnames(TreeData, "TRW", "RawRing")
#
# note that there were some duplicate "New_ID" so I had to add ID_Inst to New_ID to create a new unique ID that is called "SampleID"
TreeData[, SampleID := paste(ID_Inst,New_ID,sep="_")]
#
# Make a plot of all the raw data to see what it looks like
str(TreeData)
# ggplot(TreeData,aes(x=Year,y=RawRing))+geom_line(show.legend=FALSE)
# ggplot(TreeData,aes(x=Year,y=RawRing,group=SampleID,color=factor(SampleID)))+geom_line(show.legend=FALSE)
# +scale_color_viridis(discrete=TRUE)
#

# prepare dt.rw_long for the master chronology and dt.wide to be used for pair-wise ccf

dt.rw <- copy(TreeData)
setDT(dt.rw)

setorder(dt.rw, SampleID,Year)
dt.rw[, .N, by = .(SampleID, Year)][N > 1]
# the series to be used for ccf
dt.rw[, rw.treated:= RawRing - shift(RawRing), by = SampleID]

# SampleID.chr is the key sampleID, we use it in the functions to represent a sample.
# starting with a character to ensure the validity as column name and value of a column

dt.rw[, SampleID.chr:= paste0("d_", SampleID)]
setorder(dt.rw, SampleID.chr)
sample.lst <- sort(unique(dt.rw$SampleID.chr))

dt.trt.wide <- dcast(dt.rw[!is.na(rw.treated)], Year ~ SampleID.chr, value.var = "rw.treated")
names(dt.trt.wide)
setcolorder(dt.trt.wide, c("Year", sample.lst))

dt.trt_wide <- copy(dt.trt.wide)
# remove the column Year to faciliate pair-wise ccf
row.names(dt.trt_wide) <- dt.trt_wide$Year
dt.trt_wide$Year <- NULL
names(dt.trt_wide)
# data preparation ends here

# run the algorithm

master.trt.ccf <-find.master(dt.trt_wide = dt.trt_wide, dt.rw_long = dt.rw)

dt.trt.ccf <- master.trt.ccf$dt.ccf; masterchron <- master.trt.ccf$master; 





# create plot



# function on plotting
library(stringr)
library(ggplot2)

create_plot.series <- function(icol, data) {
  if (str_detect(deparse(substitute(data)), "trt")) label <- "treated" else label <- "raw"
   sampleID<- names(data)[icol]
  p <- ggplot(data, aes(x = Year)) +
    geom_line(aes(y = get(sampleID), color = 'Tree')) + 
    geom_point(aes(y = get(sampleID), color = 'Tree'), shape = 21, size = 2, fill = "white") +
    geom_line(aes(y = get(names(data)[2]), color = 'Master')) +
    geom_point(aes(y = get(names(data)[2]), color = 'Master'), shape = 21, size = 2, fill = "white") +
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
    geom_bar(stat = "identity") +
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




# pre data for plotting series vs year

dt.raw.series <- merge(masterchron[, c("Year", "mean.rw")], dt.o.wide, by = "Year")
dt.trt.series <- merge(masterchron[, c("Year", "mean.rw.dif")], dt.trt.wide, by = "Year")
names(dt.raw.series)
setcolorder(dt.raw.series, c("Year", "mean.rw", sample.lst))

names(dt.trt.series)
setcolorder(dt.trt.series, c("Year", "mean.rw.dif", sample.lst))


# pre data for bar plotting on ccf with master for raw and treated series
names(dt.trt.ccf)
dt.ccf.idtest <- dt.trt.ccf[ccf.ord==1, c("SampleID.chr", "test", "lag")]
dt.ccf.idtest[, idtest:= paste0(SampleID.chr,"$", test,"$", lag)]
dt.trt.ccf<- merge(dt.trt.ccf, dt.ccf.idtest[, c("SampleID.chr", "idtest")],by = "SampleID.chr")
dt.trt.ccf <- dcast(dt.trt.ccf, lag ~ idtest, value.var = "acf.trt")
names(dt.trt.ccf)
idtest.lst <- sort(unique(dt.ccf.idtest$idtest))
# test if in the same order as others
idtest.lst2 <- str_split_fixed(idtest.lst, "\\$",3)[,1]
if (!all.equal(idtest.lst2, sample.lst)) print("check the order of idtest in dt.ccf.idtest")
setcolorder(dt.trt.ccf, c("lag", idtest.lst))

# input data structure for ccf_avg Year, mean.value, sampleIDs...
ccf.raw <- lapply(3:ncol(dt.raw.series), ccf_avg, data = dt.raw.series)
dt.raw.ccf <- rbindlist(ccf.raw)
dt.raw.ccf <- dcast(dt.raw.ccf, lag ~ SampleID.chr, value.var = "acf.trt")
names(dt.raw.ccf)
setcolorder(dt.raw.ccf, c("lag", sample.lst))



# ccf of all samples with the master chronology
plot.trt.series <- lapply(3:ncol(dt.trt.series), create_plot.series, data = dt.trt.series)
plot.raw.series <- lapply(3:ncol(dt.raw.series), create_plot.series, data = dt.raw.series)

print(plot.trt.series[[2]])
print(plot.raw.series[[2]])
is.data.table(master.trt)
str(master.trt)
head(dt.ccf)


plot.trt.ccf <- lapply(2:ncol(dt.trt.ccf), create_barplot, data = dt.trt.ccf)
plot.raw.ccf <- lapply(2:ncol(dt.raw.ccf), create_barplot, data = dt.raw.ccf)

print(plot.raw.ccf[[2]])
icol <- 3; data <- dt.trt.ccf



library(ggplot2)

library(ggplot2)
library(gridExtra)

arrange_plots_to_pages <- function(plot_list1, plot_list2, plot_list3, plot_list4, output_file) {
  # Open a PDF device
  pdf(output_file, onefile = TRUE)
  
  # Loop through the plot lists and arrange them on each page
  for (i in 1:50) {
    plots <- list(plot_list1[[i]], plot_list2[[i]], plot_list3[[i]], plot_list4[[i]])
    grid.arrange(grobs = plots, ncol = 2, nrow = 2)
  }
  
  # Close the PDF device
  dev.off()
}

# Example usage
# Replace `plot_list1`, `plot_list2`, `plot_list3`, `plot_list4` with your actual plot lists
# arrange_plots_to_pages(plot_list1, plot_list2, plot_list3, plot_list4, "combined_plots.pdf")
arrange_plots_to_pages(plot.raw.series, plot.trt.series, plot.raw.ccf, plot.trt.ccf, output_file = "c:/tmp/ccf.pdf")










# display the results

library(gridExtra)
library(ggplot2)
i <-200
plot.raw.ccf[[2]]
plots <- list(plot.raw.series[[i]], plot.trt.series[[i]], plot.raw.ccf[[i]], plot.trt.ccf[[i]])
grid.arrange(grobs = plots, ncol = 2, nrow = 2)






arrange_plots_to_pages <- function(plot_list1, plot_list2, plot_list3, plot_list4) {
  for (i in 1:10) {
    plots <- list(plot_list1[[i]], plot_list2[[i]], plot_list3[[i]], plot_list4[[i]])
    grid.arrange(grobs = plots, ncol = 2, nrow = 2)
    # Save each page as a PDF
    ggsave(paste("page", i, ".pdf", sep = ""), plot = grid.arrange(grobs = plots, ncol = 2, nrow = 2))
  }
}

# Call the function to create and save the pages
arrange_plots_to_pages(plot_list1, plot_list2, plot_list3, plot_list4)



