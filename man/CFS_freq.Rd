% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_uids.R
\name{CFS_freq}
\alias{CFS_freq}
\title{data spatial distribution by species measured later than certain year}
\usage{
CFS_freq(
  tr_meta,
  crit_yr,
  resol_lat,
  resol_lon,
  uid,
  spc.lst,
  N.species = 999,
  make.plot = TRUE,
  out.csv = FALSE,
  out.dir = NULL
)
}
\arguments{
\item{tr_meta}{:  meta table from CFS_format}

\item{crit_yr}{: year criteria}

\item{resol_lat}{:latitude resolution}

\item{resol_lon}{:longitude resolution}

\item{uid}{: counts on which uid}

\item{spc.lst}{species list}

\item{N.species}{: number of top species to be output (999 for all)}

\item{make.plot}{plot the distribution (TRUE or FALSE)}

\item{out.csv}{: output csv file (TRUE or FALSE)}

\item{out.dir}{: directory of output csv files}
}
\value{
1 table species-lat * lon, values are number of uids
}
\description{
data spatial distribution by species measured later than certain year
}
