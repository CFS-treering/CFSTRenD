% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CFSformat.R
\name{CFS_format}
\alias{CFS_format}
\title{convert data to CFS-TRenD format}
\usage{
CFS_format(data, out.csv = FALSE, out.dir = NULL)
}
\arguments{
\item{data}{: a list, first is input data in wide format; second is a flat sequence referring to the column indices of meta variables}

\item{out.csv}{: output csv file (TRUE or FALSE)}

\item{out.dir}{: directory of output csv files}
}
\value{
7 tables starting with tr_
}
\description{
convert data to CFS-TRenD format
}
