% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{gamm_radius}
\alias{gamm_radius}
\title{detrending model on tree-ring width series}
\usage{
gamm_radius(data, resp_scale = "resp_gamma", m.candidates)
}
\arguments{
\item{data}{data containing all necessary columns to run the model}

\item{resp_scale}{the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is}

\item{m.candidates}{the list of formulas.}
}
\value{
list including model, fitting statistics, ptable, stable and prediction table
}
\description{
models the biological growth trends in individual tree-ring width series using mgcv::gamm
}
\details{
This function models the biological growth trends in individual tree-ring width series using mcgv::gamm.
By integrating a first-order autoregressive (AR1) component, it accounts for temporal autocorrelation.
This method can provide  “normalized” residuals, which are adjusted to reflect deviations after considering the AR1 correlation structure.
'Normalized' residuals are valuable for further analyses, such as investigating relationships with climatic variables.
If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model. Once the optimal model is identified,
it will be refitted using the restricted maximum likelihood (REML) method and output the results.

If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.
}
