% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{gamm_spatial}
\alias{gamm_spatial}
\title{spatial growth model at regional-level (multiple sites)}
\usage{
gamm_spatial(data, resp_scale = "resp_gamma", m.candidates)
}
\arguments{
\item{data}{data containing all necessary columns to run the model}

\item{resp_scale}{the scale of response variable. default is "log" for log-scale, otherwise modeling the response variable as it is}

\item{m.candidates}{the list of formulas.}
}
\value{
list including model, fitting statistics, ptable, stable, prediction table and spatial effect(moranI)
}
\description{
models the growth trend or climate-growth relationship at regional-level with multiple sites
}
\details{
This function accounts for within-site variability and temporal autocorrelation by including series identity as random effects
and a first-order autoregressive (AR1) correlation structures, respectively. Among-site variability and spatial effects are captured by incorporating site identity as random effects.
The model is refitted automatically by introducing a smooth term for latitude and longitude using the Spatial Over-Smooth ("sos") basis if significant spatial autocorrelation persists.
“Normalized” residuals are provided for future analysis.

If users specify multiple candidate models through the m.candidates argument, the function will fit each candidate model using the maximum likelihood (ML) method.
The Akaike Information Criterion (AIC) will then be compared to determine the best-fitting model.
Once the optimal model is identified, it will be refitted using the restricted maximum likelihood (REML) method and output the results.

If users specify only 1 candidate model through the m.candidates argument, the model is fitted with "REML" method.
}
