params <-
list()

## ----setup_model, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
library(data.table)
library(mgcv)
library(ggplot2)
library(dplyr)
library(stringr)
library(CFSTRenD)
library(ggeffects)
library(cowplot)


## ----warning=FALSE, message=FALSE---------------------------------------------

robj <- params$robj

## -----------------------------------------------------------------------------
if (is.list(robj$model) && all(c("gam", "lme") %in% names(robj$model))) gam_model <- robj$model$gam else gam_model <- robj$model
summary(gam_model)
term_important <- sterm_imp( gam_model)

method <- unique(term_important$method)
tab <- term_important[, c(1,2)]
colnames(tab) <- c("Term", "Score (%)")

## -----------------------------------------------------------------------------

knitr::kable(tab, caption = "Relative importance of smooth terms based on deviance partitioning.",
      align = c("l", "r"))

## ----echo=FALSE,   warning=FALSE, message=FALSE, fig.width=12, fig.height=8, fig.align='left'----
plot_list <- plot_resp(robj = robj)

combined <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  labels = LETTERS[1:length(plot_list)]
)

print(combined)


## ----echo=FALSE,   warning=FALSE, message=FALSE, fig.width=12, fig.height=8, fig.align='left'----

 # Check smoothness, residuals, K-index

 # 2x2 layout on one page
par(mfrow = c(2, 2))
  if (is.list(robj$model) && all(c("gam", "lme") %in% names(robj$model))) gam_model <- robj$model$gam else gam_model <- robj$model
gam.check(gam_model, rep = 100, verbose = TRUE)

# reset plotting layout to default
par(mfrow = c(1, 1))
 


