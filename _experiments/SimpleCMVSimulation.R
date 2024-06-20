library(tidyverse)
library(furrr)
library(lavaan)
library(simstandard)
library(semPlot)

set.seed(1337)
plan("multicore") # when using parallel processing in furrr, use all cores

source("cmv_sim_functions.R", local = knitr::knit_global())

cfa_model <- '
  # measurement model
  X =~ x1 + x2 + x3 + x4 + x5
  Y =~ y1 + y2 + y3 + y4 + y5
  Z =~ z1 + z2 + z3 + z4 + z5
  
  # covariances
  X ~~ Y
  Y ~~ Z
  X ~~ Z
'

semPlotModel_lavaanModel(cfa_model) |> 
  semPaths(layout="tree", rotation=2, nCharNodes=5, sizeMan2 = 2.5, mar = c(1,4,1,3), label.norm = "OOOOO", residuals = FALSE)

cfa_model_pars <- lavaanify(cfa_model, std.lv = TRUE) |> add_cmv_marker()

semPlotModel_lavaanModel(cfa_model_pars) |> 
  semPaths(layout="tree", rotation=2, nCharNodes=5, sizeMan2 = 2.5, mar = c(1,4,1,3), label.norm = "OOOOO", residuals = FALSE)

loadings <- 0.7
cmv_loadings = seq(0, 0.2, by = 0.1)
covariances <- 0 
cfa_model_pars <- cfa_model_pars |> fix_cmv_loadings(loadings, invert = TRUE) |> fix_covariances(covariances)
cfa_model_pars_list <- map(cmv_loadings, \(x) fix_cmv_loadings(cfa_model_pars, x)) # make a list of data frames with each frame having one of the CMV loadings
fits <- map(cfa_model_pars_list, \(x) { #TODO: change to future_map
  #cfa_sim <- sim_standardized(x, latent = FALSE, errors = FALSE, factor_scores = FALSE,  composites = FALSE, matrices = FALSE) |> 
  cfa_sim <- simulateData(x, model.type = "cfa", standardized = TRUE) |>
    # create likert scores using standard z scores for quintiles
    mutate(across(everything(), \(y) findInterval(y, vec=c(-Inf, -1.2,-0.4, 0.4, 1.2,Inf)))) 
  print(lavaanify(fixed2free(cfa_model_pars)))
  return(cfa(fixed2free(cfa_model_pars), cfa_sim, std.lv=TRUE))
}) 

test <- do.call(lavTestLRT, fits)
tibble(test)

fit_measures <- map(fits, \(x) fitmeasures(x))  |>
  # convert from list of measures to tibble
  map_dfr(~tibble(fit_measure = names(.x), value = .x), .id = "run") 