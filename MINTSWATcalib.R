library(parallel)

MINTSWATcalib <- function() {
  # Function that runs through the calibration setup.
  save(readSWAT, file = "readSWAT.R")
  change_params <- ""
  rm(change_params)
  load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
  print("DEoptim Calibration has the can currently calibrate all
         of the following parameters within the following files:\n")
  print(change_params[, 1:2])
  calib_range <- c("1999-12-31", "2021-12-31")
  params_select <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 21, 23, 24, 32, 33)
  calib_params <- change_params[params_select, ]
  print("You are calibrating the following parameters within the following files:\n")
  print(calib_params[, 1:2])

  calib_params[grep("Ksat", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.5, 1.5, 1)
  calib_params[grep("SMFMN", calib_params[, "parameter"]), c("min", "max", "current")] <- c(0, 5, 2.5)
  calib_params[grep("SMFMX", calib_params[, "parameter"]), c("min", "max", "current")] <- c(0, 5, 2.5)
  calib_params[grep("TIMP", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.01, 1, .5)
  calib_params[grep("CN2", calib_params[, "parameter"]), c("min", "max", "current")] <- c(35, 95, 70)
  calib_params[grep("Depth", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.5, 2, 1)
  calib_params[grep("Ave", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.5, 2, 1)
  calib_params[grep("ALPHA_BF", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.01, 1, .8)
  calib_params[grep("GWQMN", calib_params[, "parameter"]), c("min", "max", "current")] <- c(.1, 600, 1)
  calib_params[grep("GW_REVAP", calib_params[, "parameter"]), c("min", "max", "current")] <- c(0, .3, .02)

  setup_swatcal(calib_params)
  cl <- makeCluster(4)
  print(cl)
  clusterEvalQ(cl, c(source("SWATmodel")))
  # Test calibration
  x <- calib_params$current
  swat_objective_function_rch(x, calib_range, calib_params, flowgage, rch, save_results = F)
  outDEoptim <- DEoptim(
    swat_objective_function_rch, calib_params$min, calib_params$max,
    DEoptim.control(
      strategy = 6, NP = 16, itermax = deiter,
      cluster = cl
    ), calib_range, calib_params, flowgage, rch
  )
  x <- outDEoptim$optim$bestmem # need to save this, along with an ArcSWAT like directory structure for the basin
  swat_objective_function_rch(x, calib_range, calib_params, flowgage, rch, save_results = TRUE)
  calibrated_params <- calib_params
  calibrated_params$current <- x
  write.csv(calibrated_params, "CalibratedParams.csv", row.names = F)
}