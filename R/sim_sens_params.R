#' Simulation-based sensitivity analysis for Vmax and Km
#'
#' This function performs sensitivity analysis by testing a series of Vmax and Km values.
#' @param dat A data frame containing the pharmacokinetic data with columns for EVID and DV.
#' @param sim_vmax_list A list of Vmax values to simulate (optional).
#' @param sim_km_list A list of Km values to simulate (optional).
#' @param estvc Estimated volume of distribution.
#' @param estcl Estimated clearance (optional, required if `sim_vmax_list` and `sim_km_list` are not provided).
#' @param estcmax Estimated maximum concentration (optional, required if `sim_vmax_list` and `sim_km_list` are not provided).
#' @return A data frame containing the Vmax, Km, APE, MAPE, and time spent for each simulation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPT
#' sim_sens_vmax_km(dat = Bolus_1CPT,sim_vmax_list = c(500,1000,1500),sim_km_list = c(200,300,400),estvc = 70)
#' @export

sim_sens_vmax_km <- function(dat,
                             sim_vmax_list = NA,
                             sim_km_list = NA,
                             estvd,
                             estcl,
                             estcmax,
                             estka = NA,
                             noniv_flag = 0) {
  start.time <- Sys.time()

  if (missing(estvd)) {
    stop("Error: no volume of distribution for simulation was provided")
  }

  if (missing(estcl)) {
    stop("Error: no estimated clearance for simulation was provided")
  }

  if (missing(estcmax)) {
    stop("Error, no estimated cmax for simuation was provided")
  }

  combs_df <- data.frame(vmax = sim_vmax_list, km = sim_km_list)

  if (is.na(sim_vmax_list) == T || is.na(sim_km_list) == T) {
    # generate the default vman and km test list from linear kinetcs to nonlinear kinetics
    km_range <- c(4, 2, 1, 0.5, 0.25, 0.125, 0.1, 0.05) * estcmax
    # cl ~ range (0, vmax/km)
    # assume several concentration values.
    conc_range <- c(0.1, 0.25, 0.5, 0.75) * estcmax
    combs <- expand.grid(km_range, conc_range)
    combs_df <- data.frame(combs)
    colnames(combs_df) <- c("km", "conc")
    combs_df$vmax <- (combs_df$km + combs_df$conc) * estcl

  }

  # start simulation
  sim.vmax.km.results.all <- NULL

  for (loopmm in 1:nrow(combs_df)) {
    input.vmax <- combs_df[loopmm, ]$vmax
    input.km <- combs_df[loopmm, ]$km

    if (noniv_flag == 0) {
      sim.lst <- Fit_1cmpt_mm_iv(
        data = dat[dat$EVID != 2, ],
        est.method = "rxSolve",
        input.vmax = input.vmax,
        input.km = input.km,
        input.vd = estvd,
        input.add = 0
      )
    }

    if (noniv_flag == 1) {
      sim.lst <- Fit_1cmpt_mm_oral(
        data = dat[dat$EVID != 2, ],
        est.method = "rxSolve",
        input.ka = input.ka,
        input.vmax = input.vmax,
        input.km = input.km,
        input.vd = estvd,
        input.add = 0
      )
    }

    sim.APE <-
      round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[1], 1)
    sim.MAE <-
      round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[2], 1)
    sim.MAPE <-
      round(metrics.(pred.x = sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[3], 1)
    sim.RMSE <-
      round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[4], 1)
    sim.rRMSE <-
      round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[5], 1)

    end.time <- Sys.time()
    time.spent <- round(difftime(end.time, start.time), 4)

    sim.vmax.km.results <-
      data.frame(
        vmax = input.vmax,
        km = input.km,
        sim.mm.APE = sim.APE,
        sim.mm.MAE = sim.MAE,
        sim.mm.MAPE = sim.MAPE,
        sim.mm.RMSE = sim.RMSE,
        sim.mm.rRMSE = sim.rRMSE,
        time.spent = time.spent
      )
    sim.vmax.km.results.all <-
      rbind(sim.vmax.km.results.all, sim.vmax.km.results)

    rm(sim.lst) # To avoid occupying memory
    gc()
  }
  return(sim.vmax.km.results.all)
}


#' Simulation-based sensitivity analysis for  for a two-compartment model
#'
#' Performs sensitivity analysis by testing a series of potential ratios of vc to vp.
#' @param dat A data frame containing the pharmacokinetic data
#' @param sim_vc_list A list of central compartment volumes (Vc) to simulate (optional).
#' @param sim_vp_list A list of peripheral compartment volumes (Vp) to simulate (optional).
#' @param sim_q_list A list of inter-compartmental clearance (Q) to simulate (optional).
#' @param estka Estimated absorption rate if oral case.
#' @param estcl Estimated clearance.
#' @param estvc Estimated central volume of distribution (optional, required if `sim_vc_list` and `sim_vp_list` are not provided).
#' @return A data frame containing the Vc, Vp, APE, MAPE, and time spent for each simulation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' sim.results<-sim_sens_2cmpt(dat, estcl = 4, sim_vc_list = 70)
#' @export
#'
sim_sens_2cmpt <- function(dat,
                           sim_vc_list,
                           sim_vp_list = NA,
                           sim_q_list = NA,
                           estcl,
                           estka = NA,
                           noniv_flag = 0) {
  if (missing(estcl)) {
    stop("Error: no estimated clearance for simulation was provided")
  }

  if (missing(sim_vc_list)) {
    stop("Error: vc provided for simulation analysis")
  }

  sim_vc_list <- na.omit(sim_vc_list)
  sim_vc_list <- sort(sim_vc_list)

  # Check if the two numbers are within 20% relative difference
  if (length(sim_vc_list) == 2) {
    relative_difference <-
      abs(sim_vc_list[2] - sim_vc_list[1]) / sim_vc_list[1]
    if (relative_difference <= 0.2) {
      # If they are within 20%, keep the smaller number
      deduplicated_values <- sim_vc_list[1]
    } else {
      # If they are not within 20%, keep both numbers
      deduplicated_values <- sim_vc_list
    }
  } else {
    # If there's only one number after removing NA
    deduplicated_values <- sim_vc_list
  }
  sim_vc_list <- deduplicated_values

  if (is.na(sim_vp_list) == T) {
    # Two-compartment model
    start.time <- Sys.time()
    # from linear to nonlinear
    # 10:1, 5:1, 2:1, 1:1, 1:2, 1:5, 1:10
    vc_vp_ratio_range <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)

    combs <- expand.grid(sim_vc_list, vc_vp_ratio_range)
    combs_df <- data.frame(combs)
    combs_df$Var3 <- signif(combs_df$Var1 / combs_df$Var2)

    colnames(combs_df) <- c("estvc", "vc_vp_ratio", "estvp")

    sim_q_list <- c(1, 10, estcl)

    base_values <-
      sim_q_list[-length(sim_q_list)]  # The other base values (1 and 10)

    # Check if `estcl` is within 20% of any base value
    remove_estcl <-
      any(abs(estcl - base_values) / base_values <= 0.2)

    # Remove `estcl` if it is within 20% of any base value
    if (remove_estcl) {
      sim_q_list <- base_values  # Remove the `estcl` value
    }

  }

  sim.2cmpt.results.all <- NULL

  for (loop2cmpt in 1:length(combs_df)) {
    input.vc <-  combs_df[loop2cmpt, 1]
    input.vp <-  combs_df[loop2cmpt, 3]

    for (estq in sim_q_list) {
      if (noniv_flag == 0) {
        sim.lst <- Fit_2cmpt_iv(
          data = dat[dat$EVID != 2, ],
          est.method = "rxSolve",
          input.cl = estcl,
          input.vc2cmpt = input.vc,
          input.vp2cmpt = input.vp,
          input.q2cmpt = estq,
          input.add = 0
        )
      }

      if (noniv_flag == 1) {
        sim.lst <- Fit_2cmpt_oral(
          data = dat[dat$EVID != 2, ],
          est.method = "rxSolve",
          input.ka = estka,
          input.cl = estcl,
          input.vc2cmpt = input.vc,
          input.vp2cmpt = input.vp,
          input.q2cmpt = estq,
          input.add = 0
        )
      }

      sim.APE <-
        round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[1], 1)
      sim.MAE <-
        round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[2], 1)
      sim.MAPE <-
        round(metrics.(pred.x = sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[3], 1)
      sim.RMSE <-
        round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[4], 1)
      sim.rRMSE <-
        round(metrics.(pred.x =  sim.lst$cp , obs.y = dat[dat$EVID == 0, ]$DV)[5], 1)

      end.time <- Sys.time()
      time.spent <- round(difftime(end.time, start.time), 4)

      sim.2cmpt.results <- data.frame(
        vc = input.vc,
        vp = input.vp,
        q = estq,
        sim.2cmpt.APE = sim.APE,
        sim.2cmpt.MAE = sim.MAE,
        sim.2cmpt.MAPE = sim.MAPE,
        sim.2cmpt.RMSE = sim.RMSE,
        sim.2cmpt.rRMSE = sim.rRMSE,

        time.spent = time.spent
      )

      sim.2cmpt.results.all <-
        rbind(sim.2cmpt.results.all, sim.2cmpt.results)

      rm(sim.lst) # To avoid occupying memory
      gc()
    }
  }
  return(sim.2cmpt.results.all)
}

#' Sensitivity analysis for a three-compartment model
#'
#' This function performs sensitivity analysis by testing a series of potential ratios of vc to vp and vp to v2.
#' @param dat A data frame containing the pharmacokinetic data
#' @param sim_vc_list A list of central compartment volumes (Vc) to simulate (optional).
#' @param sim_vp_list A list of first peripheral compartment volumes (Vp) to simulate (optional).
#' @param sim_vp2_list A list of second peripheral compartment volumes (Vp2) to simulate (optional).
#' @param estcl Estimated clearance.
#' @param estvc Estimated cenntral volume of distribution (optional, required if `sim_vc_list`, `sim_vp_list`and `sim_vp2_list`  are not provided).
#' @param estq Estimated intercompartmental clearance (optional, default is `estcl`).
#' @param estq2 Estimated intercompartmental clearance between central compartment and second peripheral compartment (optional, default is `estcl`).
#' @return A data frame containing the Vc, Vp, Vp2 APE, MAPE, and time spent for each simulation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' sim_results<-sim_sens_3cmpt(dat, estcl = 4, sim_vc_list = 70)
#' sim_results
#' @export
#'

sim_sens_3cmpt <- function(dat,
                           sim_vc_list=NA,
                           sim_vp_list=NA,
                           sim_vp2_list=NA,
                           sim_q_list=NA,
                           estcl,
                           estka=NA,
                           noniv_flag=0) {

# default value
  if (missing(estcl)) {
    stop("Error: no estimated clearance for simulation was provided")
  }

  if (missing(sim_vc_list)) {
      stop("Error: No planned central volume of distribution (Vd) was specified for the simulation.")
   }

  if (is.na(sim_vp_list)==T || is.na(sim_vp2_list)==T || is.na(sim_q_list)==F) {

    sim_vc_list <- na.omit(sim_vc_list)
    sim_vc_list <- sort(sim_vc_list)

    # Check if the two numbers are within 20% relative difference
    if (length(sim_vc_list) == 2) {
      relative_difference <-
        abs(sim_vc_list[2] - sim_vc_list[1]) / sim_vc_list[1]
      if (relative_difference <= 0.2) {
        # If they are within 20%, keep the smaller number
        deduplicated_values <- sim_vc_list[1]
      } else {
        # If they are not within 20%, keep both numbers
        deduplicated_values <- sim_vc_list
      }
    } else {
      # If there's only one number after removing NA
      deduplicated_values <- sim_vc_list
    }
    sim_vc_list <- deduplicated_values


      # generate the default parameter list
      start.time <- Sys.time()
      vc_vp_ratio_range <- c(5, 2, 1, 0.5, 0.2)
      vp_vp2_ratio_range <- c(5, 2, 1, 0.5, 0.2)

      # Generate the matrix of all combinations
      combs <- expand.grid(sim_vc_list, vc_vp_ratio_range,  vp_vp2_ratio_range)
      combs_df <- data.frame(combs)

      colnames(combs_df) <- c("estvc", "vc_vp_ratio", "vp_vp2_ratio")

      combs_df$vc_vp2_ratio <-  combs_df[, 2] * combs_df[, 3]
      # combs_df$vp_ratio_range <-  combs_df[, 3] *  combs_df[, 4]
      # combs_df$vp2_ratio_range <-  combs_df[, 4]

      combs_df$estvp<-combs_df$estvc/combs_df$vc_vp_ratio
      combs_df$estvp2<-combs_df$estvc/combs_df$vc_vp2_ratio

      sim_q_list <- c(1, 10, estcl)

      base_values <-
        sim_q_list[-length(sim_q_list)]  # The other base values (1 and 10)

      # Check if `estcl` is within 20% of any base value
      remove_estcl <-
        any(abs(estcl - base_values) / base_values <= 0.2)

      # Remove `estcl` if it is within 20% of any base value
      if (remove_estcl) {
        sim_q_list <- base_values  # Remove the `estcl` value
      }


  }

  sim.3cmpt.results.all <- NULL

  for (loop3cmpt in 1:length(combs_df)) {
    input.vc <- combs_df[loop3cmpt,1]
    input.vp <- combs_df[loop3cmpt,5]
    input.vp2 <- combs_df[loop3cmpt,6]

   for (estq in sim_q_list){

    if (noniv_flag==0){
    sim.lst <- Fit_3cmpt_iv(
      data = dat[dat$EVID != 2,],
      est.method = "rxSolve",
      input.cl = estcl,
      input.vc3cmpt = input.vc,
      input.vp3cmpt = input.vp,
      input.vp23cmpt = input.vp2,
      input.q3cmpt = estq,
      input.q23cmpt = estq,
      input.add = 0
    )
    }
    if (noniv_flag==1){
    sim.lst <- Fit_3cmpt_oral(
      data = dat[dat$EVID != 2,],
      est.method = "rxSolve",
      input.ka = estka,
      input.cl = estcl,
      input.vc3cmpt = input.vc,
      input.vp3cmpt = input.vp,
      input.vp23cmpt = input.vp2,
      input.q3cmpt = estq,
      input.q23cmpt = estq,
      input.add = 0
    )
    }

     sim.APE <-  round(metrics.(pred.x =  sim.lst$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],1)
     sim.MAE <-  round(metrics.(pred.x =  sim.lst$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],1)
     sim.MAPE <- round( metrics.(pred.x = sim.lst$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],1)
     sim.RMSE <-  round(metrics.(pred.x =  sim.lst$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],1)
     sim.rRMSE <- round( metrics.(pred.x =  sim.lst$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],1)

    end.time <- Sys.time()
    time.spent <- round(difftime(end.time, start.time), 4)

    sim.3cmpt.results <- data.frame(
      vc = input.vc,
      vp = input.vp,
      vp2 = input.vp2,
      q =  estq,
      q2 = estq,

      sim.3cmpt.APE = sim.APE,
      sim.3cmpt.MAE = sim.MAE,
      sim.3cmpt.MAPE = sim.MAPE,
      sim.3cmpt.RMSE = sim.RMSE,
      sim.3cmpt.rRMSE = sim.rRMSE,

      time.spent = time.spent
    )
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results)

    rm(sim.lst) # To avoid occupying memory
    gc()

   }
  }

  return(sim.3cmpt.results.all)

}
