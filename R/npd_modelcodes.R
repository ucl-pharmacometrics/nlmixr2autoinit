#' Fit intravenous pharmacokinetic data to a one-compartment linear elimination model
#'
#' Fits intravenous (IV) pharmacokinetic data to a one-compartment model with
#' first-order elimination using the naive pooled data approach. Supports
#' multiple estimation methods provided by nlmixr2 and can optionally return
#' only predicted concentrations to support efficient simulation workflows.
#'
#' @param data A data frame of IV pharmacokinetic data formatted for nlmixr2.
#' @param est.method Estimation method to use in nlmixr2. Must be one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vd Initial estimate of volume of distribution (V).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @importFrom nlmixr2 ini model nlmixr2
#' @importFrom rxode2 rxControl
#' @importFrom nlmixr2est nlsControl nlmControl nlminbControl foceiControl
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_1CPT
#'
#' # Fit using 'nls' with default control
#' Fit_1cmpt_iv(
#'   data = dat,
#'   est.method = "nls",
#'   input.cl = 4,
#'   input.vd = 70,
#'   input.add = 1
#' )
#'
#' # Fit using 'focei' with custom control settings
#' Fit_1cmpt_iv(
#'   data = dat,
#'   est.method = "focei",
#'   input.cl = 4,
#'   input.vd = 70,
#'   input.add = 1,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200, atol = 1e-8))
#' )
#'
#' # Return only predicted concentrations
#' Fit_1cmpt_iv(
#'  data = dat,
#'  est.method = "rxSolve",
#'  input.cl = 4,
#'  input.vd = 70,
#'  input.add = 1,
#'  return.pred.only = TRUE
#')
#' }
#' @export
#'

Fit_1cmpt_iv <- function(data,
                         est.method,
                         input.cl,
                         input.vd,
                         input.add,
                         return.pred.only = FALSE,
                         ...) {
  # Assign parameters to global environment (<<-)
  input.cl  <<- input.cl
  input.vd  <<- input.vd
  input.add <<- input.add

  # Define the IV one-compartment model with linear clearance
  iv <- function() {
    ini({
      tcl     <- round(log(input.cl), 2)  # log(CL)
      tv      <- round(log(input.vd), 2)  # log(V)
      add.err <- input.add
    })
    model({
      cl <- exp(tcl)
      v  <- exp(tv)
      k  <- cl / v
      d / dt(centre) = -k * centre
      cp = centre / v
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)

    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.1cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = iv,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.1cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = iv,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.1cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.1cmpt.lst)
}


#' Fit intravenous pharmacokinetic data to a one-compartment model with Michaelis-Menten elimination
#'
#' Fits intravenous (IV) pharmacokinetic data to a one-compartment model with
#' Michaelis-Menten (nonlinear) elimination using the naive pooled data approach.
#' Supports multiple estimation methods available in nlmixr2, and optionally
#' returns only predicted concentrations to reduce memory use in simulation workflows.
#'
#' @param data A data frame of IV pharmacokinetic data formatted for nlmixr2.
#' @param est.method Estimation method to use in nlmixr2, one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.vmax Initial estimate of the maximum elimination rate (Vmax).
#' @param input.km Initial estimate of the Michaelis constant (Km).
#' @param input.vd Initial estimate of the volume of distribution (V).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Optional arguments passed to \code{nlmixr2()}, such as a custom
#'   \code{control = foceiControl(...)} or other control objects.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_1CPTMM
#' # Fit using 'nls' method with default control
#' Fit_1cmpt_mm_iv(
#'   data = dat,
#'   est.method = "nls",
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10
#' )
#' # Fit using 'focei' with custom rxControl settings
#' Fit_1cmpt_mm_iv(
#'   data = dat,
#'   est.method = "focei",
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 300, atol = 1e-8))
#' )
#' # Return only predicted concentrations
#' Fit_1cmpt_mm_iv(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export

Fit_1cmpt_mm_iv <- function(data,
                            est.method,
                            input.vmax,
                            input.km,
                            input.vd,
                            input.add,
                            return.pred.only = FALSE,
                            ...) {
  # Assign parameters to global environment (<<-)
  input.vmax <<- input.vmax
  input.km   <<- input.km
  input.vd   <<- input.vd
  input.add  <<- input.add

  # Define the IV MM model
  iv.mm <- function() {
    ini({
      lvmax <- round(log(input.vmax), 2)
      lkm   <- round(log(input.km), 2)
      tv    <- round(log(input.vd), 2)
      add.err <- input.add
    })
    model({
      vmax <- exp(lvmax)
      km   <- exp(lkm)
      v    <- exp(tv)
      d / dt(centre) <- -(vmax / (km + centre / v)) / v * centre
      cp <- centre / v
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)

    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = iv.mm,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = iv.mm,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <-
      data.frame(cp = fit.1cmpt.mm.lst$cp)
    return(sim_out)
  }

  return(fit.1cmpt.mm.lst)
}


#' Fit intravenous pharmacokinetic data to a two-compartment linear elimination model
#'
#' Fits intravenous (IV) pharmacokinetic data to a two-compartment model with
#' first-order elimination using the naive pooled data approach. Supports
#' multiple estimation methods provided by nlmixr2 and can optionally return
#' only predicted concentrations to support efficient simulation workflows.
#'
#' @param data A data frame containing IV pharmacokinetic data formatted for nlmixr2,
#' @param est.method Estimation method to use in nlmixr2. Must be one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vc2cmpt Initial estimate of central volume of distribution (V1).
#' @param input.vp2cmpt Initial estimate of peripheral volume of distribution (V2).
#' @param input.q2cmpt Initial estimate of inter-compartmental clearance (Q).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_2CPT
#' # Fit using 'nls' method with default control settings
#' Fit_2cmpt_iv(
#'   data = dat,
#'   est.method = "nls",
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 4,
#'   input.add = 10
#' )

#' # Fit using 'focei' with custom control settings
#' Fit_2cmpt_iv(
#'   data = dat,
#'   est.method = "focei",
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 4,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200))
#' )

#' # Return only predicted concentrations (cp) for all timepoints
#' Fit_2cmpt_iv(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 4,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export
Fit_2cmpt_iv <- function(data,
                         est.method,
                         input.cl,
                         input.vc2cmpt,
                         input.vp2cmpt,
                         input.q2cmpt,
                         input.add,
                         return.pred.only = FALSE,
                         ...) {
  # Assign parameters to global environment (<<-)
  input.cl       <<- input.cl
  input.vc2cmpt  <<- input.vc2cmpt
  input.vp2cmpt  <<- input.vp2cmpt
  input.q2cmpt   <<- input.q2cmpt
  input.add      <<- input.add

  # Define the IV two-compartment model
  iv2 <- function() {
    ini({
      tcl     <- round(log(input.cl), 2)       # log(CL)
      tv1     <- round(log(input.vc2cmpt), 2)   # log(V1)
      tv2     <- round(log(input.vp2cmpt), 2)   # log(V2)
      tq      <- round(log(input.q2cmpt), 2)    # log(Q)
      add.err <- input.add
    })
    model({
      cl <- exp(tcl)
      v1 <- exp(tv1)
      v2 <- exp(tv2)
      q  <- exp(tq)

      k   <- cl / v1
      k12 <- q / v1
      k21 <- q / v2

      d / dt(centre) = -k * centre - k12 * centre + k21 * A2
      d / dt(A2)     = k12 * centre - k21 * A2

      cp = centre / v1
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)
    ctrl <- switch(
      est.method,
      "rxSolve"  =  rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.2cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = iv2,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.2cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = iv2,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.2cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.2cmpt.lst)
}



#' Fit intravenous pharmacokinetic data to a three-compartment linear elimination model
#'
#' Fits intravenous (IV) pharmacokinetic data to a three-compartment model with
#' linear (first-order) elimination using the naive pooled data approach. Supports
#' multiple estimation methods provided by nlmixr2 and can optionally return
#' only predicted concentrations to support efficient simulation workflows.
#'
#' @param data A data frame containing IV pharmacokinetic data formatted for nlmixr2.
#' @param est.method Estimation method to use in nlmixr2. Must be one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vc3cmpt Initial estimate of central volume of distribution (V1).
#' @param input.vp3cmpt Initial estimate of first peripheral volume of distribution (V2).
#' @param input.vp23cmpt Initial estimate of second peripheral volume of distribution (V3).
#' @param input.q3cmpt Initial estimate of first inter-compartmental clearance (Q1).
#' @param input.q23cmpt Initial estimate of second inter-compartmental clearance (Q2).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_2CPT
#'
#' # Fit using 'nls' method
#' Fit_3cmpt_iv(
#'   data = dat,
#'   est.method = "nls",
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10
#' )

#' # Fit using 'focei' with custom control
#' Fit_3cmpt_iv(
#'   data = dat,
#'   est.method = "focei",
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200))
#' )

#' # Return only predicted concentrations
#' Fit_3cmpt_iv(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export
Fit_3cmpt_iv <- function(data,
                         est.method,
                         input.cl,
                         input.vc3cmpt,
                         input.vp3cmpt,
                         input.vp23cmpt,
                         input.q3cmpt,
                         input.q23cmpt,
                         input.add,
                         return.pred.only = FALSE,
                         ...) {
  input.cl       <<- input.cl
  input.vc3cmpt  <<- input.vc3cmpt
  input.vp3cmpt  <<- input.vp3cmpt
  input.vp23cmpt <<- input.vp23cmpt
  input.q3cmpt   <<- input.q3cmpt
  input.q23cmpt  <<- input.q23cmpt
  input.add      <<- input.add

  iv3 <- function() {
    ini({
      tcl    <- round(log(input.cl), 2)
      tv1    <- round(log(input.vc3cmpt), 2)
      tv2    <- round(log(input.vp3cmpt), 2)
      tv3    <- round(log(input.vp23cmpt), 2)
      tq     <- round(log(input.q3cmpt), 2)
      tq2    <- round(log(input.q23cmpt), 2)
      add.err <- input.add
    })
    model({
      cl <- exp(tcl)
      v1 <- exp(tv1)
      v2 <- exp(tv2)
      v3 <- exp(tv3)
      q  <- exp(tq)
      q2 <- exp(tq2)

      k    <- cl / v1
      k12  <- q / v1
      k21  <- q / v2
      k13  <- q2 / v1
      k31  <- q2 / v3

      d / dt(centre) = -k * centre - k12 * centre + k21 * A2 - k13 * centre + k31 * A3
      d / dt(A2)     = k12 * centre - k21 * A2
      d / dt(A3)     = k13 * centre - k31 * A3

      cp = centre / v1
      cp ~ add(add.err)
    })
  }

  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)
    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.3cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = iv3,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.3cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = iv3,
      data = data,
      est = est.method,
      ...
    )))
  }

  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.3cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.3cmpt.lst)
}

#' Fit oral pharmacokinetic data to a one-compartment linear elimination model
#'
#' Fits oral pharmacokinetic data to a one-compartment model with
#' first-order absorption and first-order elimination using the naive pooled data approach.
#' Supports multiple estimation methods provided by nlmixr2 and can optionally return
#' only predicted concentrations to support efficient simulation workflows.
#'
#' @param data A data frame containing oral pharmacokinetic data formatted for nlmixr2,
#' @param est.method Estimation method to use in nlmixr2. Must be one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.ka Initial estimate of the absorption rate constant (ka).
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vd Initial estimate of volume of distribution (V).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Oral_1CPT
#'
#' # Fit using 'nls'
#' Fit_1cmpt_oral(
#'   data = dat,
#'   est.method = "nls",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vd = 70,
#'   input.add = 10
#' )
#'
#' # Fit using 'focei' with custom control
#' Fit_1cmpt_oral(
#'   data = dat,
#'   est.method = "focei",
#'   input.ka = 2,
#'   input.cl = 4,
#'   input.vd = 70,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200))
#' )
#'
#' # Return only predicted concentrations
#' Fit_1cmpt_oral(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vd = 70,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export

Fit_1cmpt_oral <- function(data,
                           est.method,
                           input.ka,
                           input.cl,
                           input.vd,
                           input.add,
                           return.pred.only = FALSE,
                           ...) {
  # Assign parameters to global environment (<<-)
  input.ka  <<- input.ka
  input.cl  <<- input.cl
  input.vd  <<- input.vd
  input.add <<- input.add

  # Define the oral one-compartment model with first-order absorption and elimination
  oral <- function() {
    ini({
      tka     <- round(log(input.ka), 2)  # log(ka)
      tcl     <- round(log(input.cl), 2)  # log(CL)
      tv      <- round(log(input.vd), 2)  # log(V)
      add.err <- input.add
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v  <- exp(tv)
      k  <- cl / v

      d / dt(depot)  = -ka * depot
      d / dt(centre) =  ka * depot - k * centre

      cp = centre / v
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)

    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.1cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = oral,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.1cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = oral,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.1cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.1cmpt.lst)
}


#' Fit oral pharmacokinetic data to a one-compartment model with Michaelis-Menten elimination
#'
#' Fits oral pharmacokinetic data to a one-compartment model with
#' first-order absorption and Michaelis-Menten (nonlinear) elimination using the naive pooled data approach.
#' Supports multiple estimation methods available in nlmixr2, and optionally
#' returns only predicted concentrations to reduce memory use in simulation workflows.
#'
#' @param data A data frame of oral pharmacokinetic data formatted for nlmixr2.
#' @param est.method Estimation method to use in nlmixr2, one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.ka Initial estimate of the absorption rate constant (ka).
#' @param input.vmax Initial estimate of the maximum elimination rate (Vmax).
#' @param input.km Initial estimate of the Michaelis constant (Km).
#' @param input.vd Initial estimate of the volume of distribution (V).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Optional arguments passed to \code{nlmixr2()}, such as a custom
#'   \code{control = foceiControl(...)} or other control objects.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with columns \code{cp} (predicted concentration) and \code{DV} (observed data).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Oral_1CPTMM
#' # Fit using 'nls'
#' Fit_1cmpt_mm_oral(
#'   data = dat,
#'   est.method = "nls",
#'   input.ka = 1,
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10
#' )
#' # Fit using 'focei' with custom control
#' Fit_1cmpt_mm_oral(
#'   data = dat,
#'   est.method = "focei",
#'   input.ka = 1,
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200))
#' )
#' # Return only predicted concentrations
#' Fit_1cmpt_mm_oral(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.ka = 1,
#'   input.vmax = 1000,
#'   input.km = 250,
#'   input.vd = 70,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export

Fit_1cmpt_mm_oral <- function(data,
                              est.method,
                              input.ka,
                              input.vmax,
                              input.km,
                              input.vd,
                              input.add,
                              return.pred.only = FALSE,
                              ...) {
  # Assign parameters to global environment (<<-)
  input.ka    <<- input.ka
  input.vmax  <<- input.vmax
  input.km    <<- input.km
  input.vd    <<- input.vd
  input.add   <<- input.add

  # Define the oral one-compartment MM model
  oral.mm <- function() {
    ini({
      tka    <- round(log(input.ka), 2)
      lvmax  <- round(log(input.vmax), 2)
      lkm    <- round(log(input.km), 2)
      tv     <- round(log(input.vd), 2)
      add.err <- input.add
    })
    model({
      ka   <- exp(tka)
      vmax <- exp(lvmax)
      km   <- exp(lkm)
      v    <- exp(tv)

      d / dt(depot)  = -ka * depot
      d / dt(centre) = ka * depot - (vmax / (km + centre / v)) / v * centre

      cp <- centre / v
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)

    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = oral.mm,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = oral.mm,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.1cmpt.mm.lst$cp)
    return(sim_out)
  }

  return(fit.1cmpt.mm.lst)
}


#' Fit oral pharmacokinetic data to a two-compartment model
#'
#' Fits oral pharmacokinetic data to a two-compartment model with
#' first-order absorption and first-order elimination using the naive pooled data approach.
#' Supports multiple estimation methods available in nlmixr2, and optionally
#' returns only predicted concentrations to support simulation workflows.
#'
#' @param data A data frame containing oral pharmacokinetic data formatted for nlmixr2,
#' @param est.method Estimation method to use in nlmixr2, one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.ka Initial estimate of the absorption rate constant (ka).
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vc2cmpt Initial estimate of central volume of distribution (V1).
#' @param input.vp2cmpt Initial estimate of peripheral volume of distribution (V2).
#' @param input.q2cmpt Initial estimate of inter-compartmental clearance (Q).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Oral_2CPT

# Fit using 'nls'
#' Fit_2cmpt_oral(
#'   data = dat,
#'   est.method = "nls",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 10,
#'   input.add = 10
#' )

# Fit using 'focei' with custom rxControl settings
#' Fit_2cmpt_oral(
#'   data = dat,
#'   est.method = "focei",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 4,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200, atol = 1e-8))
#' )
# Return only predicted concentrations
#' Fit_2cmpt_oral(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc2cmpt = 70,
#'   input.vp2cmpt = 40,
#'   input.q2cmpt = 4,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export

Fit_2cmpt_oral <- function(data,
                           est.method,
                           input.ka,
                           input.cl,
                           input.vc2cmpt,
                           input.vp2cmpt,
                           input.q2cmpt,
                           input.add,
                           return.pred.only = FALSE,
                           ...) {
  # Assign parameters to global environment (<<-)
  input.ka       <<- input.ka
  input.cl       <<- input.cl
  input.vc2cmpt  <<- input.vc2cmpt
  input.vp2cmpt  <<- input.vp2cmpt
  input.q2cmpt   <<- input.q2cmpt
  input.add      <<- input.add

  # Define the oral two-compartment model
  oral2 <- function() {
    ini({
      tka     <- round(log(input.ka), 2)        # log(ka)
      tcl     <- round(log(input.cl), 2)        # log(CL)
      tv1     <- round(log(input.vc2cmpt), 2)   # log(V1)
      tv2     <- round(log(input.vp2cmpt), 2)   # log(V2)
      tq      <- round(log(input.q2cmpt), 2)    # log(Q)
      add.err <- input.add
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v1 <- exp(tv1)
      v2 <- exp(tv2)
      q  <- exp(tq)

      k   <- cl / v1
      k12 <- q / v1
      k21 <- q / v2

      d / dt(depot) = -ka * depot
      d / dt(A1)    = ka * depot - k * A1 - k12 * A1 + k21 * A2
      d / dt(A2)    = k12 * A1 - k21 * A2

      cp = A1 / v1
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)
    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.2cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = oral2,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.2cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = oral2,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.2cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.2cmpt.lst)
}


#' Fit oral pharmacokinetic data to a three-compartment linear elimination model
#'
#' Fits oral pharmacokinetic data to a three-compartment model with
#' first-order absorption and first-order elimination using the naive pooled data approach.
#' Supports multiple estimation methods provided by nlmixr2 and can optionally return
#' only predicted concentrations to support efficient simulation workflows.
#'
#' @param data A data frame containing oral pharmacokinetic data formatted for nlmixr2.
#' @param est.method Estimation method to use in nlmixr2. Must be one of:
#'   \code{"rxSolve"}, \code{"nls"}, \code{"nlm"}, \code{"nlminb"}, or \code{"focei"}.
#' @param input.ka Initial estimate of the absorption rate constant (ka).
#' @param input.cl Initial estimate of clearance (CL).
#' @param input.vc3cmpt Initial estimate of central volume of distribution (V1).
#' @param input.vp3cmpt Initial estimate of first peripheral volume of distribution (V2).
#' @param input.vp23cmpt Initial estimate of second peripheral volume of distribution (V3).
#' @param input.q3cmpt Initial estimate of first inter-compartmental clearance (Q1).
#' @param input.q23cmpt Initial estimate of second inter-compartmental clearance (Q2).
#' @param input.add Initial estimate of the additive residual error.
#' @param return.pred.only Logical; if \code{TRUE}, returns a data frame with
#'   only predicted concentrations (\code{cp}) for all observations in the input data.
#' @param ... Additional arguments passed to \code{nlmixr2()}, such as a user-defined
#'   \code{control = foceiControl(...)} or other control settings.
#'
#' @return If \code{return.pred.only = TRUE}, returns a \code{data.frame}
#'   with a single column \code{cp} (predicted concentrations).
#'   Otherwise, returns a fitted model object produced by nlmixr2.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Oral_2CPT
#'
#' # Fit using 'nls'
#' Fit_3cmpt_oral(
#'   data = dat,
#'   est.method = "nls",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10
#' )
#'
#' # Fit using 'focei' with custom control
#' Fit_3cmpt_oral(
#'   data = dat,
#'   est.method = "focei",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10,
#'   control = foceiControl(rxControl = rxControl(maxSS = 200))
#' )
#'
#' # Return only predicted concentrations
#' Fit_3cmpt_oral(
#'   data = dat,
#'   est.method = "rxSolve",
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc3cmpt = 70,
#'   input.vp3cmpt = 35,
#'   input.vp23cmpt = 35,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4,
#'   input.add = 10,
#'   return.pred.only = TRUE
#' )
#' }
#' @export

Fit_3cmpt_oral <- function(data,
                           est.method,
                           input.ka,
                           input.cl,
                           input.vc3cmpt,
                           input.vp3cmpt,
                           input.vp23cmpt,
                           input.q3cmpt,
                           input.q23cmpt,
                           input.add,
                           return.pred.only = FALSE,
                           ...) {
  # Assign parameters to global environment
  input.ka       <<- input.ka
  input.cl       <<- input.cl
  input.vc3cmpt  <<- input.vc3cmpt
  input.vp3cmpt  <<- input.vp3cmpt
  input.vp23cmpt <<- input.vp23cmpt
  input.q3cmpt   <<- input.q3cmpt
  input.q23cmpt  <<- input.q23cmpt
  input.add      <<- input.add

  # Define the oral three-compartment model
  oral3 <- function() {
    ini({
      tka    <- round(log(input.ka), 2)
      tcl    <- round(log(input.cl), 2)
      tv1    <- round(log(input.vc3cmpt), 2)
      tv2    <- round(log(input.vp3cmpt), 2)
      tv3    <- round(log(input.vp23cmpt), 2)
      tq     <- round(log(input.q3cmpt), 2)
      tq2    <- round(log(input.q23cmpt), 2)
      add.err <- input.add
    })
    model({
      ka  <- exp(tka)
      cl  <- exp(tcl)
      v1  <- exp(tv1)
      v2  <- exp(tv2)
      v3  <- exp(tv3)
      q1  <- exp(tq)
      q2  <- exp(tq2)

      k    <- cl / v1
      k12  <- q1 / v1
      k21  <- q1 / v2
      k13  <- q2 / v1
      k31  <- q2 / v3

      d / dt(depot) = -ka * depot
      d / dt(A1)    = ka * depot - k * A1 - k12 * A1 + k21 * A2 - k13 * A1 + k31 * A3
      d / dt(A2)    = k12 * A1 - k21 * A2
      d / dt(A3)    = k13 * A1 - k31 * A3

      cp = A1 / v1
      cp ~ add(add.err)
    })
  }

  # Parse additional arguments
  dot_args <- list(...)
  user_control <- "control" %in% names(dot_args)

  if (!user_control) {
    maxSSv <- 100
    rx_ctrl <- rxControl(maxSS = maxSSv)
    ctrl <- switch(
      est.method,
      "rxSolve"  = rx_ctrl,
      "nls"      = nlsControl(rxControl = rx_ctrl),
      "nlm"      = nlmControl(rxControl = rx_ctrl),
      "nlminb"   = nlminbControl(rxControl = rx_ctrl),
      "focei"    = foceiControl(rxControl = rx_ctrl),
      stop("Unsupported estimation method: ", est.method)
    )

    fit.3cmpt.lst <- suppressMessages(suppressWarnings(
      nlmixr2(
        object = oral3,
        data = data,
        est = est.method,
        control = ctrl
      )
    ))
  } else {
    fit.3cmpt.lst <- suppressMessages(suppressWarnings(nlmixr2(
      object = oral3,
      data = data,
      est = est.method,
      ...
    )))
  }

  # Return only predicted concentrations if requested
  if (return.pred.only) {
    sim_out <- data.frame(cp = fit.3cmpt.lst$cp)
    return(sim_out)
  }

  return(fit.3cmpt.lst)
}
