#' Fit intravenous pharmacokinetic data to a one-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", "focei"). The default value is "nls".
#' @param input.cl Initial estimate for clearance.
#' @param input.vd Initial estimate for volume of distribution.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPT
#' Fit_1cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vd=70,input.add=1)
#' @export

Fit_1cmpt_iv <- function(data,
                    est.method,
                    input.cl,
                    input.vd,
                    input.add) {

  input.cl <<- input.cl
  input.vd <<- input.vd
  input.add <<- input.add

  iv <- function() {
    ini({
      tcl <- round(log(input.cl), 2) # Clearance
      tv  <- round(log(input.vd), 2) # Volume of distribution
      add.err <-  input.add # Additive error
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
      # Temporary setting
      maxSSv=100

      if (est.method=="rxSolve"){
        fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
      }

      if (est.method=="nls"){
       fit.1cmpt.lst <- suppressMessages(suppressWarnings( nlmixr2( object = iv, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
      }

      if (est.method=="nlm"){
        fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
      }

       if (est.method=="nlminb"){
        fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
       }

      if (est.method=="focei"){
        fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
      }

  return(fit.1cmpt.lst)
}


#' Fit intravenous pharmacokinetic data to a one-compartment model with Michaelis-Menten Kinetics
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model with nonlinear elimination kinetics for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.vmax Initial estimate for Vmax.
#' @param input.km Initial estimate for Km.
#' @param input.vd Initial estimate for volume of distribution.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPTMM
#' Fit_1cmpt_mm_iv(dat=dat,est.method="nls",input.vmax=1000,input.km=100,input.vd=70,input.add=1)
#' @export

Fit_1cmpt_mm_iv <- function(data,
                       est.method,
                       input.vmax,
                       input.km,
                       input.vd,
                       input.add) {

  input.vmax <<- input.vmax
  input.km <<- input.km
  input.vd <<- input.vd
  input.add <<- input.add

  iv.mm <- function() {
    ini({
      lvmax  <- round(log(input.vmax), 2) # Maximum elimination rate
      lkm  <- round(log(input.km), 2) # Michaelis constant
      tv  <- round(log(input.vd), 2) # Volume of distribution
      add.err <- input.add # additive error
    })
    model({
      vmax = exp(lvmax)
      km = exp(lkm)
      v  <- exp(tv)
      d / dt(centre) <- -(vmax / (km + centre / v)) / v * centre
      cp <- centre / v
      cp ~ add(add.err)
    })
  }

  maxSSv=100

  if (est.method=="rxSolve"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv.mm, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings( nlmixr2( object = iv.mm, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv.mm, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv.mm, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = iv.mm, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  return(fit.1cmpt.mm.lst)
}


#' Fit intravenous pharmacokinetic data to a two-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl Initial estimate for clearance.
#' @param input.vc2cmpt Initial estimate for the central volume of distribution in two-compartment model.
#' @param input.vp2cmpt Initial estimate for the peripheral volume of distribution in a two-compartment model.
#' @param input.q2cmpt Initial estimate for inter-compartmental clearance in a two-compartment model.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' Fit_2cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vc2cmpt=70,input.vp2cmpt=35,input.q2cmpt=10,input.add=1)
#' @export

Fit_2cmpt_iv<-function(data,
                        est.method,
                        input.cl,
                        input.vc2cmpt,
                        input.vp2cmpt,
                        input.q2cmpt,
                        input.add){
  input.cl<<- input.cl
  input.vc2cmpt <<-input.vc2cmpt
  input.vp2cmpt<<- input.vp2cmpt
  input.q2cmpt<<-  input.q2cmpt
  input.add <<- input.add


  iv2 <- function() {
    ini({
      tcl <- round(log(input.cl),2) # Clearance
      tv1  <- round(log(input.vc2cmpt),2) # Central Volume of distribution
      tv2  <- round(log(input.vp2cmpt),2) # Peripheral volume of distribution
      tq <- round(log(input.q2cmpt),2) # Inter-compartmental clearance
      add.err <- input.add
    })
    model({
      cl <- exp(tcl)
      v1  <- exp(tv1)
      v2  <- exp(tv2)
      q  <- exp(tq)
      k  <- cl / v1

      k12<-q/v1
      k21<-q/v2

      d/dt(centre) = - k * centre   - k12 * centre + k21*A2
      d/dt(A2) =   k12 * centre -  k21 * A2
      cp = centre / v1
      cp ~ add(add.err)
    })
  }

  maxSSv=100

  if (est.method=="rxSolve"){
     fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv2, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
     fit.2cmpt.lst  <- suppressMessages(suppressWarnings( nlmixr2( object = iv2, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
     fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv2, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
     fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv2, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
     fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv2, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }
  return(fit.2cmpt.lst)
}


#' Fit intravenous pharmacokinetic data to a three-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl Initial estimate for clearance.
#' @param input.vc3cmpt Initial estimate for the central volume of distribution in a three-compartment model.
#' @param input.vp3cmpt Initial estimate for the first peripheral volume of distribution in a three-compartment model.
#' @param input.vp23cmpt Initial estimate for the second peripheral volume of distribution in a three-compartment model.
#' @param input.q3cmpt Initial estimate for first inter-compartmental clearance in a three-compartment model.
#' @param input.q23cmpt Initial estimate for second inter-compartmental clearance in a three-compartment model.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' Fit_3cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vc3cmpt=70,input.vp3cmpt=35,input.vp23cmpt=35,input.q3cmpt=10,input.q23cmpt=10,input.add=1)
#' @export

Fit_3cmpt_iv<-function(data,
                        est.method,
                        input.cl,
                        input.vc3cmpt,
                        input.vp3cmpt,
                        input.vp23cmpt,
                        input.q3cmpt,
                        input.q23cmpt,
                        input.add){

  input.cl<<-input.cl
  input.vc3cmpt <<- input.vc3cmpt
  input.vp3cmpt<<- input.vp3cmpt
  input.vp23cmpt<<- input.vp23cmpt
  input.q3cmpt<<-  input.q3cmpt
  input.q23cmpt<<-  input.q23cmpt
  input.add <<- input.add

  iv3 <- function() {
    ini({
      tcl <- round(log(input.cl),2) # Clearance
      tv1  <- round(log(input.vc3cmpt),2) # Central Volume of Distribution
      tv2  <- round(log(input.vp3cmpt),2) # First peripheral volume of distribution
      tv3  <- round(log(input.vp23cmpt),2) # Second peripheral volume of distribution
      tq <- round(log(input.q3cmpt),2) # First inter-compartmental clearance
      tq2 <- round(log(input.q23cmpt),2) # Second inter-compartmental clearance
      add.err <- input.add
    })
    model({
      cl <- exp(tcl)
      v1  <- exp(tv1)
      v2  <- exp(tv2)
      v3  <- exp(tv3)
      q  <- exp(tq)
      q2  <- exp(tq2)
      k  <- cl / v1
      k12<-q/v1
      k21<-q/v2
      k13<-q2/v1
      k31<-q2/v3
      d/dt(centre) = - k * centre   - k12 * centre + k21*A2
      d/dt(A2) =   k12 * centre -  k21 * A2
      d/dt(A3) =   k13 * centre -  k31 * A3
      cp = centre / v1
      cp ~ add(add.err)
    })
  }


  maxSSv=100

  if (est.method=="rxSolve"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv3, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.3cmpt.lst  <- suppressMessages(suppressWarnings( nlmixr2( object = iv3, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv3, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv3, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = iv3, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  return(fit.3cmpt.lst)
}


#' Fit oral pharmacokinetic data to a one-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", "focei"). The default value is "nls".
#' @param input.ka Initial estimate for absorption rate constant.
#' @param input.cl Initial estimate for clearance.
#' @param input.vd Initial estimate for volume of distribution.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Oral_1CPT
#' Fit_1cmpt_oral(dat=dat,est.method="nls",input.ka=1,input.cl=4,input.vd=70,input.add=1)
#' @export

Fit_1cmpt_oral <- function(data,
                           est.method,
                           input.ka,
                           input.cl,
                           input.vd,
                           input.add) {

  input.ka <<- input.ka
  input.cl <<- input.cl
  input.vd <<- input.vd
  input.add <<- input.add

  oral <- function() {
    ini({
      tka<-round(log(input.ka), 2) # Absorption rate
      tcl <- round(log(input.cl), 2) # Clearance
      tv  <- round(log(input.vd), 2) # Volume of distribution
      add.err <-  input.add # Additive error
    })
    model({

      ka<- exp(tka)
      cl <- exp(tcl)
      v  <- exp(tv)
      k  <- cl / v

      d / dt(depot) = - ka*depot
      d / dt(centre) = -k * centre + ka*depot

      cp = centre / v
      cp ~ add(add.err)
    })
  }

  # Temporary setting
  maxSSv=100

  if (est.method=="rxSolve"){
    fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.1cmpt.lst <- suppressMessages(suppressWarnings( nlmixr2( object = oral, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.1cmpt.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }
  return(fit.1cmpt.lst)
}


#' Fit oral pharmacokinetic data to a one-compartment model with Michaelis-Menten Kinetics
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model with nonlinear elimination kinetics for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka Initial estimate for absorption rate constant.
#' @param input.vmax Initial estimate for Vmax.
#' @param input.km Initial estimate for Km.
#' @param input.vd Initial estimate for volume of distribution.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Oral_1CPTMM
#' Fit_1cmpt_mm_oral(dat=dat,est.method="nls",input.ka=1,input.vmax=1000,input.km=100,input.vd=70,input.add=1)
#' @export

Fit_1cmpt_mm_oral <- function(data,
                              est.method,
                              input.ka,
                              input.vmax,
                              input.km,
                              input.vd,
                              input.add) {

  input.ka <<- input.ka
  input.vmax <<- input.vmax
  input.km <<- input.km
  input.vd <<- input.vd
  input.add <<- input.add

  oral.mm <- function() {
    ini({
      tka<-round(log(input.ka), 2) # Absorption rate
      lvmax  <- round(log(input.vmax), 2) # Maximum elimination rate
      lkm  <- round(log(input.km), 2) # Michaelis constant
      tv  <- round(log(input.vd), 2) # Volume of distribution
      add.err <- input.add # additive error
    })
    model({
      ka<- exp(tka)
      vmax = exp(lvmax)
      km = exp(lkm)
      v  <- exp(tv)
      d / dt(depot) = - ka*depot
      d / dt(centre) <- -(vmax / (km + centre / v)) / v * centre + ka*depot
      cp <- centre / v
      cp ~ add(add.err)
    })
  }


  maxSSv=100

  if (est.method=="rxSolve"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral.mm, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.1cmpt.mm.lst <- suppressMessages(suppressWarnings( nlmixr2( object = oral.mm, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral.mm, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral.mm, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.1cmpt.mm.lst <-  suppressMessages(suppressWarnings(nlmixr2( object = oral.mm, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  return(fit.1cmpt.mm.lst)
}


#' Fit oral pharmacokinetic data to a two-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka Initial estimate for absorption rate constant.
#' @param input.cl Initial estimate for clearance.
#' @param input.vc2cmpt Initial estimate for the central volume of distribution in two-compartment model.
#' @param input.vp2cmpt Initial estimate for the peripheral volume of distribution in a two-compartment model.
#' @param input.q2cmpt Initial estimate for inter-compartmental clearance in a two-compartment model.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Oral_2CPT
#' Fit_2cmpt_oral(dat=dat,est.method="nls",input.ka=1,input.cl=4,input.vc2cmpt=70,input.vp2cmpt=35,input.q2cmpt=10,input.add=1)
#' @export

Fit_2cmpt_oral<-function(data,
                         est.method,
                         input.ka,
                         input.cl,
                         input.vc2cmpt,
                         input.vp2cmpt,
                         input.q2cmpt,
                         input.add){

  input.ka<<-input.ka
  input.cl<<-input.cl
  input.vc2cmpt <<- input.vc2cmpt
  input.vp2cmpt<<- input.vp2cmpt
  input.q2cmpt<<-  input.q2cmpt
  input.add <<- input.add


  oral2 <- function() {
    ini({
      tka<-round(log(input.ka), 2) # Absorption rate
      tcl <- round(log(input.cl),2) # Clearance
      tv1  <- round(log(input.vc2cmpt),2) # Central Volume of distribution
      tv2  <- round(log(input.vp2cmpt),2) # Peripheral volume of distribution
      tq <- round(log(input.q2cmpt),2) # Inter-compartmental clearance
      add.err <- input.add
    })
    model({
      ka<- exp(tka)
      cl <- exp(tcl)
      v1  <- exp(tv1)
      v2  <- exp(tv2)
      q  <- exp(tq)
      k  <- cl / v1

      k12<-q/v1
      k21<-q/v2

      d / dt(depot) = - ka*depot
      d/dt(A1) = - k * A1   - k12 * A1 + k21*A2 + ka*depot
      d/dt(A2) =   k12 * A1 -  k21 * A2
      cp = A1 / v1
      cp ~ add(add.err)
    })
  }

  maxSSv=100

  if (est.method=="rxSolve"){
    fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral2, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.2cmpt.lst  <- suppressMessages(suppressWarnings( nlmixr2( object = oral2, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral2, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral2, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.2cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral2, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  return(fit.2cmpt.lst)
}


#' Fit oral pharmacokinetic data to a three-compartment model
#'
#' Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
#' @param data Intravenous pharmacokinetic data in the nlmixr2 format.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka Initial estimate for absorption rate constant.
#' @param input.cl Initial estimate for clearance.
#' @param input.vc3cmpt Initial estimate for the central volume of distribution in a three-compartment model.
#' @param input.vp3cmpt Initial estimate for the first peripheral volume of distribution in a three-compartment model.
#' @param input.vp23cmpt Initial estimate for the second peripheral volume of distribution in a three-compartment model.
#' @param input.q3cmpt Initial estimate for first inter-compartmental clearance in a three-compartment model.
#' @param input.q23cmpt Initial estimate for second inter-compartmental clearance in a three-compartment model.
#' @param input.add Initial estimate for additive error.
#' @return An object of class 'nlmixr2' containing the results of the estimation.
#' @import nlmixr2
#' @examples
#' dat <- Oral_2CPT
#' Fit_3cmpt_oral(dat=dat,est.method="nls",input.ka=1,input.cl=4,input.vc3cmpt=70,input.vp3cmpt=35,input.vp23cmpt=35,input.q3cmpt=10,input.q23cmpt=10,input.add=1)
#' @export

Fit_3cmpt_oral<-function(data,
                         est.method,
                         input.ka,
                         input.cl,
                         input.vc3cmpt,
                         input.vp3cmpt,
                         input.vp23cmpt,
                         input.q3cmpt,
                         input.q23cmpt,
                         input.add){

  input.ka<<-input.ka
  input.cl<<-input.cl
  input.vc3cmpt <<- input.vc3cmpt
  input.vp3cmpt<<- input.vp3cmpt
  input.vp23cmpt<<- input.vp23cmpt
  input.q3cmpt<<-  input.q3cmpt
  input.q23cmpt<<-  input.q23cmpt
  input.add <<- input.add

  oral3 <- function() {
    ini({
      tka <-round(log(input.ka), 2) # Absorption rate
      tcl <- round(log(input.cl),2) # Clearance
      tv1  <- round(log(input.vc3cmpt),2) # Central Volume of Distribution
      tv2  <- round(log(input.vp3cmpt),2) # First peripheral volume of distribution
      tv3  <- round(log(input.vp23cmpt),2) # Second peripheral volume of distribution
      tq <- round(log(input.q3cmpt),2) # First inter-compartmental clearance
      tq2 <- round(log(input.q23cmpt),2) # Second inter-compartmental clearance
      add.err <- input.add
    })
    model({
      ka<- exp(tka)
      cl <- exp(tcl)
      v1  <- exp(tv1)
      v2  <- exp(tv2)
      v3  <- exp(tv3)
      q  <- exp(tq)
      q2  <- exp(tq2)
      k  <- cl / v1
      k12<-q/v1
      k21<-q/v2
      k13<-q2/v1
      k31<-q2/v3
      d / dt(depot) = - ka*depot
      d/dt(A1) = - k * A1   - k12 * A1 + k21*A2 +  ka*depot
      d/dt(A2) =   k12 * A1 -  k21 * A2
      d/dt(A3) =   k13 * A1 -  k31 * A3
      cp = A1 / v1
      cp ~ add(add.err)
    })
  }

  maxSSv=100

  if (est.method=="rxSolve"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral3, data =  data,est= "rxSolve",control = rxControl(maxSS=maxSSv))))
  }

  if (est.method=="nls"){
    fit.3cmpt.lst  <- suppressMessages(suppressWarnings( nlmixr2( object = oral3, data =  data,est= "nls",control = nlsControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlm"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral3, data =  data,est= "nlm",control = nlmControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="nlminb"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral3, data =  data,est= "nlminb",control = nlminbControl(rxControl = rxControl(maxSS=maxSSv)))))
  }

  if (est.method=="focei"){
    fit.3cmpt.lst  <-  suppressMessages(suppressWarnings(nlmixr2( object = oral3, data =  data,est= "focei",control = foceiControl(rxControl = rxControl(maxSS=maxSSv)))))
  }


  return(fit.3cmpt.lst)
}
