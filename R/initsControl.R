#' Initialise control settings for pipeline and NPD compartmental analysis
#'
#' Initialise the control settings for pipeline calculation and NPD compartmental analysis, including options for half-life, number of points used for regression, trapezoidal rule method for area under the curve (AUC) calculation, and naive pooled data compartmental analysis (NPD-NCA) strategies.
#'
#' @param half_life Numeric, default is \code{NA}. User-defined half-life value used as a reference for single-point calculations.
#' @param nlastpoints Integer, default is \code{4}. Number of last points used for linear regression during the terminal elimination phase slope estimation.
#' @param trapezoidal.rule Integer, default is \code{1}. Method for applying the trapezoidal rule to calculate the area under the curve (AUC):
#' \itemize{
#'   \item \code{1}: Linear trapezoidal rule.
#'   \item \code{2}: Linear-up log-down trapezoidal rule.
#' }
#' @param nbins Integer, default is \code{8}. Number of time windows for quantile-based partitioning of the time variable in the dataset.
#' @param est.method Character, default is \code{"nls"}. The method used for PK parameter estimation, such as naive pooled data compartmental analysis (NPD-NCA).
#' @param selection.criteria Character, default is \code{"APE"}. Criteria used for evaluating and selecting appropriate parameter values.
#' @param npdcmpt.inits.strategy Integer, default is \code{0}. Strategy for initializing the naive pooled data compartmental analysis:
#' \itemize{
#'   \item \code{0}: Set as \code{1}.
#'   \item \code{1}: Set as pipeline-recommended values.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{half_life}: The user-defined half-life value.
#'   \item \code{nlastpoints}: The number of last points for regression.
#'   \item \code{trapezoidal.rule}: The selected trapezoidal rule method (1 for linear, 2 for linear-up log-down).
#'   \item \code{nbins}: The number of time windows.
#'   \item \code{est.method}: The method used for PK parameter estimation.
#'   \item \code{selection.criteria}: The selection criteria for parameter values.
#'   \item \code{npdcmpt.inits.strategy}: The strategy for initializing the NPD compartment analysis.
#'   \item \code{trapezoidal.rule.c}: A descriptive string for the selected trapezoidal rule.
#'   \item \code{npdcmpt.inits.strategy.c}: A descriptive string for the selected NPD initialisation strategy.
#' }
#'
#' @examples
#' settings <- initsControl(half_life = 5, nlastpoints = 4)
#' print(settings)
#'
#' @export
#'
initsControl<-function( half_life = NA,
                        nlastpoints = 4,
                        trapezoidal.rule = 1,
                        nbins = 8,
                        est.method = "nls",
                        selection.criteria="APE",
                        npdcmpt.inits.strategy= 0){

  if ( trapezoidal.rule==1){
  trapezoidal.rule.c<-"Linear trapezoidal rule "
  }

  if ( trapezoidal.rule==2){
    trapezoidal.rule.c<-"Linear-up log-down trapezoidal rule "
  }

  if (npdcmpt.inits.strategy==0){
   npdcmpt.inits.strategy.c<-"Set as 1"
  }


  if (npdcmpt.inits.strategy==1){
    npdcmpt.inits.strategy.c<-"Set as pipeline recommended values"
  }

  return(  list(half_life =  half_life ,
                nlastpoints =   nlastpoints,
                trapezoidal.rule = trapezoidal.rule,
                nbins = nbins,
                est.method=   est.method,
                selection.criteria=   selection.criteria,
                npdcmpt.inits.strategy=  npdcmpt.inits.strategy,
                trapezoidal.rule.c = trapezoidal.rule.c,
                npdcmpt.inits.strategy.c =  npdcmpt.inits.strategy.c))
}



