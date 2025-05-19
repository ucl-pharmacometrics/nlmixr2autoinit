#' Process Pharmacokinetic Dataset for Analysis
#'
#' Processes a pharmacokinetic (PK) dataset to derive analysis-ready variables and structure.
#'
#' @param dat A data frame containing raw pharmacokinetic data. Expected (case-insensitive) columns include:
#'   \describe{
#'     \item{ID}{Subject identifier}
#'     \item{TIME}{Time after dose}
#'     \item{DV}{Dependent variable (e.g., plasma concentration)}
#'     \item{MDV}{Missing dependent variable indicator}
#'     \item{EVID}{Event identifier (e.g., dose, observation)}
#'     \item{AMT}{Dose amount}
#'     \item{RATE}{Infusion rate (optional)}
#'     \item{DUR}{Duration of infusion (optional)}
#'     \item{ADDL}{Additional doses (optional)}
#'     \item{II}{Interdose interval (optional)}
#'     \item{SS}{Steady-state indicator (optional)}
#'     \item{CMT}{Compartment (optional)}
#'     \item{CENS}{Censoring indicator (optional)}
#'   }
#'
#' @return A list with two components:
#' \describe{
#'   \item{dat}{Processed full dataset with standardized variables and derived columns:
#'     \code{resetflag}, \code{SSflag}, \code{route}, \code{dose_number}, \code{DVstd},
#'     \code{indiv_lambda_z_eligible}, and others.}
#'   \item{Datainfo}{A formatted summary table of dataset structure, subject counts, and observation counts for
#'     first-dose and multiple-dose conditions, with contextual notes.}
#' }
#'
#' @details
#' This function performs several operations critical for PK data preprocessing:
#' \enumerate{
#'   \item \strong{Standardization}: Converts column names to uppercase; coerces numeric columns to appropriate types.
#'   \item \strong{Event Processing}: Imputes \code{EVID} from \code{MDV} if missing; handles censored observations; filters invalid \code{EVID=2} records.
#'   \item \strong{Flag Generation}:
#'     \itemize{
#'       \item \code{resetflag}: Identifies segments between distinct dose events or steady-state switches.
#'       \item \code{SSflag}: Flags rows within steady-state observation windows.
#'     }
#'   \item \strong{Infusion and Route Logic}: Calculates \code{RATE} from \code{AMT} and \code{DUR} if needed; classifies administration route (oral, infusion, bolus).
#'   \item \strong{Compartment Checks}: Ensures a maximum of two unique \code{CMT} values per subject and dosing interval.
#'   \item \strong{Dose Expansion}: Expands rows using \code{ADDL}/\code{II} if present via \code{nmpkconvert()}.
#'   \item \strong{Derived Metrics}:
#'     \itemize{
#'       \item \code{dose_number}: Count of dose administrations per subject.
#'       \item \code{tad}: Time after last dose.
#'       \item \code{DVstd}: Normalized dependent variable (DV/dose).
#'       \item \code{indiv_lambda_z_eligible}: Eligibility flag for elimination phase based on route and time of maximum DV.
#'     }
#'   \item \strong{Summary Generation}: Constructs a summary table of route, subjects, and observations across first- and multiple-dose data.
#' }
#'
#' @examples
#'
#' dat <- Bolus_1CPT
#' results <- processData(dat)
#' head(results$dat)
#' cat(results$Datainfo)
#'
#' @export

processData<-function(dat){

#-------------- STEP 1: Data Standardization -------------------------#
  # Standardize column names to uppercase
  column_names <- toupper(colnames(dat))
  colnames(dat) <- toupper(colnames(dat))

  # Convert key columns to numeric format
  column.list <- c("TIME",
                "DV",
                "MDV",
                "EVID",
                "RATE",
                "DUR",
                "AMT",
                "ADDL",
                "II",
                "SS",
                "CMT")

  for (testcolumn in colnames(dat)) {
    if (testcolumn %in% column.list) {
      dat[[testcolumn]] <- as.numeric(dat[[testcolumn]])
    }
  }

  #-------------- STEP 2: Event Flag Processing -------------------------#
  # Initialize message container
  evid_messages <- character()

  # Censored data handling
  if ("CENS" %in% colnames(dat)) {
    if (max(dat$CENS, na.rm = TRUE) == 1) {
      msg <- "Note: Censored data (CENS=1) excluded by setting EVID=2"
      evid_messages <- c(evid_messages, msg)
      # message(crayon::black(msg))

      dat <- dat %>%
        dplyr::mutate(
          EVID = dplyr::if_else(.data$CENS == 1, 2L, .data$EVID),
          CENS = 0L
        )
    }
  }

  # EVID/MDV handling
  if (!"EVID" %in% colnames(dat)) {
    if (!"MDV" %in% colnames(dat)) {
      stop(crayon::red("Error: Missing both EVID and MDV columns"))
    }

    msg <- "Note: Imputed EVID from MDV (MDV=1 to EVID=1)"
    evid_messages <- c(evid_messages, msg)
    # message(crayon::black(msg))

    dat <- dat %>%
      dplyr::mutate(
        EVID = dplyr::case_when(
          .data$MDV == 1 & !is.na(.data$DV) ~ 1L,
          .data$MDV == 0 ~ 0L,
          TRUE ~ 0L
        )
      )
  }

  # Clean EVID values
  dat <- dat %>%
    # Convert EVID=101 to 1 with logging
    dplyr::mutate(
      EVID = dplyr::if_else(.data$EVID == 101, 1L, .data$EVID)
    ) %>% {
      n_101 <- sum(.$EVID == 101, na.rm = TRUE)
      if(n_101 > 0) {
        msg <- paste("Converted", n_101, "EVID=101 records to EVID=1")
        evid_messages <<- c(evid_messages, msg)
        # message(crayon::black(msg))
      }
      .
    }

  #Reset Flag Calculation #
  dat <- dat %>%
    dplyr::group_by(ID) %>%
    # Create SS column if missing and initialize to 0
    { if (!"SS" %in% names(.)) dplyr::mutate(., SS = 0) else . } %>%
    dplyr::mutate(
      # Identify reset events (EVID=4)
      reset_event = (.data$EVID == 4 | .data$SS == 1),

      # Calculate resetflag using cumulative sum
      resetflag = 1L + cumsum(reset_event)
    ) %>%
    dplyr::select(-reset_event) %>%
    dplyr::ungroup()


  # Handle LLOQ (DV=0)
  # remove raw_EVID if it exists
  if ("raw_evid" %in% colnames(dat)) {
    dat$raw_evid <- NULL
  }

  dat <- dat %>%
    dplyr::mutate(
      raw_EVID = EVID,  # Backup original EVID
      EVID = dplyr::if_else(
        condition = .data$raw_EVID == 0 & .data$DV == 0,
        true = 2L,
        false = .data$EVID
      )
    ) %>% {
      # Count rows where EVID was changed from 0 to 2
      n_converted <- sum(.$raw_EVID == 0 & .$EVID == 2, na.rm = TRUE)
      if (n_converted > 0) {
        msg <- paste("Converted", n_converted, "rows with DV=0 from EVID=0 to EVID=2. These are excluded.")
        evid_messages <<- c(evid_messages, msg)
        # message(crayon::black(msg))
      }
      .  # Return data for further piping
    }


  # Remove EVID=2 with logging
  dat <- dat %>%
    dplyr::filter(.data$EVID != 2) %>% {
      n_removed <- sum(.$EVID == 2, na.rm = TRUE)
      if(n_removed > 0) {
        msg <- paste("Removed", n_removed, "EVID=2 records (DV=0)")
        evid_messages <<- c(evid_messages, msg)
        # message(crayon::black(msg))
      }
      .
    }

  # Final message output
  if(length(evid_messages) > 0) {
    message(crayon::black(paste(evid_messages, collapse = "\n")))
  }


 #-------------- STEP 3: DUR/RATE/SS Processing -----------------------#

  if ("DUR" %in% colnames(dat)) {
    # Scenario 1: RATE Column Missing
    if (!"RATE" %in% colnames(dat)) {
      dat <- dat %>%
        dplyr::mutate(RATE = dplyr::case_when(# Calculate rate for positive duration values
          .data$DUR > 0 ~ .data$AMT / .data$DUR,
          # Default to zero for non-positive durations
          TRUE ~ 0))
    }

    # Scenario 2: RATE Exists with Special Flag Value
    else {
      # Check if RATE needs reset (using -2 as indicator)
      if (min(dat$RATE, na.rm = TRUE) == -2) {
        dat <- dat %>%
          dplyr::mutate(# Recalculate rate using same logic as Scenario 1
            RATE = dplyr::case_when(.data$DUR > 0 ~ .data$AMT / .data$DUR,
                                    TRUE ~ 0)) %>%
          # Remove DUR column after recalculation
          dplyr::select(-"DUR")
      }
    }
  }

  # Handle mandatory RATE column creation
  if (!"DUR" %in% colnames(dat) && !"RATE" %in% colnames(dat)) {
    dat <- dat %>%
      dplyr::mutate(RATE = 0)  # Initialize with numeric zeros
  }

#-------------- STEP 4: Reset and SS Flags ---------------------------#
  dat <- dat %>%
    # Create II column if missing and initialize to 0
    { if (!"II" %in% names(.)) dplyr::mutate(., II = 0) else . } %>%
    # Initialize SSflag column with 0 for all rows
    dplyr::mutate(SSflag = 0) %>%
    # Process within each subject ID
    dplyr::group_by(ID) %>%
    dplyr::mutate(
      # Create list of SS event windows: start/end times for SS=1 events
      SS_events = list(
        data.frame(
          start = TIME[SS == 1],
          end = TIME[SS == 1] + II[SS == 1]
        )
      )
    ) %>%
    # Check each row individually
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Determine if current time falls within any SS event window
      SSflag = as.integer(
        any(
          TIME >= SS_events$start & TIME <= SS_events$end
        )
      )
    ) %>%
    dplyr::ungroup() %>%
    # Remove temporary event tracking column
    dplyr::select(-SS_events)


#------------- STEP 4: Compartment Handling & Administration Route ----------#
# Ensure CMT column exists with default value 1
for (testcolumn in colnames(dat)) {
    if (testcolumn %in% column.list) {
      dat[[testcolumn]] <- as.numeric(dat[[testcolumn]])
    }
  }


# Process compartment logic per ID/reset group
dat <- dat %>%
  # Ensure CMT exists: if missing, initialize with 1
  { if (!"CMT" %in% colnames(.)) dplyr::mutate(., CMT = 1) else . } %>%
  dplyr::group_by(ID, resetflag) %>%
  dplyr::mutate(
    # Get observation compartment (EVID=0) for current reset group
    obs_cmt = dplyr::first(CMT[EVID == 0], na_rm = TRUE),

    # Count unique compartments in current group
    n_unique_cmt = dplyr::n_distinct(CMT)
  ) %>%
  dplyr::ungroup() %>%
  # Route determination logic
  dplyr::mutate(
    route = dplyr::case_when(
      # Single compartment logic
      n_unique_cmt == 1 &
        EVID %in% c(1, 4) & RATE > 0 ~ "infusion",
      n_unique_cmt == 1 & EVID %in% c(1, 4) & RATE == 0 ~ "bolus",

      # Dual compartment logic
      n_unique_cmt == 2 &
        EVID %in% c(1, 4) & CMT != obs_cmt ~ "oral",
      n_unique_cmt == 2 &
        EVID %in% c(1, 4) & CMT == obs_cmt & RATE > 0 ~ "infusion",
      n_unique_cmt == 2 &
        EVID %in% c(1, 4) & CMT == obs_cmt & RATE == 0 ~ "bolus",

      # Default case
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-obs_cmt,-n_unique_cmt)

# Error handling for unsupported configurations
invalid_groups <- dat %>%
  dplyr::group_by(ID, resetflag) %>%
  dplyr::filter(dplyr::n_distinct(CMT) > 2) %>%
  dplyr::ungroup()

if (nrow(invalid_groups) > 0) {
  stop(
    paste(
      "Package requires ≤2 distinct CMT values. Found",
      length(unique(invalid_groups$ID)),
      "subjects with >2 CMT values",
      "(not suitable for metabolites/multi-site analysis)."
    )
  )
}

#------- STEP 5: Dose/tad/DVstd(standardised concentration) processing ------#

# Convert pharmacokinetic dataset to depreciated format with additional dosing
if ("ADDL" %in% column_names) {
  if (!is.numeric(dat$ADDL)) {
    stop("ADDL column must be numeric. Found type: ", class(dat$ADDL))
  }
  if (sum(dat$ADDL)>0){
  dat <- nmpkconvert(dat)
  }
}

# Mark the dose number
dat <- mark_dose_number(dat)

# Calculate time after the last dose (tad)
dat<-calculate_tad(dat)

# The duration inherited from the most recent dose, applied to observation rows.
dat <- dat %>%
  dplyr::mutate(
    durationobs = dplyr::if_else(
      condition = (rateobs != 0),
      true      = dose / rateobs,
      false     = 0.0,
      missing   = 0.0
    )
  )

dat$DVstd <- dat$DV / dat$dose

#-------------- STEP 6: individual lambda-z eligibility -----------------#
# Flag observations eligible for elimination phase (lambdaz) calculations
dat <- dat %>%
  dplyr::group_by(ID, dose_number) %>%                # Group by subject and dose interval
  # 1. Calculate Tmax (time of maximum observed concentration)
  dplyr::mutate(
    Tmax = dplyr::if_else(
      any(EVID == 0),
      TIME[which.max(dplyr::if_else(EVID == 0, DV, -Inf))],
      NA_real_
    )
  ) %>%
  # 2. Count eligible post-Tmax points based on administration route
  dplyr::mutate(
    n_post_Tmax = dplyr::case_when(
      routeobs %in% c("bolus", "infusion") ~
        sum(EVID == 0 & TIME > Tmax, na.rm = TRUE),   # IV: strict post-Tmax
      routeobs == "oral" ~
        sum(EVID == 0 & TIME >= Tmax, na.rm = TRUE),  # Oral: include Tmax
      TRUE ~ 0L                                        # Fallback counter
    ),
    # 3. Create eligibility flag (>=3 points required)
    indiv_lambdaz_eligible = as.integer(n_post_Tmax >= 3)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-Tmax, -n_post_Tmax)                  # Cleanup temporary columns

#------------------ STEP 7: Summarise data --------------------#

# 1. Process first-dose data (excluding steady-state cases) ----
fd_data <- dat %>%
  dplyr::filter(
    dose_number == 1,
    iiobs == 0,
    SS != 1 | is.na(SS)  # Exclude SS=1 pseudo-first doses
  ) %>%
  dplyr::mutate(DVstd = DV / dose)  # Dose-normalize concentrations

# Observation subset for first-dose data
fd_data_obs <- fd_data %>%
  dplyr::filter(EVID == 0)  # Keep observation records only

# 2. Process multi-dose data (including steady-state cases) ----
md_data <- dat %>%
  dplyr::filter(
    dose_number > 1 |
      (dose_number == 1 & (iiobs > 0 | SS == 1))  # Include SS=1 special case
  ) %>%
  # Order records by:
  dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT)) %>%
  dplyr::mutate(DVstd = DV / dose)  # Dose-normalize concentrations

# Observation subset for multi-dose data (derived from processed data)
md_data_obs <- md_data %>%
  dplyr::filter(EVID == 0)  # Keep observation records only

####################Summary output#############################
# Calculate summary metrics ----
metrics <- list(
  # 1. Total metrics
  total = dat %>%
    dplyr::summarize(
      route = toString(unique(route[EVID %in% c(1, 4)])),  # Collapse multiple routes
      n_ids = dplyr::n_distinct(ID),
      n_obs = sum(EVID == 0)
    ),

  # 2. First-dose metrics
  fd = fd_data_obs %>%
    dplyr::summarize(
      n_ids = dplyr::n_distinct(ID),
      n_obs = dplyr::n()
    ),

  # 3. Multi-dose metrics
  md = md_data_obs %>%
    dplyr::summarize(
      n_ids = dplyr::n_distinct(ID),
      n_obs = dplyr::n()
    )
)

# Determine data type based on dosing structure ----
has_first_dose <-
  any(dat$dose_number == 1 & dat$EVID == 0 & dat$iiobs == 0)
has_repeated_dose <- any(dat$dose_number > 1 & dat$EVID == 0)

dose_type <- dplyr::case_when(
  has_first_dose & !has_repeated_dose ~ "first_dose",
  has_first_dose & has_repeated_dose  ~ "combined_doses",
  TRUE                                ~ "repeated_doses"
)

# Build summary table
summary_doseinfo <- data.frame(
  Infometrics = c("Dose Route",
                  "Dose Type",
                  "Total Number of Subjects",
                  "Total Number of Observations",
                  "Subjects with First-Dose Interval Data",
                  "Observations in the First-Dose Interval",
                  "Subjects with Multiple-Dose Data",
                  "Observations after Multiple Doses"),
  Value = c(
    metrics$total$route,
    dose_type,
    metrics$total$n_ids,
    metrics$total$n_obs,
    metrics$fd$n_ids,
    metrics$fd$n_obs,
    metrics$md$n_ids,
    metrics$md$n_obs
  ),
  stringsAsFactors = FALSE    # Prevent automatic factor conversion
)
# Format output ----
summary_doseinfo_output <-
  paste(capture.output(knitr::kable(summary_doseinfo, format = "simple")),
        collapse = "\n")

# Footnote generation
if (length(evid_messages) > 0) {
  # Format messages as bullet points with proper line wrapping
  formatted_messages <- paste0("* ", evid_messages, collapse = "\n")

  # Construct final output with separation between table and footnotes
  complete_output <- paste(summary_doseinfo_output,
                           # Main summary table
                           "\n\n[Notes]\n",
                           # Section header for validation messages
                           formatted_messages,
                           # Formatted bullet points
                           sep = "")
} else {
  # Maintain clean output when no messages exist
  complete_output <- summary_doseinfo_output
}

# Display the complete output with footnote in the console
message(crayon::magenta(complete_output))

message(crayon::magenta(paste0(
  "----------------------------------------  ------"
)))

return(list(dat=dat,Datainfo =summary_doseinfo ))
}







