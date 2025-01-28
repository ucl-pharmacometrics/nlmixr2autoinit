#' Process Pharmacokinetic Dataset for Analysis
#'
#' Processes a pharmacokinetic dataset to prepare it for analysis.
#'
#' @param dat A data frame containing pharmacokinetic data. Expected columns include:
#'   - `ID`: Identifier for each subject.
#'   - `DV`: Dependent variable (e.g., concentration).
#'   - `MDV`: Column indicating missing dependent variable values, used to infer `EVID` if absent.
#'   - `DOSE` or `dose`: Dose amount administered.
#'   - `DUR`: Duration of infusion (optional).
#'   - `RATE`: Infusion rate (optional).
#'   - `SS`: Steady-state indicator (optional).
#'   - `CMT`: Compartment for dosing or sampling (optional).
#'
#' @return A list containing three data frames:
#'   \itemize{
#'     \item `dat`: The complete processed dataset with additional columns for `resetflag`, `SSflag`,
#'       `route`, `dose_number`, and `DVstd` (normalized concentration).
#'     \item `fd_data`: Data between the first and second doses.
#'     \item `md_data`: Data after at least two doses have been administered.
#'   }
#'
#' @details
#' The function performs the following processing steps:
#'   - **Column Verification**: Checks for essential columns (`EVID`, `MDV`, `DOSE`, `RATE`, `SS`, `CMT`)
#'     and imputes `EVID` based on `MDV` if necessary. Issues warnings if `EVID=2` rows are found
#'     and removes these rows.
#'   - **Dosing and Steady-State Flags**: Calculates `resetflag` for each `ID` group based on `EVID=4`
#'     occurrences and sets `SSflag` based on steady-state dosing intervals within each `ID`.
#'   - **Route Identification**: Identifies dosing route as "oral", "infusion", or "bolus" based on
#'     `CMT` and `RATE`.
#'   - **Dataset Splits**: Creates subsets for first-dose interval data (`fd_data`) and multiple-dose
#'     data (`md_data`), normalizing `DV` by dose in each subset.
#'   - **Summary Table**: Generates a summary of the dataset, including total subjects, observations,
#'     and unique counts in each subset. Outputs the summary table and any warnings as console messages.
#'
#' @examples
#' # Example usage:
#' dat <- Bolus_1CPT
#' processed_data <- processData(dat)
#' # Access processed full data, first-dose interval, and multiple-dose data
#' dat <- processed_data$dat
#' Datainfo<-processed_data$Datainfo
#'
#' dat <- Bolus_1CPT
#' dat[1,]$SS=1
#' dat[1,]$II=24
#' processed_data <- processData(dat)
#' dat <- processed_data$dat
#'
#' dat <-theo_sd
#' processed_data <- processData(dat)
#' dat <- processed_data$dat
#'
#' @importFrom dplyr group_by mutate ungroup filter select arrange case_when
#' @importFrom knitr kable
#' @importFrom crayon magenta red black
#' @importFrom stats lm
#' @export
#'
processData<-function(dat){

column_names <- toupper(colnames(dat))
colnames(dat) <- toupper(colnames(dat))

# Convert the values in key columns of input data frame to numeric format
column.list<-c("TIME","DV","MDV","EVID","RATE","DUR", "AMT","ADDL","II")

for (testcolumn in colnames(dat)) {
  if (testcolumn %in% column.list) {
    dat[[testcolumn]] <- as.numeric(dat[[testcolumn]])
  }
}

evid_message1<-NULL
evid_message2<-NULL


if ("CENS" %in% colnames(dat)) {
  if (max(dat$CENS)==1){
  message(black(paste("Warning: 'cens' column is detected, the current version will ignore censored data for initial estimate analysis.")))
    dat[dat$CENS==1,]$EVID=2
   dat$CENS<-0
  }
}

if (!"EVID" %in% colnames(dat)) {

  if (!"MDV" %in% colnames(dat)) {
    evid_message1<-stop(red("Error, no EVID or MDV column found"))
  }

  if ("MDV" %in% colnames(dat)) {
    evid_message1<- paste0("Note: No EVID column detected in the dataset. Impute EVID based on the MDV column.")

    dat$EVID<-0
    dat[dat$MDV==1 & !is.na(dat$DV),]$EVID=1
    dat[dat$MDV==0,]$EVID= 0

  }
}

# Check if other format of EVID
  if (101 %in% dat$EVID) {
  dat[dat$EVID==101,]$EVID<-1
  }

# Check if any row contains EVID=2 and issue a warning if found
if (2 %in% dat$EVID) {
  evid_message2<-paste0("Note: EVID=2 found in the dataset. Rows with EVID=2 are removed for the entire analysis.")
  dat <- dat %>% filter(EVID != 2)
}

if ("DOSE" %in% column_names) {
  dat$DOSE_PRE <- dat$DOSE
  dat$dose<-NULL
}

if ("DUR" %in% colnames(dat)) {
  if (!"RATE" %in% colnames(dat)) {
    dat$RATE <- 0
    dat[dat$DUR > 0,]$RATE <-
      dat[dat$DUR > 0,]$AMT / dat[dat$DUR > 0,]$DUR
  }

  else{
    if (min(dat$RATE)==-2){
      dat$RATE <- 0
      dat[dat$DUR > 0,]$RATE <-
        dat[dat$DUR > 0,]$AMT / dat[dat$DUR > 0,]$DUR
      dat <- dat %>%
        select(-DUR)
    }
  }
}

if (!"DUR" %in% colnames(dat)) {
  if (!"RATE" %in% colnames(dat)) {
    dat$RATE <- 0
  }
}

if (!"SS" %in% colnames(dat)) {
    dat$SS<- 0
}

# CMT

if (!"CMT" %in% colnames(dat)) {
  dat$CMT<- 1
}

if ("CMT" %in% column_names) {
  if (length(unique(dat$CMT)) == 1) {
      if (is.character(dat$CMT)){
        dat$CMT= 1
        message(black(
          paste0(
            "CMT value is a number and has been set to 1 by default."
          )
        ))
      }
  }

  if (length(unique(dat$CMT)) == 2) {
    message(black(
      paste0(
        "Administration site detected to differ from measurement site; extravascular (oral) administration assumed."
      )
    ))

    if (is.character(dat$CMT)){

      cpcmptname<-dat[dat$EVID==0,]$CMT[1]
      dat[dat$EVID==0,]$CMT=2
      # iv case
      if (nrow( dat[dat$EVID %in% c(1,4) & dat$CMT== cpcmptname,])>0){
      dat[dat$EVID %in% c(1,4) & dat$CMT== cpcmptname,]$CMT=2
      }

      dat[dat$EVID %in% c(1,4) & dat$CMT!= cpcmptname,]$CMT=1

      message(black(
        paste0(
          "CMT value is NOT a number and has been set to 1 (depot) 2 (centre) by default."
        )
      ))
    }
  }

  if (length(unique(dat$CMT)) > 2) {
    message(black(
      paste0(
        "Administration site detected to differ from measurement site; extravascular (oral) administration assumed."
      )
    ))

    if (is.character(dat$CMT)){
       stop("Error: Only cases with CMT values of up to two are supported currently")
    }
  }
}

# Analyse EVID=4 reset situation and add an reseflag (similar with rxode2 resetno)
# Apply the resetflag only within each ID group
# Initialise the resetflag column with 1
dat$resetflag <- 1

# Loop through each unique ID to apply the resetflag condition separately within each group
for (id in unique(dat$ID)) {
  # Subset rows for the current ID
  id_rows <- which(dat$ID == id)

  for (i in id_rows[-1]) {  # Start from the second row within the ID subset
    # If EVID equals 4, increment resetflag by 1 from the previous row's value
    if (dat$EVID[i] == 4) {
      dat$resetflag[i] <- dat$resetflag[i - 1] + 1
    } else {
      # Otherwise, keep resetflag the same as the previous row's value
      dat$resetflag[i] <- dat$resetflag[i - 1]
    }
  }
}

# Apply the SSflag condition only within each ID group
# Initialize the SSflag column with 0
dat$SSflag <- 0

for (id in unique(dat$ID)) {
  # Identify rows for the current ID
  id_rows <- which(dat$ID == id)

  for (i in id_rows) {
    # Check if SS is 1 in the current row
    if (dat$SS[i] == 1) {
      # Define the end time by adding hours_to_flag to the current TIME
      end_time <- dat$TIME[i] + dat$II[i]

      # Find the rows within the same ID that fall within the end_time
      rows_to_flag <- id_rows[dat$TIME[id_rows] <= end_time & dat$TIME[id_rows] >= dat$TIME[i]]

      # Set SSflag to 1 for the specified range within the same ID
      dat$SSflag[rows_to_flag] <- 1
    }
  }
}

# Convert the values in key columns of input data frame to numeric format
column.list<-c("TIME","DV","MDV","EVID","RATE","DUR","AMT","CMT", "ADDL","II")

for (testcolumn in colnames(dat)) {
  if (testcolumn %in% column.list) {
    dat[[testcolumn]] <- as.numeric(dat[[testcolumn]])
  }
}

# Identify route
dat <- dat %>%
  group_by(ID, resetflag) %>%
  mutate(
    # Identify the CMT of EVID=0 within the group (assume first occurrence of EVID=0 is representative)
    evid0_cmt = first(CMT[EVID == 0], default = NA_real_),

    # Determine route based on CMT and RATE
    route = case_when(
      EVID %in% c(1,4) & CMT != evid0_cmt ~ "oral",  # If CMT does not match the CMT of EVID=0, set route to "oral"
      EVID %in% c(1,4) &  CMT == evid0_cmt & RATE > 0 ~ "infusion",  # If CMT matches and RATE > 0, set route to "infusion"
      EVID %in% c(1,4) &  CMT == evid0_cmt & RATE == 0 ~ "bolus",  # If CMT matches and RATE == 0, set route to "bolus"
      EVID == 0 ~ NA # Keep existing value in `route` for other rows
    )
  ) %>%
  ungroup() %>%
  select(-evid0_cmt)  # Remove temporary column


# Convert pharmacokinetic dataset to depreciated format with additional dosing
# Must run before marking dose number
if ("ADDL" %in% column_names) {
  if (sum(dat$ADDL)>0){
  dat <- nmpkconvert(dat)
  }
}

# Mark the dose number
dat <- mark_dose_number(dat)
dat<-calculate_tad(dat)

dat$duration_obs<-0

if (nrow(dat[dat$rateobs!=0,])>0){
  dat$duration_obs <- dat[dat$rateobs!=0,]$dose / dat[dat$rateobs!=0,]$rateobs
}

dat <- dat %>%
  mutate(original_EVID = EVID) %>%  # Create a backup column for EVID
  group_by(ID, dose_number) %>%  # Group by ID and dose_number to work within each dose cycle
  mutate(
    # Calculate Tmax within each ID and dose_number group where EVID == 0
    Tmax = ifelse(any(EVID == 0), TIME[which.max(DV[EVID == 0])], NA_real_),
    # Set EVID to 2 for rows where EVID == 0, TIME > Tmax, and DV == 0 within each ID and dose_number group
    # Only set if Tmax is not NA
    EVID = if_else(EVID == 0 & !is.na(Tmax) & TIME > Tmax & DV == 0, 2, EVID),
    # Set DV to NA for rows where EVID is set to 2
    DV = if_else(EVID == 2, 0, DV)
  ) %>%
  ungroup()

dat <- select(dat, -Tmax)
dat<-as.data.frame(dat)


# Check if any EVID was changed from 0 to 2
if (any(dat$original_EVID == 0 & dat$EVID == 2)) {
  message(black("Some rows have been updated: EVID set to 2 for DV=0 after Tmax, and they will be removed for the entire analysis."))
  # Check if any row contains EVID=2 and issue a warning if found
  if (2 %in% dat$EVID) {
    dat <- dat %>% filter(EVID != 2)
  }
}

dat$DVstd <- dat$DV / dat$dose

############################# Group data####################################
# First dose data (Data between the first and second doses,dose number=1)
fd_data<-dat[dat$dose_number==1 & dat$iiobs==0,]
# Normalise concentration by dose
fd_data$DVstd <- fd_data$DV / fd_data$dose
fd_data_obs<-dat[dat$dose_number==1 & dat$EVID==0& dat$iiobs==0,]

# Non-first dose data (Data after at least two doses have been administered,dose number>1)
md_data1 <-dat[dat$dose_number>1,]
md_data2 <-dat[dat$dose_number==1 & dat$iiobs>0,]
md_data<-rbind(md_data1,md_data2)
md_data<-md_data[with(md_data, order(ID, resetflag, TIME, -AMT)), ]

md_data$DVstd <-md_data$DV / md_data$dose
md_data_obs <-dat[dat$dose_number>1 & dat$EVID==0,]

####################Summary output#############################

# total subjects and observations for each dataset
nids <- nrow(dat[!duplicated(dat$ID),])
nobs <- nrow(dat[dat$EVID == 0,])

nids_fd <- nrow(fd_data_obs[!duplicated(fd_data_obs$ID),])   # Number of unique IDs in FD
nobs_fd <- nrow(fd_data_obs[fd_data_obs$EVID == 0,])         # Number of observations in FD

nids_md <- nrow(md_data_obs[!duplicated(md_data_obs$ID),])   # Number of unique IDs in MD
nobs_md <- nrow(md_data_obs[md_data_obs$EVID == 0,])         # Number of observations in MD


summary_doseinfo <- data.frame(
  Infometrics = c("Dose Route",
                  "Total Number of Subjects",
                  "Total Number of Observations",
                  "Subjects with First-Dose Interval Data",
                  "Observations in the First-Dose Interval",
                  "Subjects with Multiple-Dose Data",
                  "Observations after Multiple Doses"),
  Values = c(unique(dat[dat$EVID %in% c(1,4),]$route),
            nids,
            nobs,
            nids_fd,
            nobs_fd,
            nids_md,
            nobs_md)
)


# message(magenta(
#   paste(capture.output(knitr::kable(summary_doseinfo, format = "simple")), collapse = "\n")
# ))


# Create the formatted table output
summary_doseinfo_output <- paste(capture.output(knitr::kable(summary_doseinfo, format = "simple")), collapse = "\n")
complete_output <- paste(summary_doseinfo_output, sep = "\n")

if (!is.null(evid_message1)||!is.null(evid_message2)){
# Define the footnote text
footnote_evid<- paste(evid_message1, evid_message2, sep = "\n")
# Combine table output and footnote, separated by a line for clarity
complete_output <- paste(summary_doseinfo_output, footnote_evid, sep = "\n")
}
# Display the complete output with footnote in the console
message(magenta(complete_output))

message(magenta(
  paste0("------------------------------------  ------")
))

return(list(dat=dat,Datainfo =complete_output ))
}






