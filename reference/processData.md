# Process time–concentration dataset for pharmacokinetic analysis

Processes a time–concentration dataset to derive analysis-ready
variables and structured output for pharmacokinetic evaluation.

## Usage

``` r
processData(dat, verbose = TRUE)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format. The following columns are required (case-insensitive)
  and must be present:

  ID

  :   Subject identifier (required)

  TIME

  :   Nominal or actual time after dose (required)

  DV

  :   Observed concentration (dependent variable) (required)

  EVID

  :   Event ID indicating observation (0) or dosing event (1) (required)

  AMT

  :   Dose amount for dosing records (required)

  RATE

  :   Infusion rate (optional)

  DUR

  :   Infusion duration (optional)

  MDV

  :   Missing dependent variable flag (optional)

  CMT

  :   Compartment number (optional)

  ADDL

  :   Number of additional doses (optional)

  II

  :   Interdose interval (optional)

  SS

  :   Steady-state indicator (optional)

  CENS

  :   Censoring indicator (optional)

- verbose:

  Logical (default = TRUE). When TRUE, the function prints detailed
  processing messages and summary tables to the console, including notes
  on data cleaning and event handling. When FALSE, these messages are
  suppressed and only the returned list is produced.

## Value

A list with two elements:

- dat: A data frame containing the processed time–concentration dataset
  with standardized and derived pharmacokinetic variables, including
  resetflag, SSflag, route, dose_number, DVstd, indiv_lambda_z_eligible,
  and others.

- Datainfo: A data frame summarizing the dataset structure, including
  subject counts and observation counts for first-dose and multiple-dose
  conditions, with contextual notes.

## Details

This function standardizes and preprocesses time–concentration data to
ensure compatibility with pharmacokinetic modeling workflows in nlmixr2.
The operations follow these steps:

1.  Standardize data

    - convert column names to uppercase

    - coerce key columns (TIME, DV, EVID, AMT, RATE, etc.) to numeric

2.  Process events and observations

    - impute EVID from MDV if missing

    - handle censored data (CENS) by converting them to excluded records

    - remove or recode invalid EVID values (e.g., DV = 0 observations)

3.  Expand dose events

    - expand dosing records using nmpkconvert() when ADDL and II are
      present

    - assign dose occasions using mark_dose_number()

4.  Determine administration route and infusion logic

    - derive RATE and DUR when needed

    - identify route (bolus, infusion, oral) based on compartment and
      rate

5.  Generate derived variables

    - calculate time after dose using calculate_tad()

    - compute dose-normalized concentration (DVstd)

    - flag eligible records for terminal elimination phase

6.  Summarize dataset

    - classify dataset as first-dose, repeated-dose, or mixed

    - generate summary metrics for nlmixr2 analysis

    Classification of dosing context is based on pharmacokinetic
    observation records (EVID equal to 0), determining whether observed
    concentrations occur after the first dose, during repeated dosing,
    or across both contexts. The categories are:

    - first_dose: observations occur only after the initial
      administration, without repeated-dose or steady-state intervals.

    - repeated_doses: observations occur only after multiple
      administrations or under steady-state conditions.

    - combined_doses: observations include both first-dose and
      repeated-dose intervals and are analyzed together.

## See also

[`nmpkconvert`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/nmpkconvert.md),
[`mark_dose_number`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/mark_dose_number.md),
[`calculate_tad`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/calculate_tad.md)

## Author

Zhonghui Huang

## Examples

``` r
dat <- Bolus_1CPT
processData(dat)
#> 
#> 
#> Infometrics                               Value          
#> ----------------------------------------  ---------------
#> Dose Route                                bolus          
#> Dose Type                                 combined_doses 
#> Number of Subjects                        120            
#> Number of Observations                    6951           
#> Subjects with First-Dose Interval Data    120            
#> Observations in the First-Dose Interval   2276           
#> Subjects with Multiple-Dose Data          120            
#> Observations after Multiple Doses         4675           
#> ----------------------------------------  ------
#> $dat
#> # A tibble: 7,911 × 28
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1  0       0   0        1 60000     1    60000  65.1  4.07    99     0
#>  2     1  0.25 1126.  7.03     0     0     0    60000  65.1  4.07    99     0
#>  3     1  0.5   870.  6.77     0     0     0    60000  65.1  4.07    99     0
#>  4     1  0.75  884.  6.78     0     0     0    60000  65.1  4.07    99     0
#>  5     1  1    1244   7.13     0     0     0    60000  65.1  4.07    99     0
#>  6     1  1.5   995.  6.90     0     0     0    60000  65.1  4.07    99     0
#>  7     1  2     946.  6.85     0     0     0    60000  65.1  4.07    99     0
#>  8     1  2.5   589.  6.38     0     0     0    60000  65.1  4.07    99     0
#>  9     1  3     754.  6.63     0     0     0    60000  65.1  4.07    99     0
#> 10     1  4    1061.  6.97     0     0     0    60000  65.1  4.07    99     0
#> # ℹ 7,901 more rows
#> # ℹ 16 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>
#> 
#> $Datainfo
#>                               Infometrics          Value
#> 1                              Dose Route          bolus
#> 2                               Dose Type combined_doses
#> 3                      Number of Subjects            120
#> 4                  Number of Observations           6951
#> 5  Subjects with First-Dose Interval Data            120
#> 6 Observations in the First-Dose Interval           2276
#> 7        Subjects with Multiple-Dose Data            120
#> 8       Observations after Multiple Doses           4675
#> 
```
