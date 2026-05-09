# Perform extended single-point pharmacokinetic calculations

Extended single-point pharmacokinetic analysis for deriving CL, Vd, and
ka

## Usage

``` r
run_single_point_extra(
  dat = NULL,
  half_life = NULL,
  single_point_base.lst = NULL,
  route = c("bolus", "oral", "infusion"),
  dose_type = NULL,
  pooled_ctrl = pooled_control(),
  ssctrl = ss_control()
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- half_life:

  Optional numeric value representing the elimination half-life of the
  drug. If not provided, half-life is estimated within
  [`run_single_point_base()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md)
  using
  [`get_hf()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_hf.md)
  applied to the pharmacokinetic observations.

- single_point_base.lst:

  A list object returned from the base single-point calculation that
  includes input data, preprocessing information, and initial estimates
  of CL and Vd. If not supplied, the function will internally call the
  base routine using dat.

- route:

  Character string specifying the route of administration. Must be one
  of "bolus", "oral", or "infusion".

- dose_type:

  Classified as "first_dose", "repeated_doses", or "combined_doses"
  based on whether observed concentrations occur following the first
  administration, during repeated dosing, or across both contexts. This
  parameter is passed directly to
  [`run_single_point_base()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md).

- pooled_ctrl:

  A list of pooled-analysis control options created using
  [`pooled_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md).
  These control time binning and time-after-dose rounding when pooled
  summaries are required.

- ssctrl:

  A list of steady-state control options created using
  [`ss_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md).
  These govern assumptions and thresholds used when interpreting
  steady-state behavior in single-point logic.

## Value

A list containing:

- singlepoint.results: A data frame with estimated ka, CL, Vd,
  computation start time, processing time, and a descriptive message of
  the applied logic.

- dat: The dataset used for processing.

- single_point.ka.out: Output used for estimating the absorption rate
  constant (for oral administration), if applicable.

- single_point_cl_df: Data used for clearance estimation.

- single_point_vd_df: Data used for volume of distribution estimation.

- approx.vc.out: Data used for estimating the volume of distribution
  from maximum concentration (Cmax) and dose when needed.

## Details

The function derives pharmacokinetic parameters using the following
logic:

- When both clearance (CL) and volume of distribution (Vd) are available
  from the base calculation, these values are directly used without
  modification.

- If Vd is missing but CL and elimination half-life are provided, Vd is
  calculated using: \$\$V_d = \frac{CL \cdot t\_{1/2}}{\ln(2)}\$\$

- If CL is missing but Vd and half-life are available, CL is calculated
  using: \$\$CL = \frac{V_d \cdot \ln(2)}{t\_{1/2}}\$\$

- If both CL and Vd are unavailable but half-life is provided, Vd is
  estimated using dose and Cmax: \$\$V_d =
  \frac{\mathrm{Dose}}{C\_{\mathrm{max}}}\$\$ and CL is subsequently
  derived: \$\$CL = \frac{V_d \cdot \ln(2)}{t\_{1/2}}\$\$

- For oral administration, the absorption rate constant (ka) is
  estimated using a solution-based approach implemented in
  [`run_ka_solution()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_ka_solution.md).
  Only concentration–time data from the absorption phase are used,
  defined as time after dose (tad) less than 20% of the terminal
  half-life and occurring prior to Tmax, where absorption is the
  dominant process.

The function supports bolus, oral, and infusion administration routes.
For oral dosing, only data within the absorption phase are used to
estimate the absorption rate constant (ka), specifically using
concentration-time observations prior to the maximum concentration
(Tmax).

## See also

[`run_single_point_base`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md),
[`run_single_point`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point.md),
[`run_ka_solution`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_ka_solution.md)

## Author

Zhonghui Huang

## Examples

``` r
dat <- Bolus_1CPT
out <- processData(dat)
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
fdat <- out$dat
froute <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
half_life <- get_hf(dat = fdat)$half_life_median
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 11 h
run_single_point_extra(
  dat = fdat,
  half_life = half_life,
  single_point_base.lst = run_single_point_base(
    dat = fdat,
    half_life = half_life,
    route = froute
  ),
  route = froute
)
#> $singlepoint.results
#>   ka   cl   vd           starttime time.spent
#> 1 NA 3.04 67.6 2026-05-09 22:07:59      1.557
#>                                                         single_point.message
#> 1 CL and Vd were calculated directly from steady-state and single-dose data.
#> 
#> $dat
#> # A tibble: 7,911 × 43
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
#> # ℹ 31 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>, SteadyState <lgl>,
#> #   recent_ii <dbl>, doses_req_fixed <dbl>, last_two_doses_interval <dbl>,
#> #   doses_req_hl <dbl>, doses_required <dbl>, dose_count_before_obs <int>, …
#> 
#> $single_point.ka.out
#> [1] NA
#> 
#> $single_point_cl_df
#> # A tibble: 842 × 55
#>       ID  TIME     DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1  144.  287.   5.66     0     0     0    60000  65.1  4.07    99     0
#>  2     1  145  1690.   7.43     0     0     0    60000  65.1  4.07    99     0
#>  3     1  164   243.   5.49     0     0     0    60000  65.1  4.07    99     0
#>  4     1  192.  242.   5.49     0     0     0    60000  65.1  4.07    99     0
#>  5     1  216.  281.   5.64     0     0     0    60000  65.1  4.07    99     0
#>  6     1  216. 1500.   7.31     0     0     0    60000  65.1  4.07     0     0
#>  7     1  240   268.   5.59     0     0     0    60000  65.1  4.07     0     0
#>  8     2  144.   72    4.28     0     0     0    10000  58.6  2.94    99     0
#>  9     2  146.  273.   5.61     0     0     0    10000  58.6  2.94    99     0
#> 10     2  168.   72.2  4.28     0     0     0    10000  58.6  2.94    99     0
#> # ℹ 832 more rows
#> # ℹ 43 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>, SteadyState <lgl>,
#> #   recent_ii <dbl>, doses_req_fixed <dbl>, last_two_doses_interval <dbl>,
#> #   doses_req_hl <dbl>, doses_required <dbl>, dose_count_before_obs <int>, …
#> 
#> $single_point_vd_df
#> # A tibble: 120 × 45
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1  0.25 1126.  7.03     0     0     0    60000  65.1  4.07    99     0
#>  2     2  0.25  220.  5.39     0     0     0    10000  58.6  2.94    99     0
#>  3     3  0.25 1061.  6.97     0     0     0    60000  66.6  3.82    99     0
#>  4     4  0.25  609.  6.41     0     0     0    30000  85.4  3.39    99     0
#>  5     5  0.25 2251.  7.72     0     0     0   120000  55.8  3.45    99     0
#>  6     6  0.25  596.  6.39     0     0     0    30000  52.2  4.26    99     0
#>  7     7  0.25  140.  4.94     0     0     0    10000  66.5  4.16    99     0
#>  8     8  0.25 1282.  7.16     0     0     0   120000  95.1  2.93    99     0
#>  9     9  0.25 2429.  7.80     0     0     0   120000  56.0  6.21    99     0
#> 10    10  0.25  134.  4.90     0     0     0    10000  77.6  7.08    99     0
#> # ℹ 110 more rows
#> # ℹ 33 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>, SteadyState <lgl>,
#> #   recent_ii <dbl>, doses_req_fixed <dbl>, last_two_doses_interval <dbl>,
#> #   doses_req_hl <dbl>, doses_required <dbl>, dose_count_before_obs <int>, …
#> 
#> $approx.vc.out
#> $approx.vc.out$approx.vc.value
#> [1] 48.33771
#> 
#> $approx.vc.out$approx.vc.dat
#>      ID dose_number    vd
#> 1     1           1  40.0
#> 2     2           1  32.4
#> 3     3           1  37.6
#> 4     4           1  45.2
#> 5     5           1  35.6
#> 6     6           1  39.0
#> 7     7           1  37.2
#> 8     8           1  35.0
#> 9     9           1  45.3
#> 10   10           1  50.9
#> 11   12           1  43.5
#> 12   13           1  48.0
#> 13   14           1  55.1
#> 14   15           1  49.3
#> 15   16           1  31.0
#> 16   17           1  41.8
#> 17   18           1  49.1
#> 18   19           1  37.8
#> 19   21           1  69.6
#> 20   23           1  27.7
#> 21   24           1  23.6
#> 22   25           1  42.6
#> 23   26           1  38.8
#> 24   27           1  41.5
#> 25   28           1  30.5
#> 26   29           1  39.6
#> 27   30           1  43.7
#> 28   31           1  38.1
#> 29   33           1  29.6
#> 30   35           1  40.0
#> 31   36           1  34.4
#> 32   38           1  47.2
#> 33   39           1  50.9
#> 34   41           1  42.2
#> 35   42           1  44.3
#> 36   43           1  32.8
#> 37   45           1  34.8
#> 38   46           1  80.4
#> 39   47           1  36.2
#> 40   48           1  56.8
#> 41   49           1  34.3
#> 42   50           1  50.5
#> 43   51           1  48.0
#> 44   52           1  41.5
#> 45   53           1  26.6
#> 46   54           1  31.5
#> 47   56           1  34.2
#> 48   57           1  47.4
#> 49   59           1  31.1
#> 50   61           1  62.1
#> 51   62           1  28.9
#> 52   63           1  29.4
#> 53   64           1  40.6
#> 54   65           1  34.9
#> 55   66           1  31.1
#> 56   67           1  29.8
#> 57   68           1  30.2
#> 58   69           1  21.5
#> 59   70           1  49.8
#> 60   71           1  41.3
#> 61   72           1  32.2
#> 62   73           1  29.8
#> 63   74           1  44.5
#> 64   75           1  40.8
#> 65   76           1  37.4
#> 66   77           1  36.8
#> 67   78           1  37.8
#> 68   79           1  37.5
#> 69   80           1  50.4
#> 70   81           1  26.8
#> 71   82           1  46.4
#> 72   83           1  69.6
#> 73   84           1  72.0
#> 74   85           1  55.2
#> 75   87           1  45.5
#> 76   88           1  35.8
#> 77   89           1  48.5
#> 78   90           1  35.7
#> 79   91           1  19.5
#> 80   92           1  42.3
#> 81   93           1  27.7
#> 82   94           1  37.9
#> 83   95           1  35.6
#> 84   97           1  49.7
#> 85   99           1  52.0
#> 86  100           1  34.7
#> 87  101           1  39.9
#> 88  102           1  49.1
#> 89  103           1  44.1
#> 90  104           1  47.2
#> 91  105           1  57.9
#> 92  107           1  45.8
#> 93  108           1  49.6
#> 94  109           1  47.2
#> 95  110           1  57.3
#> 96  111           1  37.9
#> 97  112           1  36.7
#> 98  113           1  48.0
#> 99  115           1  50.5
#> 100 116           1  45.5
#> 101 117           1  36.1
#> 102 118           1  38.2
#> 103 119           1  35.0
#> 104 120           1  44.1
#> 105   1           1  51.3
#> 106   1           5  45.5
#> 107   2           1  41.6
#> 108   2           5  47.0
#> 109   3           1  48.2
#> 110   3           5  29.8
#> 111   4           1  57.9
#> 112   5           1  45.7
#> 113   5           5  48.3
#> 114   6           1  50.1
#> 115   6           5  43.9
#> 116   7           1  47.7
#> 117   7           5  41.5
#> 118   8           1  44.9
#> 119   8           5  55.1
#> 120   9           1  58.1
#> 121   9           5  79.1
#> 122  10           1  65.2
#> 123  10           5  75.3
#> 124  11           5  58.4
#> 125  12           1  72.2
#> 126  12           5  64.4
#> 127  13           1  61.5
#> 128  13           5  66.4
#> 129  14           1  70.7
#> 130  14           5  74.2
#> 131  15           1  63.2
#> 132  15           5  52.5
#> 133  16           1  39.8
#> 134  16           5  50.7
#> 135  17           1  53.6
#> 136  17           5  48.9
#> 137  18           1  63.0
#> 138  18           5  63.7
#> 139  19           1  48.4
#> 140  19           5  62.6
#> 141  20           5  71.0
#> 142  21           1  89.3
#> 143  22           5  38.7
#> 144  23           1  35.6
#> 145  23           5  38.6
#> 146  24           1  30.3
#> 147  24           5  37.4
#> 148  25           1  54.7
#> 149  25           5  56.3
#> 150  26           1  49.8
#> 151  26           5  45.7
#> 152  27           1  53.2
#> 153  28           1  39.1
#> 154  29           1  50.8
#> 155  29           5  40.5
#> 156  30           1  56.1
#> 157  30           5  62.2
#> 158  31           1  48.8
#> 159  31           5  48.4
#> 160  32           5  42.2
#> 161  33           1  38.0
#> 162  33           5  40.1
#> 163  34           5  53.5
#> 164  35           1  51.4
#> 165  36           1  44.1
#> 166  36           5  41.7
#> 167  37           5  55.5
#> 168  38           1  60.6
#> 169  38           5  57.4
#> 170  39           1  65.2
#> 171  39           5  62.2
#> 172  40           5  71.7
#> 173  41           1  54.1
#> 174  41           5  47.0
#> 175  42           1  56.9
#> 176  42           5  76.0
#> 177  43           1  42.0
#> 178  43           5  40.3
#> 179  44           5  69.6
#> 180  45           1  44.6
#> 181  45           5  41.5
#> 182  46           1 103.0
#> 183  46           5  95.8
#> 184  47           1  46.4
#> 185  47           5  41.2
#> 186  48           1  72.9
#> 187  48           5  79.2
#> 188  49           1  44.0
#> 189  49           5  58.3
#> 190  50           1  64.8
#> 191  50           5  53.7
#> 192  51           1  61.6
#> 193  51           5  55.1
#> 194  52           1  53.2
#> 195  53           1  34.1
#> 196  53           5  38.8
#> 197  54           1  40.3
#> 198  55           5  46.6
#> 199  56           1  43.9
#> 200  56           5  44.7
#> 201  57           1  60.8
#> 202  57           5  68.3
#> 203  58           5  42.4
#> 204  59           1  40.0
#> 205  59           5  39.6
#> 206  60           5  34.9
#> 207  61           1  79.7
#> 208  61           5  75.3
#> 209  62           1  37.1
#> 210  62           5  43.8
#> 211  63           1  37.6
#> 212  63           5  36.7
#> 213  64           1  52.1
#> 214  64           5  74.6
#> 215  65           1  44.8
#> 216  65           5  55.8
#> 217  66           1  39.9
#> 218  66           5  45.7
#> 219  67           1  39.7
#> 220  67           5  37.6
#> 221  68           1  38.7
#> 222  68           5  34.3
#> 223  69           1  27.6
#> 224  69           5  34.4
#> 225  70           1  63.9
#> 226  70           5  75.9
#> 227  71           1  52.9
#> 228  72           1  41.3
#> 229  72           5  49.0
#> 230  73           1  38.2
#> 231  73           5  35.1
#> 232  74           1  57.1
#> 233  74           5  46.5
#> 234  75           1  52.4
#> 235  75           5  51.6
#> 236  76           1  48.0
#> 237  76           5  52.6
#> 238  77           1  47.2
#> 239  77           5  46.2
#> 240  78           1  48.5
#> 241  79           1  48.1
#> 242  79           5  57.3
#> 243  80           1  64.7
#> 244  80           5  53.9
#> 245  81           1  34.4
#> 246  82           1  59.5
#> 247  83           1  89.3
#> 248  84           1  92.4
#> 249  84           5  53.9
#> 250  85           1  70.9
#> 251  85           5  70.9
#> 252  86           5  38.6
#> 253  87           1  58.4
#> 254  87           5  54.7
#> 255  88           1  45.9
#> 256  88           5  36.0
#> 257  89           1  62.2
#> 258  89           5  68.5
#> 259  90           1  45.8
#> 260  90           5  58.0
#> 261  91           1  25.0
#> 262  91           5  36.9
#> 263  92           1  54.2
#> 264  92           5  70.1
#> 265  93           1  35.5
#> 266  93           5  37.1
#> 267  94           1  48.6
#> 268  95           1  45.7
#> 269  95           5  38.9
#> 270  96           5  61.9
#> 271  97           1  63.8
#> 272  97           5  47.7
#> 273  98           5  57.7
#> 274  99           1  66.7
#> 275 100           1  44.5
#> 276 100           5  45.1
#> 277 101           1  51.2
#> 278 101           5  67.4
#> 279 102           1  62.9
#> 280 102           5  48.9
#> 281 103           1  56.6
#> 282 103           5  46.9
#> 283 104           1  60.6
#> 284 105           1  74.3
#> 285 106           5  66.1
#> 286 107           1  58.7
#> 287 107           5  63.2
#> 288 108           1  63.6
#> 289 108           5  80.1
#> 290 109           1  75.2
#> 291 109           5  56.0
#> 292 110           1  73.6
#> 293 110           5  65.5
#> 294 111           1  48.6
#> 295 111           5  51.2
#> 296 112           1  47.1
#> 297 112           5  50.2
#> 298 113           1  61.5
#> 299 113           5  69.4
#> 300 114           5  67.4
#> 301 115           1  64.8
#> 302 115           5  67.5
#> 303 116           1  58.4
#> 304 116           5  54.4
#> 305 117           1  46.3
#> 306 117           5  46.8
#> 307 118           1  49.1
#> 308 118           5  53.5
#> 309 119           1  50.1
#> 310 119           5  53.9
#> 311 120           1  56.5
#> 312 120           5  52.7
#> 
#> $approx.vc.out$approx.vc.dat.sd
#> # A tibble: 104 × 44
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1  216. 1500.  7.31     0     0     0    60000  65.1  4.07     0     0
#>  2     2  216.  309.  5.73     0     0     0    10000  58.6  2.94     0     0
#>  3     3  216. 1595.  7.37     0     0     0    60000  66.6  3.82     0     0
#>  4     4  216.  664.  6.50     0     0     0    30000  85.4  3.39     0     0
#>  5     5  218. 3369.  8.12     0     0     0   120000  55.8  3.45     0     0
#>  6     6  216.  769.  6.64     0     0     0    30000  52.2  4.26     0     0
#>  7     7  217   269.  5.59     0     0     0    10000  66.5  4.16     0     0
#>  8     8  216. 3428.  8.14     0     0     0   120000  95.1  2.93     0     0
#>  9     9  217. 2651.  7.88     0     0     0   120000  56.0  6.21     0     0
#> 10    10  217   197.  5.28     0     0     0    10000  77.6  7.08     0     0
#> # ℹ 94 more rows
#> # ℹ 32 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>, SteadyState <lgl>,
#> #   recent_ii <dbl>, doses_req_fixed <dbl>, last_two_doses_interval <dbl>,
#> #   doses_req_hl <dbl>, doses_required <dbl>, dose_count_before_obs <int>, …
#> 
#> $approx.vc.out$approx.vc.dat.md
#> # A tibble: 208 × 57
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1  216. 1500.  7.31     0     0     0    60000  65.1  4.07     0     0
#>  2     1  145  1690.  7.43     0     0     0    60000  65.1  4.07    99     0
#>  3     2  216.  309.  5.73     0     0     0    10000  58.6  2.94     0     0
#>  4     2  146.  273.  5.61     0     0     0    10000  58.6  2.94    99     0
#>  5     3  216. 1595.  7.37     0     0     0    60000  66.6  3.82     0     0
#>  6     3  145. 2586.  7.86     0     0     0    60000  66.6  3.82    99     0
#>  7     4  216.  664.  6.50     0     0     0    30000  85.4  3.39     0     0
#>  8     5  218. 3369.  8.12     0     0     0   120000  55.8  3.45     0     0
#>  9     5  144. 3187.  8.07     0     0     0   120000  55.8  3.45    99     0
#> 10     6  216.  769.  6.64     0     0     0    30000  52.2  4.26     0     0
#> # ℹ 198 more rows
#> # ℹ 45 more variables: SD <int>, CMT <int>, resetflag <int>, raw_EVID <dbl>,
#> #   RATE <dbl>, SSflag <int>, route <chr>, dose_number <int>, tad <dbl>,
#> #   dose <dbl>, iiobs <dbl>, rateobs <dbl>, routeobs <chr>, durationobs <dbl>,
#> #   DVstd <dbl>, indiv_lambdaz_eligible <int>, SteadyState <lgl>,
#> #   recent_ii <dbl>, doses_req_fixed <dbl>, last_two_doses_interval <dbl>,
#> #   doses_req_hl <dbl>, doses_required <dbl>, dose_count_before_obs <int>, …
#> 
#> 
#> attr(,"class")
#> [1] "single.point.lst"
```
