# Calculate metrics for model predictive performance evaluation

Computes common error metrics that quantify the predictive performance
of pharmacometric models by comparing predicted (pred.x) and observed
(obs.y) concentration values.

## Usage

``` r
metrics.(pred.x, obs.y)
```

## Arguments

- pred.x:

  Numeric vector of model-predicted values.

- obs.y:

  Numeric vector of corresponding observed values.

## Value

A numeric vector with named elements:

- APE: absolute prediction error

- MAE: mean absolute error

- MAPE: mean absolute percentage error

- RMSE: root mean squared error

- rRMSE1: relative RMSE (type 1)

- rRMSE2: relative RMSE (type 2)

## Details

The function stops with an error if pred.x and obs.y have unequal
lengths. The following metrics are calculated:

\$\$APE = \sum \|pred.x - obs.y\|\$\$ Absolute prediction error (APE) is
the sum of absolute differences.

\$\$MAE = \frac{1}{n} \sum \|pred.x - obs.y\|\$\$ Mean absolute error
(MAE) expresses the average absolute deviation.

\$\$MAPE = \frac{100}{n} \sum \left\| \frac{pred.x - obs.y}{obs.y}
\right\|\$\$ Mean absolute percentage error (MAPE) normalizes the error
by observed values.

\$\$RMSE = \sqrt{\frac{1}{n} \sum (pred.x - obs.y)^2}\$\$ Root mean
squared error (RMSE) penalizes larger deviations.

\$\$rRMSE1 = \frac{RMSE}{\bar{obs.y}} \times 100\$\$ Relative RMSE type
1 is the RMSE normalized by the mean observed value.

\$\$rRMSE2 = 100 \times \sqrt{\frac{1}{n} \sum \left( \frac{pred.x -
obs.y}{(pred.x + obs.y)/2} \right)^2}\$\$ Relative RMSE type 2 is
symmetric and normalizes by the mean of each predicted–observed pair.

## Examples

``` r

obs.y  <- rnorm(100, mean = 100, sd = 10)
pred.x <- obs.y + rnorm(100, mean = 0, sd = 5)
metrics.(pred.x = pred.x, obs.y = obs.y)
#>    metrics.ape    metrics.mae   metrics.mape   metrics.rmse metrics.rrmse1 
#>     391.314477       3.913145       3.962923       4.779439       4.791252 
#> metrics.rrmse2 
#>       4.853167 
```
