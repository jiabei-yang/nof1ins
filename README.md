# nof1gen (N-of-1 study general analysis tool)

This repo was created primarily to serve as an analysis tool for n-of-1 studies. Currently, the code is tested on the following outcomes and models, and the remaining models are under development or being tested.

| Models |  Continuous  |     Count    |    Binary    |   Ordinal    |   Nominal    |
| ------ | :----------: | :----------: | :----------: | :----------: | :----------: | 
| Trend and correlation off (Mean-only model) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Trend on, correlation off | :heavy_check_mark: |
| Trend off, correlation on | :heavy_check_mark: |
| Trend and correlation on  | :heavy_check_mark: |

## To install and load the package

```{r}
library(devtools)
install_github("jiabei-yang/nof1gen", force = TRUE)
library(nof1gen)
```

## Prepare data
```{r}
data.used <- data.frame(Y     = c(9, 9, 10, 11, 9, 7, 6, 
                                  6, 6, 6, 5, 5, 6, 4, 
                                  6, 6, 7, 7, 9, 8, 8, 
                                  7, NA, 6, 6, 7, NA, 8),
				Treat = rep(c(1, 2, 1, 2), each = 7),
				Day   = 1:28,
				Block = rep(c(1, 2), each = 14))
```


## Run models
### Trend and correlation off (Mean-only model)
```{r}
nof1_part1_bsF_corrF <- nof1.data(Y              = observations$Y, 
                                  Treat          = observations$Treat, 
                                  response       = "normal", 
                                  ncat           = NULL, 
                                  bs.trend       = F,
                                  y.time         = NULL, 
                                  knots.bt.block = NULL,
                                  block.no       = NULL,
                                  corr.y         = F,
                                  alpha.prior    = NULL, 
                                  beta.prior     = NULL, 
                                  eta.prior      = NULL, 
                                  dc.prior       = NULL, 
                                  c1.prior       = NULL,
                                  rho.prior      = NULL, 
                                  hy.prior       = NULL)
```


Thank you,    
Jiabei Yang
