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
nof1_bsF_corrF <- nof1.data(Y              = data.used$Y, 
                            Treat          = data.used$Treat, 
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

result_bsF_corrF <- nof1.run(nof1_bsF_corrF, 
                             inits           = NULL, 
                             n.chains        = 3,
                             max.run         = 100000, 
                             setsize         = 10000, 
                             n.run           = 50000,
                             conv.limit      = 1.05, 
                             extra.pars.save = NULL)
```

### Trend on, correlation off
```{r}
nof1_bsT_corrF <- nof1.data(Y              = data.used$Y, 
                            Treat          = data.used$Treat, 
                            response       = "normal", 
                            ncat           = NULL, 
                            bs.trend       = T,
                            y.time         = data.used$Day, 
                            knots.bt.block = T,
                            block.no       = data.used$Block,
                            corr.y         = F,
                            alpha.prior    = NULL, 
                            beta.prior     = NULL, 
                            eta.prior      = NULL, 
                            dc.prior       = NULL, 
                            c1.prior       = NULL,
                            rho.prior      = NULL, 
                            hy.prior       = NULL)

result_bsT_corrF <- nof1.run(nof1_bsT_corrF, 
                             inits           = NULL, 
                             n.chains        = 3,
                             max.run         = 100000, 
                             setsize         = 10000, 
                             n.run           = 50000,
                             conv.limit      = 1.05, 
                             extra.pars.save = NULL)
```

### Trend off, correlation on 
```{r}
nof1_bsF_corrT <- nof1.data(Y              = data.used$Y, 
                            Treat          = data.used$Treat, 
                            response       = "normal", 
                            ncat           = NULL, 
                            bs.trend       = F,
                            y.time         = NULL, 
                            knots.bt.block = NULL,
                            block.no       = NULL,
                            corr.y         = T,
                            alpha.prior    = NULL, 
                            beta.prior     = NULL, 
                            eta.prior      = NULL, 
                            dc.prior       = NULL, 
                            c1.prior       = NULL,
                            rho.prior      = NULL, 
                            hy.prior       = NULL)

result_bsF_corrT <- nof1.run(nof1_bsF_corrT, 
                             inits           = NULL, 
                             n.chains        = 3,
                             max.run         = 100000, 
                             setsize         = 10000, 
                             n.run           = 50000,
                             conv.limit      = 1.05, 
                             extra.pars.save = NULL)
```
### Trend on, correlation on
```{r}
nof1_bsT_corrT <- nof1.data(Y              = data.used$Y, 
                            Treat          = data.used$Treat, 
                            response       = "normal", 
                            ncat           = NULL, 
                            bs.trend       = T,
                            y.time         = data.used$Day, 
                            knots.bt.block = T,
                            block.no       = data.used$Block,
                            corr.y         = T,
                            alpha.prior    = NULL, 
                            beta.prior     = NULL, 
                            eta.prior      = NULL, 
                            dc.prior       = NULL, 
                            c1.prior       = NULL,
                            rho.prior      = NULL, 
                            hy.prior       = NULL)

result_bsT_corrT <- nof1.run(nof1_bsT_corrT, 
                             inits           = NULL, 
                             n.chains        = 3,
                             max.run         = 100000, 
                             setsize         = 10000, 
                             n.run           = 50000,
                             conv.limit      = 1.05, 
                             extra.pars.save = NULL)
```




Thank you,    
Jiabei Yang
