# N-Of-1 (Single Subject Design)

This repo was created primarily to serve as an analysis tool for the PCORI project. For public use, please see nof1.pdf for more details. You can run Bayesian linear regression, ordinal/logistic regression, and poisson regression using this package.

# To Install and Load the package

```{r}
library(devtools)
install_github("jiabei-yang/nof1", force = TRUE)
library(nof1)
```

# PCORI Code

```{r}
library(jsonlite)

##### PRODUCE
json.file <- fromJSON("sample input.json")
#json.file <- fromJSON("sample input-small.json")
#json.file <- fromJSON("sample input-no mscd.json")
#json.file <- fromJSON("sample input-no scd.json")
result <- do.call(wrap, json.file)
output <- toJSON(result, pretty = TRUE, UTC = TRUE, auto_unbox = TRUE, na = NULL)
output

######## afib study
json.file2 <- fromJSON("afib sample input.json")
result2 <- do.call(wrap2, json.file2)
output2 <- toJSON(result2, pretty = TRUE, UTC = TRUE, auto_unbox = TRUE, na = NULL)
output2
```


Thank you,    
Michael Seo, Jiabei Yang
