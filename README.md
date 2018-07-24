[![CRAN](https://www.r-pkg.org/badges/version/BootValidation)](https://cran.r-project.org/web/packages/BootValidation/index.html)

# BootValidation
This package provides internal bootstrap validation for lineal, logistic, cox and multinomial 'glmnet' models, as well as lm and glm (binomial) regression.

## Usage
```
vboot(glmnet_fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1))
```
## Arguments
 * `glmnet_fit`: Object from glmnet fit 
 * `x`: A matrix of the predictors, each row is an observation vector.
 * `y`  A vector of response variable.
 * `s`  Value of the penalty parameter "lambda" selected from the original cv.glmnet
 * `nfolds` Number of folds for cross validation as in cv.glmnet
 * `B` Number of bootsrap samples
 * `cv_replicates` Number of replicates for the cross-validation step
 * `n_cores` number of cores to use in parallel. Default detectCores()-1

## Details

Main objective of a predictive model is to provide accurated predictions of a new observations. Unfortunately we don´t know how well the model performs. In addition, at the current era of omic data where *p >> n* is not reasonable applying internal validation using data-splitting. Under this background a good method to assessing model performance is applying internal bootstrap validation.                                                                                             
The followed approach is described in Harrel et al. (1996) and on the fantastic [blog](http://thestatsgeek.com/2014/10/04/adjusting-for-optimismoverfitting-in-measures-of-predictive-ability-using-bootstrapping/) written by Jonathan Bartlett. The bootstrap validation procedure consists of the following steps.

   1. Fit the model to original data, and estimate the measure of predictive accuracy *A* (for example AUC from the ROC curve in case of binary outcome or R^2 for numeric outcome). Denote this as *A{orig}*
   2. Repeat this process almost B = 100 or 200 times
      *   Make a bootstrap sample from the original data
      *   Fit the model to the bootstrap sample, and estimate *A* using the fitted model on the bootstrap sample. Denote this as *A_b*
      *   Estimate *A* by applying the fitted bootstrap model on the original dataset. Denote this as *A{b,orig}*
      
   3. Calculate the estimate of optimism *O = B^{-1} \sum_{b=1}^B A_b - A_{b,orig}*
   4. Calculate the optimism adjusted as *A_{orig} - O* 

## Examples

The following example applies internal bootstrap validation on `glmnet` logistic regression. 

 On one hand, we create `data.frame` (`x`) storing the predictors, and `y` the binary outcome. `model.matrix` is required to `glmnet` function. 
```{r}
 x <- data.frame(matrix(rnorm(200), ncol = 200, nrow = 20), factor = factor(rep(c("A", "B"), each = 10)))
 y <- rep(c("Yes", "No"), each = 10)
 x <- model.matrix(~., data = x)
 ```
 On second hand, we fit a penalized elastic net logistic model with alpha = 0.5. 
 
 ```{r}
 library(glmnet)
 cv.lognet <- cv.glmnet(x, y, alpha = 0.5, nfolds = 5, family = "binomial")
 l <- cv.lognet$lambda.1se
 fit_lognet <- glmnet(x, y, alpha = 0.5, family = "binomial")
 ```
 Finally apply `vboot` function to get an original AUC, Optimism and adjusted AUC. 
 
 ```{r}
 vboot(fit_lognet, x, y, nfolds = 3, B = 100, s=l, cv_replicates = 20, n_cores = 3)
```

# Further lines

`BootValidation` package is still on development. The next aim is to set new features to vboot function like internal validation for non penalized regressions.
