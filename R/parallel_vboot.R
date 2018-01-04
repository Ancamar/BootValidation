#' Internal bootstraping validation logistic glmnet model
#'
#' @description Validate glmnet logistic regression using bootstrap.
#' @param glmnet_fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be a factor with two levels
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param nfolds Number of folds for cross validation as in cv.glmnet
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step in 'cv.glmnet'
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.lognet <- function(glmnet_fit, x, y, s, nfolds, B, cv_replicates, n_cores){
    orig_predict <- glmnet::predict.glmnet(glmnet_fit, newx = x, s = s, type = "response")
    orig_auc <- pROC::roc(y, as.numeric(orig_predict))$auc

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha,
                                                                           family = "binomial", nfolds = nfolds)$lambda.1se, error = function(e) NA))
        l = median(cv.glmnet_b, na.rm = TRUE)
        boot_fit <- glmnet::glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha, family = "binomial")
        boot_predict <- glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response")
        Cb_boot <- pROC::roc(yboot, as.numeric(boot_predict), direction = "<")$auc

        # fit bootstrap model to the original dataset
        bootorig_predict <-  glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response")
        Cb_orig <- pROC::roc(y, as.numeric(bootorig_predict), direction = "<")$auc
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig))
    }
    if(n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "glmnet_fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds)))
    }

    # Optimist
    O <- B^-1 * sum(unlist(CBOOT[1, ]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- orig_auc - O
    output <- round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj), 3)
    return(output)
}


#' Internal bootstraping validation lineal glmnet model
#'
#' @description Validate glmnet linear regression using bootstrap.
#' @param glmnet_fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be numeric
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.elnet <- function(glmnet_fit, x, y, s, nfolds, B, cv_replicates, n_cores){
    orig_predict <- glmnet::predict.glmnet(glmnet_fit, newx = x, s=s, type = "response")
    orig_r2 <- 1 - (var(orig_predict - y) / var(y))
    # Create index to bootstrap
    bootstrap <- function(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- pbapply::pbreplicate(cv_replicates, glmnet::cv.glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha, family = "gaussian", nfolds = nfolds)$lambda.1se)
        l <- median(cv.glmnet_b)
        boot_fit <- glmnet::glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha, family = "gaussian")
        boot_predict <- glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response")
        Cb_boot <- 1 - (var(boot_predict - yboot) / var(yboot))

        # fit bootstrap model to the original dataset
        bootorig_predict <-  glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response")
        Cb_orig <- 1 - (var(bootorig_predict - y) / var(y))
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig))
    }
    if( n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "glmnet_fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds)))
    }
    # Optimist
    O <- B^-1 * sum(unlist(CBOOT[1, ]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- as.numeric(orig_r2 - O)
    output <- round(data.frame(Original_R2 = as.numeric(orig_r2), Optimism = O, Validated_R2 = Oadj),3)
    return(output)
}

#' Internal bootstraping validation cox glmnet model
#'
#' @description Validate glmnet cox regression using bootstrap.
#' @param glmnet_fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y Should be a two-column matrix with columns named 'time' and 'status' as in 'glmnet'
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.coxnet
#' @importFrom survAUC AUC.cd
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.coxnet <- function(glmnet_fit, x, y, s, nfolds, B, cv_replicates, n_cores){
    times <- times <- seq(1,max(as.numeric(y)),1)
    orig_predict <- glmnet::predict.coxnet(glmnet_fit, newx = x, s = s, type = "response")
    orig_auc <- survAUC::AUC.cd(y,y, orig_predict, orig_predict, times)$iauc

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha,
                                                                           family = "cox", nfolds = nfolds)$lambda.1se, error = function(e) NA))
        l = median(cv.glmnet_b, na.rm = TRUE)
        boot_fit <- glmnet::glmnet(xboot, yboot, alpha = glmnet_fit$call$alpha, family = "cox")
        boot_predict <- glmnet::predict.coxnet(boot_fit, newx = xboot, s = l, type = "response")
        Cb_boot <- survAUC::AUC.cd(y, yboot,orig_predict, boot_predict, times)$iauc

        # fit bootstrap model to the original dataset
        bootorig_predict <-  glmnet::predict.coxnet(boot_fit, newx = x, s = l, type = "response")
        Cb_orig <- survAUC::AUC.cd(y, y, bootorig_predict, bootorig_predict, times)$iauc
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig))
    }
    if(n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "glmnet_fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = glmnet_fit$call$alpha, nfolds = nfolds)))
    }

    # Optimist
    O <- 10^-1 * sum(unlist(CBOOT[1, ]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- orig_auc - O
    output <- round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj), 3)
    return(output)
}

#' Generic function for bootstrap validation
#'
#' @description Validate 'glmnet' linear, logistic or cox regression using bootstrap.
#' @param glmnet_fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @references Frank Harrell (2015). Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom survAUC AUC.cd
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
#' @examples
#' # Create the data
#' set.seed(25)
#' x <- matrix(rnorm(80),ncol=4)
#' y <- x[,4]*0.8+x[,3]*0.4+rnorm(20)
#' # Fit glmnet model
#' fit_enet <- glmnet::glmnet(x, y, alpha = 0.5)
#' # Bootstrap validation
#' vboot(fit_enet, x, y, nfolds = 3, B = 2, s = 0.5, cv_replicates = 5, n_cores = 1)

vboot <- function(glmnet_fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    n_cores <- max(1, n_cores)
    UseMethod("vboot")
}

