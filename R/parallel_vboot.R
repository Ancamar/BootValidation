#' Internal bootstraping validation logistic model
#'
#' @description Validate logistic regression using bootstrap.
#' @param fit Object from glm fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom stats predict glm formula
#' @export
vboot.glm <- function(fit, x = NULL, y = NULL, s = NULL, nfolds = NULL, B = 200, cv_replicates = NULL, n_cores = max(1, parallel::detectCores() - 1)){
    orig_predict <- predict(fit,type = "response")
    orig_auc <- pROC::roc(fit$y, as.numeric(orig_predict))$auc

    x <- fit$data[,!names(fit$data) %in% names(fit$model[1])]
    y <- fit$model[1]
    # Making index to bootstrap
    bootstrap <- function(x, y, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index,]
        databoot <- data.frame(yboot, xboot)
        names(databoot)[1] <- names(fit$model)[1]
        # Fit the model using bootstrap dataset
        boot_fit <- glm(formula(fit), data = databoot, family = fit$family)
        boot_predict <- predict.glm(boot_fit, newdata = databoot ,type = "response")
        Cb_boot <- pROC::roc(yboot, as.numeric(boot_predict), direction = "<")$auc

        # fit bootstrap model to the original dataset
        bootorig_predict <-  predict.glm(boot_fit, newdata = fit$data, type = "response")
        Cb_orig <- pROC::roc(fit$y, as.numeric(bootorig_predict), direction = "<")$auc
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig))
    }
    if(n_cores > 1 & grepl("ns\\(" ,deparse(formula(fit), width.cutoff = 500L))){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "bootstrap", "ns"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        if(n_cores > 1){
            cl <- parallel::makeCluster(n_cores)
            parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "bootstrap"), envir = environment())
            CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y))
            parallel::stopCluster(cl)
            closeAllConnections()
        } else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y)))
        }
    }

    # Optimist
    O <- B^-1 * sum(unlist(CBOOT[1, ]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- orig_auc - O
    output <- round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj), 3)
    return(output)
}


#' Internal bootstraping validation linear model
#'
#' @description Validate linear regression using bootstrap.
#' @param fit Object from lm fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom stats var predict lm formula
#' @export
vboot.lm <- function(fit, x = NULL, y = NULL, s = NULL, nfolds = NULL, B = 200, cv_replicates = NULL, n_cores = max(1, parallel::detectCores() - 1)){
    x <- eval(fit$call$data)[, !names(eval(fit$call$data)) %in% names(fit$model)[1]]
    y <- fit$model[,1]
    orig_predict <- predict.lm(fit, type = "response")
    orig_r2 <- 1 - (var(orig_predict - y) / var(y))
    # Create index to bootstrap
    bootstrap <- function(x, y, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]
        databoot <- data.frame(yboot, xboot)
        names(databoot)[1] <- names(fit$model)[1]
        # Fit the model using bootstrap dataset
        boot_fit <- lm(formula(fit), data = databoot)
        boot_predict <- predict(boot_fit, newdata = databoot, type = "response")
        Cb_boot <- 1 - (var(boot_predict - yboot) / var(yboot))

        # fit bootstrap model to the original dataset
        bootorig_predict <-  predict(boot_fit, newdata = x, type = "response")
        Cb_orig <- 1 - (var(bootorig_predict - y) / var(y))
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig))
    }
    if(n_cores > 1 & grepl("ns\\(" ,deparse(formula(fit), width.cutoff = 500L))){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "bootstrap", "ns"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        if(n_cores > 1){
            cl <- parallel::makeCluster(n_cores)
            parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "bootstrap"), envir = environment())
            CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y))
            parallel::stopCluster(cl)
            closeAllConnections()
        } else{
            CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y)))
        }
    }
    # Optimist
    O <- B^-1 * sum(unlist(CBOOT[1,]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- as.numeric(orig_r2 - O)
    #output
    output <- round(data.frame(Original_R2 = as.numeric(orig_r2), Optimism = O, Validated_R2 = Oadj),3)
    return(output)

}



#' Internal bootstraping validation logistic glmnet model
#'
#' @description Validate glmnet logistic regression using bootstrap.
#' @param fit Object from glmnet fit
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
#' @importFrom stats median coef
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.lognet <- function(fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    orig_predict <- glmnet::predict.glmnet(fit, newx = x, s = s, type = "response")
    orig_auc <- pROC::roc(y, as.numeric(orig_predict))$auc

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- suppressWarnings(pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha,
                                                                                                       family = "binomial", nfolds = nfolds)$lambda.1se, error = function(e) NA)))
        l <- median(cv.glmnet_b, na.rm = TRUE)
        boot_fit <- suppressWarnings(tryCatch(glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "binomial"), error = function(e) NA))
        boot_predict <- tryCatch(glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response"), error = function(e) NA)
        Cb_boot <- tryCatch(pROC::roc(yboot, as.numeric(boot_predict), direction = "<")$auc, error = function(e) NA)
        #selected variables
        var_sel <- tryCatch(names(coef(boot_fit, s = l)[coef(boot_fit, s = l)[,1] != 0,])[-1], error = function(e) NA)
        # fit bootstrap model to the original dataset
        bootorig_predict <-  tryCatch(glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response"), error = function(e) NA)
        Cb_orig <- tryCatch(pROC::roc(y, as.numeric(bootorig_predict), direction = "<")$auc, error = function(e) NA)
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if(n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) tryCatch(bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds), error = function(e) NA))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
    }

    # Optimist
    Optimisms <- as.numeric(stats::na.omit(unlist(CBOOT[1, ])) - stats::na.omit(unlist(CBOOT[2, ])))
    eff_B <- length( stats::na.omit(unlist(CBOOT[2, ])))
    O <- eff_B ^-1 * sum(Optimisms)
    # Adjusted Optimist
    if(length( stats::na.omit(unlist(CBOOT[2, ]))) != B) warning(paste("Due to sample size effective bootstraps samples are",eff_B))
    Oadj <- orig_auc - O
    output <- list(round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj),3), varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B, Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)
}


#' Internal bootstraping validation linear glmnet model
#'
#' @description Validate glmnet linear regression using bootstrap.
#' @param fit Object from glmnet fit
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
#' @importFrom stats median var coef
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.elnet <- function(fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    orig_predict <- glmnet::predict.glmnet(fit, newx = x, s=s, type = "response")
    orig_r2 <- 1 - (var(orig_predict - y) / var(y))
    # Create index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- pbapply::pbreplicate(cv_replicates, glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha, family = "gaussian", nfolds = nfolds)$lambda.1se)
        l <- median(cv.glmnet_b)
        boot_fit <- glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "gaussian")
        boot_predict <- glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response")
        Cb_boot <- 1 - (var(boot_predict - yboot) / var(yboot))
        #selected variables
        var_sel <- names(coef(fit, s = s)[coef(fit, s = s)[,1] != 0,])[-1]

        # fit bootstrap model to the original dataset
        bootorig_predict <-  glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response")
        Cb_orig <- 1 - (var(bootorig_predict - y) / var(y))
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if( n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds)))
    }
    # Optimist
    Optimisms <- as.numeric(stats::na.omit(unlist(CBOOT[1, ])) - stats::na.omit(unlist(CBOOT[2, ])))
    eff_B <- length( stats::na.omit(unlist(CBOOT[2, ])))
    O <- B^-1 * sum(unlist(CBOOT[1,]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- as.numeric(orig_r2 - O)
    #output
    output <- list(round(data.frame(Original_R2 = as.numeric(orig_r2), Optimism = O, Validated_R2 = Oadj),3), varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B, Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)

}


#' Internal bootstraping validation cox glmnet model
#'
#' @description Validate glmnet cox regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y Should be a two-column matrix with columns named 'time' and 'status' as in 'glmnet'
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.coxnet
#' @importFrom risksetROC CoxWeights IntegrateAUC
#' @importFrom survival survfit Surv
#' @importFrom pbapply pbreplicate
#' @importFrom stats median coef
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.coxnet <- function(fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    orig_predict <- glmnet::predict.coxnet(fit, newx = x, s = s, type = "link")
    surv.prob <- unique(survival::survfit(survival::Surv(y[,1], y[,2])~1)$surv)
    utimes <- unique( y[,1][ y[,2] == 1 ])
    utimes <- utimes[order(utimes)]
    ## find AUC at unique failure times
    AUC <- rep( NA, length(utimes) )
    for( j in 1:length(utimes) )
    {
        out <- risksetROC::CoxWeights( orig_predict, y[,1], y[,2],utimes[j])
        AUC[j] <- out$AUC
    }
    ## integrated AUC to get concordance measure
    orig_auc <- risksetROC::IntegrateAUC( AUC, utimes, surv.prob, tmax=Inf)

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index, ]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- suppressWarnings(pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha,
                                                                                                       family = "cox", nfolds = nfolds)$lambda.1se, error = function(e) NA)))
        l <- median(cv.glmnet_b, na.rm = TRUE)
        #add trychatch to return NA in case of too many event size in the bootstrap sample
        boot_fit <- suppressWarnings(tryCatch(glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "cox"), error = function(e) NA))
        boot_predict <- tryCatch(glmnet::predict.coxnet(boot_fit, newx = xboot, s = l, type = "link"),error = function(e) NA)

        surv.probb <- unique(survival::survfit(survival::Surv(yboot[,1], yboot[,2])~1)$surv)
        utimes <- unique( yboot[,1][ yboot[,2] == 1 ])
        utimes <- utimes[order(utimes)]
        ## find AUC at unique failure times
        AUCb <- rep( NA, length(utimes) )
        for( j in 1:length(utimes) )
        {
            outb <- tryCatch(risksetROC::CoxWeights( boot_predict, yboot[,1], yboot[,2],utimes[j]), error = function(e) NA)
            if(is.logical(outb)){
                AUCb[j] <- NA
            }else{
                AUCb[j] <- outb$AUC
            }
        }
        ## integrated AUC to get concordance measure
        Cb_boot  <- tryCatch(risksetROC::IntegrateAUC( AUCb, utimes, surv.probb, tmax=Inf), error = function(e) NA)
        #selected variables
        var_sel <- tryCatch(names(coef(fit, s = s)[coef(fit, s = s)[,1] != 0,])[-1],error = function(e) NA)
        # fit bootstrap model to the original dataset
        bootorig_predict <-  tryCatch(glmnet::predict.coxnet(boot_fit, newx = x, s = l, type = "link"),error = function(e) NA)

        surv.probc <- unique(survival::survfit(survival::Surv(yboot[,1], yboot[,2])~1)$surv)
        utimes <- unique( yboot[,1][ yboot[,2] == 1 ])
        utimes <- utimes[order(utimes)]
        ## find AUC at unique failure times
        AUCc <- rep( NA, length(utimes) )
        for( j in 1:length(utimes) )
        {
            outc <- tryCatch(risksetROC::CoxWeights( bootorig_predict, y[,1], y[,2],utimes[j]), error = function(e) NA)
            if(is.logical(outc)){
                AUCc[j] <- NA
            }else{
                AUCc[j] <- outc$AUC
            }
        }
        ## integrated AUC to get concordance measure
        Cb_orig  <- tryCatch(risksetROC::IntegrateAUC( AUCc, utimes, surv.probc, tmax=Inf), error = function(e) NA)

        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if(n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds)))
    }

    # Optimist
    Optimisms <- as.numeric(stats::na.omit(unlist(CBOOT[1, ])) - stats::na.omit(unlist(CBOOT[2, ])))
    eff_B <- length( stats::na.omit(unlist(CBOOT[2, ])))
    O <- eff_B^-1 * sum(Optimisms)
    # Adjusted Optimist
    if(length( stats::na.omit(unlist(CBOOT[2, ]))) != B) warning(paste("Due to sample size effective bootstraps samples are",eff_B))
    Oadj <- orig_auc - O
    output <- list(round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj), 3), varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B, Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)
}



#' Internal bootstraping validation multinomial glmnet model
#'
#' @description Validate glmnet logistic regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be a factor with two levels
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param nfolds Number of folds for cross validation as in cv.glmnet
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step in 'cv.glmnet'
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet predict.multnet
#' @importFrom pROC multiclass.roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median coef
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.multnet <- function(fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    if(sort(table(y))[1] < 3) stop("one multinomial or binomial class has 2 or less observations")
    orig_predict <- glmnet::predict.multnet(fit, newx = x, s = s, type = "class")
    orig_auc <- pROC::multiclass.roc(y, as.numeric(as.factor(orig_predict)))$auc

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- suppressWarnings(pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha,
                                                                                                       family = "multinomial", nfolds = nfolds)$lambda.1se, error = function(e) NA)))
        l = median(cv.glmnet_b, na.rm = TRUE)
        #add trychatch to return NA in case of too many event size in the bootstrap sample
        boot_fit <- suppressWarnings(tryCatch(glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "multinomial"), error = function(e) NA))
        boot_predict <- tryCatch(glmnet::predict.multnet(boot_fit, newx = xboot, s = l, type = "class"), error = function(e) NA)
        Cb_boot <- tryCatch(pROC::multiclass.roc(yboot, as.numeric(as.factor(boot_predict)), direction = "<")$auc, error = function(e) NA)
        #selected variables
        var_sel <- tryCatch(names(coef(fit, s = s)[coef(fit, s = s)[,1] != 0,])[-1], error = function(e) NA)
        # fit bootstrap model to the original dataset
        bootorig_predict <-  tryCatch(glmnet::predict.multnet(boot_fit, newx = x, s = l, type = "class"), error = function(e) NA)
        Cb_orig <- tryCatch(pROC::multiclass.roc(y, as.numeric(as.factor(bootorig_predict)), direction = "<")$auc, error = function(e) NA)
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if(n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds)))
    }

    # Optimist
    Optimisms <- as.numeric(stats::na.omit(unlist(CBOOT[1, ])) - stats::na.omit(unlist(CBOOT[2, ])))
    eff_B <- length( stats::na.omit(unlist(CBOOT[2, ])))
    O <- eff_B^-1 * sum(Optimisms)
    # Adjusted Optimist
    if(length( stats::na.omit(unlist(CBOOT[2, ]))) != B) warning(paste("Due to sample size effective bootstraps samples are",eff_B))
    Oadj <- orig_auc - O
    output <- list(round(data.frame(Original_AUC = as.numeric(orig_auc), Optimism = O, Validated_AUC = Oadj), 3), varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B, Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)
}


#' Print function
#'
#' @description internal function to print vboot object
#' @param x A vboot object
#' @param ... Further arguments of generic function
#' @return validation
#' @export
print.bootVal <- function(x, ...){
    print(x[[1]])
}


#' Plot for repeat_cv results
#'
#' @description Plots a grid of slices from the estimates of the repeat_cv function
#' @param x A vboot object
#' @param order order plot by importance
#' @param n Variables to be displayed
#' @param ... further arguments passed to plot
#' @importFrom graphics barplot
#' @return validation
#' @export
plot.bootVal <- function(x, order = TRUE, n = length(unique(unlist(x$varImportance))), ...){
    raw_table <- table(unlist(x[[2]]))[1:n]/x$Effective_Bootstraps
    if(order) raw_table <- raw_table[order(unlist(x[[2]]), decreasing = TRUE)]
    barplot(raw_table, horiz = TRUE,
            las = 1, xlim = c(0,1), main = "Importance variable", xlab = "Inclusion porportion", ...)
}


#' Generic function for bootstrap validation
#'
#' @description Validate 'glmnet' linear, logistic or cox regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @references Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software, 39(5), 1-13. URL http://www.jstatsoft.org/v39/i05/.
#' @references Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
#' @references Gordon C.S. Smith, Shaun R. Seaman, Angela M. Wood, Patrick Royston, Ian R. White (2014). Correcting for Optimistic Prediction in Small Data Sets, American Journal of Epidemiology, Volume 180, Issue 3, 1 August 2014, Pages 318-324, https://doi.org/10.1093/aje/kwu140
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet predict.multnet
#' @importFrom pROC roc multiclass.roc
#' @importFrom risksetROC CoxWeights IntegrateAUC
#' @importFrom survival survfit Surv
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var predict.glm predict.lm glm lm coef formula
#' @importFrom graphics barplot
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

vboot <- function(fit, x, y, s, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    UseMethod("vboot")
}


