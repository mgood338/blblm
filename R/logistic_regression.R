#' @import purrr
#' @import furrr
#' @import stats
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Logistic Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Computes Logistic Regression for Bag of Little Bootstraps
#'
#'
#' Give a formula, data, value for m, and value for B. The user should
#' run plan(multisession, workers = 4) if they want to use parallelization
#' for example, parallel = TRUE, and they can change the number of workers to
#' any numeric value they wish to use.
#'
#' @param formula a formula
#'
#' @param data dataframe
#' @param m an integer giving  the number of subsets for the data
#' @param B an integer giving the number of bootstraps
#' @param parallel logical value indicating TRUE or FALSE for parallelization
#'
#' @export
#' @return formula
blblog <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (parallel == TRUE) {
    mpa_func <- future_map
  } else {
    mpa_func <- map
  }
  estimates <- mpa_func(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
  )
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblog"
  invisible(res)
}


#' Split data into m parts of approximated equal sizes
#'
#' @param data dataframe
#' @param m integer
#' @return split dataframe
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' Compute the estimates
#'
#' @param formula logistic regression formula
#' @param data dataframe
#' @param n numeric
#' @param B numeric
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}


#' Compute the regression estimates for a blb dataset
#'
#' @param formula a formula
#' @param data dataframe
#' @param n integer
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, n)
}


#' Fitting Linear Model using Bag of Little Bootstraps
#'
#' Estimate the regression estimates based on given the number of repetitions
#'
#' @param formula logistic regression formula
#' @param data dataframe
#' @param freqs numeric value
#' @param n numeric value
#'
#' @return object of class glm
glm1 <- function(formula, data, freqs, n) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- suppressWarnings(glm(formula, family = binomial, data, weights = freqs, maxit = 100))
  while(!fit$converged){
    freqs <- rmultinom(1, n, rep(1, nrow(data)))
    fit <- suppressWarnings(glm(formula, family= binomial, data, weights = freqs, maxit=100))
  }
  list(coef = blblogcoef(fit), sigma = blblogsigma(fit))
}


#' Compute the coefficients from fit
#'
#' @param fit fitted blblm model
#' @return numeric value
blblogcoef <- function(fit) {
  coef(fit)
}


#' Compute sigma from fit
#'
#' @param fit fitted blblm model
#' @return numeric value
blblogsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' Prints the Logistic Model
#'
#' @param x an object
#'
#' @param ... additional arguments affecting the result
#'
#' @export
#' @return formula of blblog
#' @method print blblog
print.blblog <- function(x, ...) {
  cat("blblog model:", capture.output(x$formula))
  cat("\n")
}

#' Computing Sigma
#'
#' Computes the standard deviation/error of Blblm using Polymorphism
#'
#' @param object a fitted model object
#'
#' @param confidence logical value indicating whether we construct the confidence interval
#' @param level integer between 0 and 1 indicating the confidence level
#' @param ... additional arguments affecting the result
#'
#' @export
#' @return numeric
#' @method sigma blblog
sigma.blblog <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Computing the Coefficients for the Model
#'
#' @param object a fitted model object
#'
#' @param ... additional arguments affecting the result
#'
#' @export
#' @method coef blblog
coef.blblog <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Confidence Interval
#'
#' @param object a fitted model object
#'
#' @param parm numeric value indicating the parameters you want to construct confidence interval for
#' @param level integer between 0 and 1 indicating the confidence level
#' @param ... additional arguments affecting the result
#'
#' @export
#' @return numeric confidence interval
#' @method confint blblog
confint.blblog <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

invlogit <- function(x) {
  1 / (1 + exp(-x))
}

#' Prediction Function
#'
#' Computes the predictions and the confidence interval for the predictions
#'
#' @param object a fitted model object
#'
#' @param new_data data frame
#' @param confidence logical value indicating whether we construct the confidence interval
#' @param level integer between 0 and 1 indicating the confidence level
#' @param ... additional arguments affecting the predictions produced
#'
#' @export
#' @return A matrix, vector, or list with the fit, upper and lower value
#' @method predict blblog
predict.blblog <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
