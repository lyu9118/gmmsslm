#' Bayes' rule of allocation
#'
#' Bayes' rule of allocation
#'
#' Classifier specified by Bayes' rule
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param params A list containing the required parameters (p, g, pi, mu, and sigma) from sslparlist() or provided manually.
#' @return
#' \item{clust}{Class membership for the ith entity}
#' @details
#' The classifier/Bayes rule of allocation \eqn{R(y_j;\theta)} assigns an entity with observation \eqn{y_j} to  class \eqn{C_k }(that is, \eqn{R(y_j;\theta)=k}) if
#' \eqn{	k=\arg\max_i \tau_i(y_j;\theta),}
#' @export
#' @import mvtnorm
#'
#' @examples
#'n <- 150
#'pi <- c(0.25, 0.25, 0.25, 0.25)
#'sigma <- array(0, dim = c(3, 3, 4))
#'sigma[, , 1] <- diag(1, 3)
#'sigma[, , 2] <- diag(2, 3)
#'sigma[, , 3] <- diag(3, 3)
#'sigma[, , 4] <- diag(4, 3)
#'mu <- matrix(c(0.2, 0.3, 0.4, 0.2, 0.7, 0.6, 0.1, 0.7, 1.6, 0.2, 1.7, 0.6), 3, 4)
#'dat <- rmix(n = n, pi = pi, mu = mu, sigma = sigma)
#'params <- list(p = 3, g = 4, pi = pi, mu = mu, sigma = sigma)
#'clust <- bayesclassifier(dat = dat$Y, params = params)

bayesclassifier <- function(dat, params){
  n <- dim(dat)[1]
  p <- params$p
  g <- params$g
  pi <- params$pi
  mu <- params$mu
  sigma <- params$sigma

  logdens <- matrix(0, n, g)
  ncov <- ifelse(is.na(dim(sigma)[3]), 1, 2)

  if (ncov == 1) {
    for (i in 1:g) {
      logdens[, i] <- mvtnorm::dmvnorm(dat, mean = mu[, i], sigma = as.matrix(sigma), log = TRUE)
    }
  } else {
    for (i in 1:g) {
      logdens[, i] <- mvtnorm::dmvnorm(dat, mean = mu[, i], sigma = as.matrix(sigma[, , i]), log = TRUE)
    }
  }

  logprobs <- t(t(logdens) + log(pi))
  clusprobs <- t(apply(logprobs, 1, normalise_logprob))
  clust <- max.col(clusprobs)
  return(clust)
}
