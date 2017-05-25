#' Estimate the Bayesian posterior probability of each alternative
#' being the best binomial bandit at time t (quadrature version)
#'
#' @param x number of success
#' @param n number of trials
#' @param .alpha alpha of beta distribution
#' @param .beta beta of beta distribution
#' @param ncp prior center
#' @return optimal probabilities of each arm w_{a, t}
#' @seealso Eq (11)
#' @export
binomial_bandit_quad <- function(x,
                                 n,
                                 .alpha = 1,
                                 .beta = 1,
                                 ncp = 0) {
  # calculate number of bandits
  k <- length(x)
  prob <- rep(0, k)

  for (i in 1:k) {
    # integration function
    # z \in (0, 1)
    # dbeta (density)
    # pbeta (pdf from -\inf to z) where z is z-score
    # Eq (11) in the paper
    f <- function(z) {
      r <- dbeta(z,
                 x[i] + .alpha,
                 n[i] - x[i] + .beta,
                 ncp)

      indices <- c(1:k)[-i]
      for (j in indices) {
        r <- r * pbeta(z,
                       x[j] + .alpha,
                       n[j] - x[j] + .beta,
                       ncp)
      }
      r
    }

    # Calculate posterior based on x, n, and ncp
    prob[i] = integrate(f, 0, 1)$value
  }

  prob
}

#' Estimate the Bayesian posterior probability of each alternative
#' being the best binomial bandit at time t (simulation version)
#'
#' @param x number of success
#' @param n number of trials
#' @param .alpha alpha of beta distribution
#' @param .beta beta of beta distribution
#' @param ncp prior center
#' @return optimal probabilities of each arm w_{a, t}
#' @seealso Eq (11)
#' @export
binomial_bandit_sim <- function(x,
                                n,
                                .alpha = 1,
                                .beta = 1,
                                ndraws = 5000,
                                ncp = 0) {
  # calculate number of bandits
  k <- length(x)
  prob <- matrix(nrow = ndraws,
                 ncol = k)
  no = n - x
  for (i in 1:k) {
    # Calculate posterior based on x, n, and ncp
    # rbeta: random numbers from beta distribution
    prob[, i] = rbeta(ndraws,
                      x[i] + .alpha,
                      no[i] + .beta,
                      ncp)
  }

  # TODO: Calculate credible interval
  w = table(factor(max.col(prob),
            levels = 1:k))
  w/sum(w)
}
