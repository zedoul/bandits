context("binom_bandits")

test_that("binom_bandits", {
  # number of success
  x <- c(5,6,7)

  # number of trials
  n <- c(10,20,15)

  .alpha <- 23
  .beta <- 30

  # calc
  binomial_bandit_quad(x, n, .alpha, .beta)
  binomial_bandit_sim(x, n, .alpha, .beta)
})
