#' objective function for size at age
#'
#' @param pars parameter list
#'
#' @export
f_saa <- function(pars) {

  # compiler::enableJIT(0)

  RTMB::getAll(pars, dat_saa)
  linf = exp(log_linf)
  k = exp(log_k)
  alpha = exp(log_alpha)
  beta = exp(log_beta)

  pred = linf * (1 - exp(-k * (age - t0)))
  yvar = log(1. + sd^2 / lbar^2)
  yconst = log(2.0 * pi * yvar * lbar^2)
  rss = 0.5 * (yconst + (log(pred) - log(lbar))^2 / yvar)

  lpred = alpha * log(age) + beta
  lnll = sqrt(n) * (log(lpred) - log(sd))^2
  nll = sum(rss) + sum(lnll)

  RTMB::REPORT(linf)
  RTMB::REPORT(k)
  RTMB::REPORT(t0)
  RTMB::REPORT(alpha)
  RTMB::REPORT(beta)
  return(nll)
}
