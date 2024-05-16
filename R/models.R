f_lw <- function(pars) {
  RTMB::getAll(pars, dat_lw)
  alpha = exp(log_alpha)
  beta = exp(log_beta)
  n = 1:length(ages)

  pred = alpha * length^beta
  yvar = log(1. + sd^2 / wbar^2)
  yconst = log(2.0 * pi * yvar * wbar^2)
  rss = sum(0.5 * (yconst + (log(pred) - log(wbar))^2 / yvar))

  wpred = alpha * log(length) + beta

  RTMB::REPORT(alpha)
  RTMB::REPORT(beta)
  return(rss)
}

f_waa <- function(pars) {
  RTMB::getAll(pars, dat_waa)
  winf = exp(log_winf)
  k = exp(log_k)
  beta = exp(log_beta)

  pred = winf * (1 - exp(-1.0 * k * (age - t0)))^beta;
  yvar = log(1. + sd^2 / wbar^2)
  yconst = log(2.0 * pi * yvar * wbar^2)
  rss = sum(0.5 * (yconst + (log(pred) - log(wbar))^2 / yvar))

  RTMB::REPORT(winf)
  RTMB::REPORT(k)
  RTMB::REPORT(t0)
  RTMB::REPORT(beta)
  return(rss)
}

f_saa <- function(pars, dat_saa) {
  # RTMB::getAll(pars, dat_saa)
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
