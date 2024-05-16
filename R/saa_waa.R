#' calculate size and weight at age
#'
#' @param age_data goa trawl survey specimen data
#' @param length_data goa trawl survey length data
#' @param age_error age error matrix (since rockfish have a rectangle rather than square)
#' @param len_bins length bins used in assessment model
#' @param rec_age recruitment age used in assessment
#'

#' @export
#'
#' @examples
#' \dontrun{
#' data("age_data")
#' data("length_data")
#' data("age_error")
#' saa_waa(age_data, length_data, age_error, len_bins=15:45, rec_age=2)
#' }
#'
saa_waa <- function(age_data, length_data, age_error, len_bins, rec_age) {
  data <- prep_alw_data(age_data, length_data, age_error, len_bins, rec_age)

  # run size at age
  dat_saa <- data$dat_saa
  par <- list(log_linf = log(max(dat_saa$lbar)),
              log_k = log(0.2),
              t0 = -0.1,
              log_alpha = log(0.07),
              log_beta = log(2.07))

  f_saa <- function(pars) {
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

  obj <- RTMB::MakeADFun(f_saa, par)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
  sd <- RTMB::sdreport(obj)
  report <- obj$report(obj$env$last.par.best)

  # run length-weight
  dat_lw <- data$dat_lw
  par <- list(log_alpha = log(0.07),
              log_beta = log(2.07))
  f_lw <- function(pars) {
    RTMB::getAll(pars, dat_lw)
    alpha = exp(log_alpha)
    beta = exp(log_beta)

    pred = alpha * length^beta
    yvar = log(1. + sd^2 / wbar^2)
    yconst = log(2.0 * pi * yvar * wbar^2)
    rss = sum(0.5 * (yconst + (log(pred) - log(wbar))^2 / yvar))

    wpred = alpha * log(length) + beta

    RTMB::REPORT(alpha)
    RTMB::REPORT(beta)
    return(rss)
  }

  obj <- RTMB::MakeADFun(f_lw, par)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
  sd1 <- RTMB::sdreport(obj)
  report1 <- obj$report(obj$env$last.par.best)
  alpha_lw <- report1$alpha
  beta_lw <- report1$beta

  # run weight at age
  dat_waa <- data$dat_waa
  par = list(log_winf = log(800),
             log_k = log(0.1),
             t0 = 0,
             log_beta = log(beta_lw))

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

  obj <- RTMB::MakeADFun(f_waa, par)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
  sd2 <- RTMB::sdreport(obj)
  report2 <- obj$report(obj$env$last.par.best)


  # Compute Sz@A transition matrix

  Linf = report$linf
  k = report$k
  t0 = report$t0
  alpha = report$alpha
  beta = report$beta
  Lbar = dat_saa$lbar
  sd_Lbar = dat_saa$sd
  ages = dat_saa$age
  Winf = report2$winf
  wk = report2$k
  wt0 = report2$t0
  wbeta = report2$beta



  expand.grid(age = ages,
              length = len_bins) %>%
    dplyr::mutate(Lbar = Linf * (1 - exp(-k * (age - t0))),
                  Lbar = ifelse(age == max(age), 0.5 * (Lbar + Linf), Lbar),
                  sd_Lbar = alpha * log(age) + beta,
                  prob = ifelse(length == min(length),
                                pnorm(length + 0.5, Lbar, sd_Lbar),
                                pnorm(length + 0.5, Lbar, sd_Lbar) -
                                  pnorm(length -0.5, Lbar, sd_Lbar)),
                  prob = round(prob, digits = 4)) %>%
    dplyr::select(age, length, prob) %>%
    tidyr::pivot_wider(names_from = length, values_from = prob) %>%
    # dplyr::mutate(!!rev(names(.))[1] := 1 - rowSums(.[2:(ncol(.) - 1)])) %>%
    dplyr::mutate_at(2:ncol(.), round, 4) -> saa

  data.frame(age = ages) %>%
    dplyr::mutate(wbar = Winf * (1 - exp(-wk * (age - wt0)))^wbeta,
                  wbar = ifelse(age==max(age), 0.5 * (wbar + Winf), wbar),
                  wbar = round(wbar, 1)) -> waa


  list(laa_stats = data$laa_stats,
       lbar_params = data.frame(Linf = Linf,
                                k = k,
                                t0 = t0,
                                a = alpha,
                                b = beta),
       saa = saa,
       wal_stats = data$lw_data,
       alpha_beta_lw = data.frame(alpha = alpha_lw,
                                 beta = beta_lw),
       waa_stats = data$waa_stats,
       wbar_params = data.frame(Winf = Winf,
                                k = wk,
                                t0 = wt0,
                                beta = wbeta),
       waa = waa
       )
}


