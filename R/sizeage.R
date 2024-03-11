#' size at age matrix
#'
#' @param age_data from ...
#' @param length_data from ...
#' @param len_bins length bins
#'
#' @export
#'
saa <- function(age_data, length_data, len_bins){
  age_data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::select(year, age, length) %>%
    dplyr::filter(year>=1990, !is.na(age))  %>%
    dplyr::select(-year) %>%
    dplyr::group_by(age) %>%
    dplyr::filter(dplyr::n()>1) %>%
    dplyr::group_by(length) %>%
    dplyr::mutate(n_l = dplyr::n()) %>%
    dplyr::arrange(age, length) %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(sample_size =  dplyr::n()) -> inter

  length_data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year>=1990, !is.na(length)) %>%
    dplyr::group_by(length) %>%
    dplyr::summarise(tot = dplyr::n()) -> dat

  dat %>%
    dplyr::left_join(inter, .) %>%
    dplyr::group_by(age, length) %>%
    dplyr::mutate(prop =  dplyr::n() / n_l * tot) %>%
    dplyr::distinct() %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(sample_size = mean(sample_size),
                     Lbar = sum(prop * length) / sum(prop) * 0.1,
                     sd_Lbar = sqrt(1 / (sum(prop) - 1) * sum(prop * (length / 10 - Lbar)^2))) %>%
    dplyr::filter(sd_Lbar>=0.01) -> laa_stats

  df <- list(ages = laa_stats$age,
             Lbar = laa_stats$Lbar,
             sd_Lbar = laa_stats$sd_Lbar)

  par <- list(log_Linf = log(max(laa_stats$Lbar)),
              log_K = log(0.2),
              t0 = -0.1,
              log_alpha = log(0.07),
              log_beta = log(2.07))


  f <- function(par) {

    RTMB::getAll(par, df)
    Linf = exp(log_Linf)
    k = exp(log_K)
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    n = 1:length(ages)

    pred = Linf * (1 - exp(-k * (ages - t0)))
    yvar = log(1. + sd_Lbar^2 / Lbar^2)
    yconst = log(2.0 * pi * yvar * Lbar^2)
    rss = 0.5 * (yconst + (log(pred) - log(Lbar))^2 / yvar)

    lpred = alpha * log(ages) + beta
    lnll = sqrt(n) * (log(lpred) - log(sd_Lbar))^2
    nll = sum(rss) + sum(lnll)
    RTMB::REPORT(pred)
    RTMB::REPORT(Linf)
    RTMB::REPORT(k)
    RTMB::REPORT(t0)
    RTMB::REPORT(alpha)
    RTMB::REPORT(beta)

    return(nll)

  }
  # run model
  obj <- RTMB::MakeADFun(f, par)
  fit <- nlminb(obj$par, obj$fn, obj$gr)
  sd <- RTMB::sdreport(obj)
  report <- obj$report(obj$env$last.par.best)

  Linf = report$Linf
  k = report$k
  t0 = report$t0
  alpha = report$alpha
  beta = report$beta
  Lbar = df$Lbar
  sd_Lbar = df$sd_Lbar
  ages = df$ages

  # Compute Sz@A transition matrix

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

  list(laa_stats = laa_stats, sd = sd, report = report)

}
