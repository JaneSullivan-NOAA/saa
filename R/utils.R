prep_alw_data <- function(age_data, length_data, age_error, len_bins, rec_age) {
  ages = rec_age:(rec_age + nrow(age_error) - 1)
  age_data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year>=1990, age %in% rec_age:nrow(age_error))  %>%
    dplyr::select(age, length) %>%
    tidytable::filter(.N>1, .by = age) %>%
    tidytable::mutate(n_l = .N, .by = length) %>%
    tidytable::mutate(sample_size = .N, .by = age) %>%
    tidytable::arrange(age, length) -> inter

  length_data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year>=1990, !is.na(length)) -> dat

  if(!("frequency" %in% names(dat))){
    dat %>%
      tidytable::summarise(tot = .N, .by = length) -> dat
  } else {
    dat %>%
      tidytable::summarise(tot = sum(frequency), .by = length) -> dat
  }

  # size at age
  dat %>%
    dplyr::left_join(inter, .) %>%
    tidytable::mutate(prop =  dplyr::n() / n_l * tot, .by = c(age, length)) %>%
    tidytable::distinct() %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(sample_size = mean(sample_size),
                     lbar = sum(prop * length) / sum(prop) * 0.1,
                     sd = sqrt(1 / (sum(prop) - 1) * sum(prop * (length / 10 - lbar)^2))) %>%
    tidytable::filter(sd>=0.01) -> laa_stats


  age_data %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(year >= 1990, !is.na(age)) %>%
    dplyr::select(age, length, weight) %>%
    tidytable::filter(.N>1, .by = age) -> a1

  # length-weight
  a1 %>%
    dplyr::filter(length > 0, !is.na(weight)) %>%
    tidytable::summarise(wbar = mean(weight, na.rm = T),
                         sd = sd(weight, na.rm = T), .by = length) %>%
    tidytable::drop_na() -> lw_data

  # weight at age
  a1 %>%
    tidytable::mutate(sample_size = .N, .by = age) %>%
    tidytable::mutate(n_al = .N, .by = c(age, length)) %>%
    tidytable::mutate(n_l = .N, .by = length) %>%
    tidytable::filter(n_l>1) %>%
    tidytable::arrange(length) %>%
    tidytable::left_join(dat %>%
                dplyr::mutate(alpha_l = tot / sum(tot))) %>%
    tidytable::mutate(N_al = n_al / n_l * tot) %>%
    tidytable::filter(!is.na(weight)) %>%
    tidytable::mutate(wbar_la = mean(weight, na.rm=T),
                      v_wbar_la = var(weight, na.rm=T) / dplyr::n(),
                      .by = c(age, length)) %>%
    tidytable::mutate(theta_la = n_al / n_l,
                      r_la = alpha_l * theta_la) %>%
    tidytable::select(-weight) %>%
    tidytable::distinct() %>%
    tidytable::mutate(theta_a = sum(r_la), .by = age) -> inter1

  inter1 %>%
    tidytable::distinct(length, tot) %>%
    tidytable::summarise(L = sum(tot)) %>%
    tidytable::pull(L) -> L

  inter1 %>%
    dplyr::mutate(v_r_la = alpha_l^2 * theta_la * (1 - theta_la) / (n_l - 1) + alpha_l * (theta_la - theta_a)^2 / L) %>%
    tidytable::mutate(wbar = sum(r_la * wbar_la, na.rm=TRUE) / sum(r_la),
                      sd = sqrt(sum(r_la^2 * v_wbar_la + (wbar_la - wbar)^2 * v_r_la, na.rm=TRUE) /
                                  theta_a^2) *
                        sqrt(sample_size), .by = age) %>%
    dplyr::select(age, sample_size, wbar, sd) %>%
    dplyr::distinct() %>%
    tidytable::filter(sample_size >= 30) %>%
    tidytable::arrange(age) -> waa_stats

  list(dat_saa = list(age = laa_stats$age,
                      n = laa_stats$sample_size,
                      lbar = laa_stats$lbar,
                      sd = laa_stats$sd),
       dat_lw = list(length = lw_data$length,
                     wbar = lw_data$wbar,
                     sd = lw_data$sd),
       dat_waa = list(age = waa_stats$age,
                      wbar = waa_stats$wbar,
                      sd = waa_stats$sd),
       laa_stats = laa_stats,
       lw_data = lw_data,
       waa_stats = waa_stats)


}
