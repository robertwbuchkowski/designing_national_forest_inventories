# Functions

bwvr_alex <- function(data, vars, unit_col) {
  purrr::map_dfr(vars, function(v) {
    d <- data %>%
      dplyr::select(unit = !!rlang::sym(unit_col), value = !!rlang::sym(v)) %>%
      dplyr::filter(is.finite(value))

    # per-unit stats
    agg <- d %>%
      dplyr::group_by(unit) %>%
      dplyr::summarise(
        n_obs_per_unit  = dplyr::n(),
        mean_val  = mean(value),
        sd = sd(value),
        ss = sum((value - mean(value))^2), # sum of squares
        .groups = "drop"
      )

    if (nrow(agg) < 2) {
      return(tibble::tibble(
        Variable = v,
        Overall_Mean_raw = NA,
        Overall_Mean_weighted = NA,
        Between_SD_Raw    = NA,
        Within_SD_Raw    = NA,
        Between_SD_Pooled = NA,
        Within_SD_Pooled  = NA
      ))
    }

    # Raw variances
    raw_between <- sd(agg$mean_val, na.rm = TRUE)
    raw_within  <- mean(agg$sd, na.rm = TRUE) # simple average of within-unit SD

    # Pooled variances
    k <- nrow(agg)# total number of units
    N <- sum(agg$n_obs_per_unit) # total number of replicates
    mraw <- mean(agg$mean_val)
    mbar <- sum(agg$n_obs_per_unit * agg$mean_val) / N # weighted mean of all replicates (if unequal N)


    SSE <- sum(agg$ss) # sum  of squares error (within): total within unit variation
    SSB <- sum(agg$n_obs_per_unit * (agg$mean_val - mbar)^2) # sum of squares between (how far each units mean is from grand mean, weighted by unit sample size)
    MSW <- SSE / (N - k) # mean square within (pooled within-unit variance)
    MSB <- SSB / (k - 1) # mean square between (variance of unit means)
    n_eff <- (N - sum(agg$n_obs_per_unit^2) / N) / (k - 1) # effetive sample size per group (important when groups are unbalanced)

    pooled_within  <- sqrt(MSW) # pooled within unit sd
    pooled_between <- sqrt(pmax(0, (MSB - MSW) / n_eff)) # de-noised between-unit SD// THIS is SAME AS RANDOM-INTERCEPT MODEL (Supposed to be)

    tibble::tibble(
      Variable       = v,
      Overall_Mean_Raw = mraw,
      Overall_Mean_Weighted = mbar,
      Between_SD_Raw    = raw_between,
      Within_SD_Raw    = raw_within,
      Between_SD_Pooled = pooled_between,
      Within_SD_Pooled  = pooled_within,
      number_plots = k,
      number_micrplots = N
    )
  })
}


# MDC Simple --------------------------------------------------------------

# This function uses a temporal dataset with exactly two time points
# Computes the delta, simply gets the sd(delta) and uses that in equation for MDC

# Ignores any within site variance

simple_mdc <- function(df,
                       unit_col = "unit",
                       time_col = "time",
                       value_col = "value",
                       alpha = 0.05) {
  stopifnot(all(c(unit_col, time_col, value_col) %in% names(df)))

  # Pivot to wide: one row per unit, columns for T0 and T1
  wide <- df %>% select(!!sym(unit_col), !!sym(time_col),!!sym(value_col)) %>%
    tidyr::pivot_wider(names_from = !!sym(time_col), values_from = !!sym(value_col)) %>%
    drop_na()

  # Get time names (assuming just 2)
  time_names <- setdiff(names(wide), unit_col)
  if (length(time_names) != 2) stop("Need exactly two time points")

  # Calculate change per unit
  delta <- wide[[time_names[2]]] - wide[[time_names[1]]]

  # MDC = z * SD(delta) / sqrt(n)
  n <- length(delta)
  z_alpha <- qnorm(1 - alpha / 2)
  sd_delta <- sd(delta)

  mdc <- z_alpha * sd_delta / sqrt(n)  # Option 1: basic
  mdc_alt <- z_alpha * sd_delta * sqrt(2 / n)  # Option 2: from some papers, not sure why exactly

  tibble::tibble(
    n_units = n,
    mean_change = mean(delta),
    sd_delta = sd_delta,
    alpha = alpha,
    z_alpha = z_alpha,
    MDC = mdc,
    MDC_sqrt2 = mdc_alt
  )
}


# MDC from single time point (within and between) -------------------------

# with calcualted values from bwvr_alex
mdc_from_bwvr <- function(bwvr_out,
                                 alpha = 0.05,
                                 power = 0.80,
                                 use_t = TRUE) {

  k <- as.integer(bwvr_out$number_plots)
  N <- as.integer(bwvr_out$number_micrplots)
  if (k < 2) stop("n_units must be >= 2.")

  purrr::pmap_dfr(
    list(
      Variable = bwvr_out$Variable,
      grand_mean = bwvr_out$Overall_Mean_Raw,
      sigma_b = bwvr_out$Between_SD_Pooled,
      sigma_w = bwvr_out$Within_SD_Pooled
    ),
    function(Variable, grand_mean, sigma_b, sigma_w) {

      within_term <- if (is.finite(sigma_w)) (sigma_w^2) else 0
      var_mean <- (sigma_b^2)/k + within_term/N
      SE_pop <- sqrt(pmax(var_mean, 0))

      crit_alpha <- if (use_t) stats::qt(1 - alpha/2, df = k - 1) else stats::qnorm(1 - alpha/2)
      crit_beta  <- if (is.null(power)) 0 else stats::qnorm(power)

      tibble::tibble(
        Variable            = Variable,
        m_units             = k,
        grand_mean          = grand_mean,
        within_sd_pooled    = sigma_w,
        between_sd_pooled   = sigma_b,
        SE_pop_change       = SE_pop,
        MDC_pop_alpha       = crit_alpha * SE_pop,
        MDD_pop_alpha_power = (crit_alpha + crit_beta) * SE_pop,
        alpha               = alpha,
        power               = power
      )
    }
  )
}


# MDC propagate uncertainty  ----------------------------------------------

# temporal within and between variation is not scaled by year...


# This is a more complicated version, a bit stitched together for incorperating within site variance
# calculates within-unit variance
# also calculates temporal correlations of paired unit (site) or core using a pooling method which weights by sample size
  # The rho can be used for simulations or other stuff; also incorporated into more comprehensive MDC calculation; rho = 0 is more conservative


# the main formula is
  # MDC = z_alpha * SE_pop
  # where,
  # SE_pop <- sqrt((sigma2_delta_true + meas_var_bar) / m)
  # sigma2_delta_tru is variance of unit-level deltas
  # meas_var_bar is the average variance per delta

# meas_var_bar is computed with the full variance equation of
  # two correlated means: Var(X-Y) = Var(X) + Var(Y) - 2*Cov(X,Y)


mdc_pop_full <- function(df,
                         unit_col   = "nfi_plot",
                         time_col   = "meas_num",
                         value_col  = "TC",
                         repl_col   = NULL,     # core/pit id to pair cores across time to calculate rho
                         alpha = 0.05,
                         power = 0.80) {

  stopifnot(all(c(unit_col, time_col, value_col) %in% names(df)))

  d0 <- df %>%
    dplyr::transmute(
      unit  = .data[[unit_col]],
      time  = .data[[time_col]],
      value = .data[[value_col]],
      repl  = if (!is.null(repl_col) && repl_col %in% names(df)) .data[[repl_col]] else NA
    )

  # Exactly two time points
  tt <- sort(stats::na.omit(unique(d0$time)))
  if (length(tt) != 2) stop("Need exactly two time points in ", time_col)
  t0 <- tt[1]; t1 <- tt[2]

  # For each unit, calculate mean, and within variance
  ut <- d0 %>%
    dplyr::group_by(unit, time) %>%
    dplyr::summarise(
      n_core = dplyr::n(),
      mean_y = mean(value, na.rm = TRUE),
      var_w  = stats::var(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = time,
      values_from = c(n_core, mean_y, var_w),
      names_sep = "_T"
    )

  req_cols <- c(paste0("n_core_T", t0), paste0("n_core_T", t1),
                paste0("mean_y_T", t0), paste0("mean_y_T", t1),
                paste0("var_w_T",  t0), paste0("var_w_T",  t1))
  if (!all(req_cols %in% names(ut))) {
    stop("Some units are missing one of the time points; cannot pair.")
  }

  paired <- ut %>%
    dplyr::filter(!is.na(.data[[paste0("mean_y_T", t0)]]) & !is.na(.data[[paste0("mean_y_T", t1)]]))
  m <- nrow(paired)
  if (m < 2) stop("Need at least two paired units.")

  n0 <- paired[[paste0("n_core_T", t0)]]
  n1 <- paired[[paste0("n_core_T", t1)]]
  y0 <- paired[[paste0("mean_y_T", t0)]]
  y1 <- paired[[paste0("mean_y_T", t1)]]
  v0 <- paired[[paste0("var_w_T",  t0)]]
  v1 <- paired[[paste0("var_w_T",  t1)]]

  delta <- y1 - y0 # calculate delta

    # Plot-level correlation of means (always available)
  rho_plot <- suppressWarnings(stats::cor(y0, y1, use = "complete.obs")) # this should be like rho_unit to keep names consistent

  # Core-level rho (only if real core IDs exist at both times)
  rho_core <- NA_real_
  n_core_pairs <- 0L
  n_units_with_pairs <- 0L

  if (!all(is.na(d0$repl))) {

    # paired cores within each unit
    core_pairs <- d0 %>%
      dplyr::filter(time %in% tt) %>%
      tidyr::pivot_wider(names_from = time, values_from = value, names_prefix = "T_") %>%
      dplyr::filter(is.finite(.data[[paste0("T_", t0)]]) & is.finite(.data[[paste0("T_", t1)]]))

    if (nrow(core_pairs) >= 3) {
      # Per-unit (pooled) correlations, Fisher-z pooled
      per_unit <- core_pairs %>%
        dplyr::group_by(unit) %>%
        dplyr::summarise(
          n_pairs = dplyr::n(),
          r       = suppressWarnings(stats::cor(.data[[paste0("T_", t0)]],
                                                .data[[paste0("T_", t1)]],
                                                use = "complete.obs")),
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(r), n_pairs >= 2)

      if (nrow(per_unit) > 0) {
        # Fisher z pooling weighted by (n_pairs - 3), clamp r to avoid infinities
        r_clamp <- pmin(pmax(per_unit$r, -0.999999), 0.999999)
        z <- atanh(r_clamp)
        w <- pmax(per_unit$n_pairs - 3, 1)  # weights
        rho_core <- tanh(stats::weighted.mean(z, w))
        n_core_pairs <- sum(per_unit$n_pairs)
        n_units_with_pairs <- nrow(per_unit)
      }
    }
  }

  # Which rho to use in measurement variance (0 if no paired cores)
  rho_for_meas <- if (is.finite(rho_core)) rho_core else 0

  # variance of delta per unit; difference between two means ()
  # formula for variance of two correlated means:
    # Var(X-Y) = Var(X) + Var(Y) - 2*Cov(X,Y)
    # Where, Covariance of sample means Cov(X,Y) = rho*sqrt(v1/n1*v0/no)

  meas_var_i <- (v1 / n1) + (v0 / n0) - 2 * rho_for_meas * sqrt((v1 / n1) * (v0 / n0))
  meas_var_i[meas_var_i < 0] <- 0
  meas_var_bar <- mean(meas_var_i, na.rm = TRUE) # this is what goes directly into the MDC/MDD

  # variance
  s2_delta_obs <- stats::var(delta, na.rm = TRUE)
  sigma2_delta_true <- max(0, s2_delta_obs - meas_var_bar) # this goes into for SE of population.



  # Using lmm for estimating effect size
  # also using it as a simpler way to get decomposed within and between SD

  d_long <- d0 %>%
    dplyr::filter(time %in% tt) %>%
    dplyr::mutate(time = factor(time, levels = tt))

  fit <- lme4::lmer(value ~ time + (1 | unit), data = d_long, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(fit))
  within_sd  <- sqrt(vc$vcov[vc$grp == "Residual"])
  between_sd <- sqrt(vc$vcov[vc$grp == "unit" & vc$var1 == "(Intercept)"])

  # Fixed effect for second time (T1 vs T0)
  fx <- broom.mixed::tidy(fit, effects = "fixed", conf.int = FALSE)
  term_T1 <- paste0("time", tt[2])
  eff_row <- fx[fx$term == term_T1, , drop = FALSE]
  if (nrow(eff_row) == 0) {
    # fallback: last time term
    eff_row <- fx[grepl("^time", fx$term), , drop = FALSE]
    eff_row <- eff_row[nrow(eff_row), , drop = FALSE]
  }
  est  <- eff_row$estimate[1]
  se   <- eff_row$std.error[1]
  pval <- if ("p.value" %in% names(eff_row)) eff_row$p.value[1] else 2 * stats::pnorm(-abs(est / se))

  # Population level estimate of MDC and MDD
  SE_pop <- sqrt((sigma2_delta_true + meas_var_bar) / m)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- if (is.null(power)) 0 else stats::qnorm(power)

  tibble::tibble(
    m_units             = m,
    time0               = t0,
    time1               = t1,
    mean_change         = mean(delta, na.rm = TRUE),
    lmm_time_est        = est,
    lmm_time_se         = se,
    lmm_time_p          = pval,
    within_sd_pooled    = within_sd,
    between_sd_pooled   = between_sd,
    sd_delta_obs        = sqrt(s2_delta_obs),
    meas_var_bar        = meas_var_bar,
    # explicitly outputting variation of change within plots
    within_sd_delta = sqrt(meas_var_bar), # previously: meas_var_bar
    between_sd_delta       = sqrt(sigma2_delta_true), # previously: sd_delta_true
    SE_pop_change       = SE_pop,
    MDC_pop_alpha       = z_alpha * SE_pop,
    MDD_pop_alpha_power = (z_alpha + z_beta) * SE_pop,
    rho_core            = rho_core,              # NA if not enough paired cores
    rho_plot            = rho_plot,
    rho_used_for_meas   = rho_for_meas,         # 0 if cores not remeasured
    n_core_pairs        = n_core_pairs,
    n_units_with_pairs  = n_units_with_pairs,
    alpha               = alpha,
    power               = power
  )
}


# MDD with per year scaling for temporal ----------------------------------

mdc_pop_full_yr <- function(df,
                            unit_col   = "nfi_plot",
                            time_col   = "meas_num",
                            date_col   = "meas_date",
                            value_col  = "TC",
                            repl_col   = NULL,      # core/pit id to pair cores across time
                            alpha = 0.05,
                            power = 0.80,
                            horizon_years = NULL) {  # e.g., 10 for “10-year MDD”
  stopifnot(all(c(unit_col, time_col, date_col, value_col) %in% names(df)))

  d0 <- df %>%
    dplyr::transmute(
      unit  = .data[[unit_col]],
      time  = .data[[time_col]],
      date  = as.Date(.data[[date_col]]),
      value = .data[[value_col]],
      repl  = if (!is.null(repl_col) && repl_col %in% names(df)) .data[[repl_col]] else NA
    )

  # Exactly two time points overall (NFI style); paired within unit
  tt <- sort(stats::na.omit(unique(d0$time)))
  if (length(tt) != 2) stop("Need exactly two time points in ", time_col)
  t0 <- tt[1]; t1 <- tt[2]

  # Per-unit summaries at each time (including within-plot variance)
  ut <- d0 %>%
    dplyr::group_by(unit, time) %>%
    dplyr::summarise(
      n_core = dplyr::n(),
      mean_y = mean(value, na.rm = TRUE),
      var_w  = stats::var(value, na.rm = TRUE),
      date   = max(date, na.rm = TRUE),  # if multiple cores, they share date; take max safely
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = time,
      values_from = c(n_core, mean_y, var_w, date),
      names_sep = "_T"
    )

  req_cols <- c(paste0("n_core_T", t0), paste0("n_core_T", t1),
                paste0("mean_y_T", t0), paste0("mean_y_T", t1),
                paste0("var_w_T",  t0), paste0("var_w_T",  t1),
                paste0("date_T",   t0), paste0("date_T",   t1))
  if (!all(req_cols %in% names(ut))) {
    stop("Some units are missing one of the time points; cannot pair.")
  }

  paired <- ut %>%
    dplyr::filter(!is.na(.data[[paste0("mean_y_T", t0)]]) & !is.na(.data[[paste0("mean_y_T", t1)]]) &
                    !is.na(.data[[paste0("date_T", t0)]])   & !is.na(.data[[paste0("date_T", t1)]]))

  m <- nrow(paired)
  if (m < 2) stop("Need at least two paired units.")

  n0 <- paired[[paste0("n_core_T", t0)]]
  n1 <- paired[[paste0("n_core_T", t1)]]
  y0 <- paired[[paste0("mean_y_T", t0)]]
  y1 <- paired[[paste0("mean_y_T", t1)]]
  v0 <- paired[[paste0("var_w_T",  t0)]]
  v1 <- paired[[paste0("var_w_T",  t1)]]
  d0_date <- paired[[paste0("date_T", t0)]]
  d1_date <- paired[[paste0("date_T", t1)]]

  # interval in decimal years
  dt_years <- as.numeric(difftime(d1_date, d0_date, units = "days")) / 365.25
  if (any(dt_years <= 0 | !is.finite(dt_years))) stop("Non-positive or invalid remeasurement interval found.")

  delta <- y1 - y0
  rate  <- delta / dt_years  # per-year change for each unit

  # Plot-level correlation of means (optional diagnostic)
  rho_plot <- suppressWarnings(stats::cor(y0, y1, use = "complete.obs"))

  # Core-level rho (only if real core IDs exist at both times)
  rho_core <- NA_real_
  n_core_pairs <- 0L
  n_units_with_pairs <- 0L

  if (!all(is.na(d0$repl))) {
    core_pairs <- d0 %>%
      dplyr::filter(time %in% tt, !is.na(repl)) %>%
      dplyr::group_by(unit, repl, time) %>%
      dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(
        id_cols     = c(unit, repl),
        names_from  = time,
        values_from = value,
        names_prefix = "T_",
        values_fill = NA_real_
      ) %>%
      dplyr::filter(
        is.finite(.data[[paste0("T_", t0)]]) &
          is.finite(.data[[paste0("T_", t1)]])
      )

    if (nrow(core_pairs) >= 3) {
      per_unit <- core_pairs %>%
        dplyr::group_by(unit) %>%
        dplyr::summarise(
          n_pairs = dplyr::n(),
          r       = suppressWarnings(stats::cor(.data[[paste0("T_", t0)]],
                                                .data[[paste0("T_", t1)]],
                                                use = "complete.obs")),
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(r), n_pairs >= 2)

      if (nrow(per_unit) > 0) {
        r_clamp <- pmin(pmax(per_unit$r, -0.999999), 0.999999)
        z <- atanh(r_clamp)
        w <- pmax(per_unit$n_pairs - 3, 1)
        rho_core <- tanh(stats::weighted.mean(z, w))
        n_core_pairs <- sum(per_unit$n_pairs)
        n_units_with_pairs <- nrow(per_unit)
      }
    }
  }

  # Which rho to use in measurement variance (0 if no paired cores)
  rho_for_meas <- if (is.finite(rho_core)) rho_core else 0

  # Measurement variance of the DIFFERENCE, then convert to rate (divide by dt^2)
  meas_var_delta_i <- (v1 / n1) + (v0 / n0) - 2 * rho_for_meas * sqrt((v1 / n1) * (v0 / n0))
  meas_var_delta_i[meas_var_delta_i < 0] <- 0
  meas_var_rate_i <- meas_var_delta_i / (dt_years^2)
  meas_var_rate_bar <- mean(meas_var_rate_i, na.rm = TRUE)

  # Observed variance of rates across units; back out true between-unit variance in rates
  s2_rate_obs <- stats::var(rate, na.rm = TRUE)
  sigma2_rate_true <- max(0, s2_rate_obs - meas_var_rate_bar)

  # LMM on irregular timing: value ~ years_since_baseline + (1|unit)
  d_long <- d0 %>%
    dplyr::filter(time %in% tt) %>%
    dplyr::group_by(unit) %>%
    dplyr::mutate(t_years = as.numeric(difftime(date, min(date, na.rm = TRUE), units = "days"))/365.25) %>%
    dplyr::ungroup()

  fit <- lme4::lmer(value ~ t_years + (1 | unit), data = d_long, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(fit))
  within_sd  <- sqrt(vc$vcov[vc$grp == "Residual"])
  between_sd <- sqrt(vc$vcov[vc$grp == "unit" & vc$var1 == "(Intercept)"])

  fx <- broom.mixed::tidy(fit, effects = "fixed", conf.int = FALSE)
  est_rate  <- fx$estimate[fx$term == "t_years"][1]
  se_rate   <- fx$std.error[fx$term == "t_years"][1]
  pval_rate <- if ("p.value" %in% names(fx)) fx$p.value[fx$term == "t_years"][1] else 2 * stats::pnorm(-abs(est_rate / se_rate))

  # SE of population mean rate
  SE_pop_rate <- sqrt((sigma2_rate_true + meas_var_rate_bar) / m)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- if (is.null(power)) 0 else stats::qnorm(power)

  out <- tibble::tibble(
    m_units               = m,
    time0                 = t0,
    time1                 = t1,
    mean_dt_years         = mean(dt_years),
    median_dt_years       = stats::median(dt_years),
    min_dt_years          = min(dt_years),
    max_dt_years          = max(dt_years),
    mean_change           = mean(delta, na.rm = TRUE),
    mean_rate             = mean(rate, na.rm = TRUE),        # per-year mean change
    lmm_slope_est_rate    = est_rate,                         # per-year slope from LMM
    lmm_slope_se          = se_rate,
    lmm_slope_p           = pval_rate,
    within_sd_pooled      = within_sd,
    between_sd_pooled     = between_sd,
    sd_rate_obs           = sqrt(s2_rate_obs),
    meas_var_rate_bar     = meas_var_rate_bar,
    within_sd_rate        = sqrt(meas_var_rate_bar),          # measurement SD component (per year)
    between_sd_rate       = sqrt(sigma2_rate_true),           # true between-unit SD of rate
    SE_pop_rate           = SE_pop_rate,
    MDC_rate_pop_alpha    = z_alpha * SE_pop_rate,            # annual MDC
    MDD_rate_pop_a_p      = (z_alpha + z_beta) * SE_pop_rate, # annual MDD
    rho_core              = rho_core,
    rho_plot              = rho_plot,
    rho_used_for_meas     = rho_for_meas,
    n_core_pairs          = n_core_pairs,
    n_units_with_pairs    = n_units_with_pairs,
    alpha                 = alpha,
    power                 = power
  )

  if (!is.null(horizon_years) && is.finite(horizon_years) && horizon_years > 0) {
    out <- out %>%
      dplyr::mutate(
        MDC_horizon_alpha    = MDC_rate_pop_alpha * horizon_years,
        MDD_horizon_a_p      = MDD_rate_pop_a_p   * horizon_years,
        horizon_years        = horizon_years
      )
  }

  out
}



# Simulation Function -----------------------------------------------------


simulate_two_time <- function(n_sites = 100,
                              n_rep   = 4,            # cores per site per time
                              sigma_b = 15,           # between-site SD
                              sigma_w = 8,            # within-site SD
                              delta_mean = 0,         # true mean change T2 - T1
                              sigma_delta = 2,
                              paired = TRUE,          # TRUE = remeasure same cores (positive rho)
                              rho_core = 0.5,         # correlation of core-level error across time if paired
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # True site means at T1
  mu_site <- rnorm(n_sites, mean = 0, sd = sigma_b)

  delta_i <- rnorm(n_sites, mean = delta_mean, sd = sigma_delta)
  # Correlated core errors across time if paired; otherwise independent
  # We'll generate for each site x core: (e1, e2) with Corr = rho_core
  # Cholesky of 2x2 corr matrix:
  L <- matrix(c(1, 0, rho_core, sqrt(1 - rho_core^2)), 2, 2, byrow = TRUE)

  # Build data
  df_list <- vector("list", n_sites)
  for (i in seq_len(n_sites)) {
    # core-level errors
    if (paired) {
      z <- matrix(rnorm(2 * n_rep), nrow = 2) # 2 x n_rep
      e12 <- L %*% z                          # induces correlation
      e1  <- e12[1, ] * sigma_w
      e2  <- e12[2, ] * sigma_w
    } else {
      e1 <- rnorm(n_rep, sd = sigma_w)
      e2 <- rnorm(n_rep, sd = sigma_w)
    }



    y1 <- mu_site[i] + e1                     # T1 cores
    y2 <- (mu_site[i] + delta_i[i]) + e2      # T2 cores

    df_list[[i]] <- tibble(
      site  = i,
      core  = seq_len(n_rep),
      T1    = y1,
      T2    = y2
    )
  }

  bind_rows(df_list)
}

print("One time piont parameter function: bwvr_alex(data, vars, unit_col)")
print("Resample Parameter Function: mdc_pop_full(df, unit_col, time_col, value_col, repl_col, alpha, power)")
print("Simulation function: simulate_two_time(n_sites, n_rep, sigma_b, sigma_w, delta_mean, paired, rho_core, seed)")
