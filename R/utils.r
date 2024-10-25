process_extracted_data <- function(n, d, se, sp, ppv, npv, name = NA_character_) {
    tp <- round(se * d)
    tn <- round(sp * (n - d))
    fn <- d - tp
    fp <- round((1 - ppv) * (se * d) / ppv)
    se_hat <- binom::binom.exact(tp, d)
    sp_hat <- binom::binom.exact(tn, n - d)
    ppv_hat <- binom::binom.exact(tp, tp + fp)
    npv_hat <- binom::binom.exact(tn, tn + fn)
    out <- data.frame(
        name = name, n = n, d = d,
        tp = tp, tn = tn, fn = fn, fp = fp,
        p = round(d / n * 100, 1),
        se_reported = round(se * 100, 1),
        se_hat = round(se_hat$mean * 100, 1),
        se_lower = round(se_hat$lower * 100, 1),
        se_upper = round(se_hat$upper * 100, 1),
        sp_reported = round(sp * 100, 1),
        sp_hat = round(sp_hat$mean * 100, 1),
        sp_lower = round(sp_hat$lower * 100, 1),
        sp_upper = round(sp_hat$upper * 100, 1),
        ppv_reported = round(ppv * 100, 1),
        ppv_hat = round(ppv_hat$mean * 100, 1),
        ppv_lower = round(ppv_hat$lower * 100, 1),
        ppv_upper = round(ppv_hat$upper * 100, 1),
        npv_reported = round(npv * 100, 1),
        npv_hat = round(npv_hat$mean * 100, 1),
        npv_lower = round(npv_hat$lower * 100, 1),
        npv_upper = round(npv_hat$upper * 100, 1)
    )
    return(out)
}

generate_data_from_counts <- function(n, d, tp, tn, fp, fn) {
    # check if things match up
    total <- tp + tn + fp + fn
    if (n != total) {
        msg <- stringr::str_glue(
            "n={n} but tp + tn + fp + fn = {total}"
        )
        stop(msg)
    }
    total_d <- tp + fn
    if (d != total_d) {
        msg <- stringr::str_glue(
            "d={d} but tp + fn = {total_d}"
        )
        stop(msg)
    }

    # generate the data for each outcome type
    tp_data <- data.frame(
        outcomes = rep(1, tp),
        mced_test = rep(1, tp)
    )
    tn_data <- data.frame(
        outcomes = rep(0, tn),
        mced_test = rep(0, tn)
    )
    fp_data <- data.frame(
        outcomes = rep(0, fp),
        mced_test = rep(1, fp)
    )
    fn_data <- data.frame(
        outcomes = rep(1, fn),
        mced_test = rep(0, fn)
    )

    # combine all data
    out <- dplyr::bind_rows(
        tp_data,
        tn_data,
        fp_data,
        fn_data
    )
    return(out)
}

compute_nb_untreated <- function(
    bdca_fit,
    decision_strategy = "mced_test",
    return_posterior = FALSE) {
    se <- bdca_fit$fit$distributions$sensitivity[[decision_strategy]]
    sp <- bdca_fit$fit$distributions$specificity[[decision_strategy]]
    p <- bdca_fit$fit$distributions$prevalence

    trade_off_ut <- delta_ut <- treat_none_ut <- nb_mced_ut <- matrix(
        nrow = nrow(se),
        ncol = ncol(se)
    )

    for (i in seq_along(bdca_fit$thresholds)) {
        .t <- bdca_fit$thresholds[i]
        w_t <- (1 - .t) / .t
        # Net TN = Sp * (1 - p) - 1/odds(t) * (1 - Se) * p
        nb_mced_ut[, i] <- sp[, i] * (1 - p) - w_t * (1 - se[, i]) * p
        treat_none_ut[, i] <- 1 * (1 - p) - w_t * (1 - 0) * p
        treat_all_ut <- 0.00
        delta_ut[, i] <- nb_mced_ut[, i] - pmax(treat_none_ut[, i], treat_all_ut)
        trade_off_ut[, i] <- 1 / delta_ut[, i]
    }

    df <- data.frame(
        thr = bdca_fit$thresholds,
        estimate = colMeans(nb_mced_ut),
        lower = matrixStats::colQuantiles(nb_mced_ut, probs = 0.025),
        upper = matrixStats::colQuantiles(nb_mced_ut, probs = 0.975),
        estimate_none = colMeans(treat_none_ut),
        lower_none = matrixStats::colQuantiles(treat_none_ut, probs = 0.025),
        upper_none = matrixStats::colQuantiles(treat_none_ut, probs = 0.975),
        estimate_delta = colMeans(delta_ut),
        lower_delta = matrixStats::colQuantiles(delta_ut, probs = 0.025),
        upper_delta = matrixStats::colQuantiles(delta_ut, probs = 0.975),
        p_useful = colMeans(delta_ut > 0),
        estimate_trade_off_negative = matrixStats::colQuantiles(trade_off_ut, probs = 0.5),
        lower_trade_off_negative = matrixStats::colQuantiles(trade_off_ut, probs = 0.025),
        upper_trade_off_negative = matrixStats::colQuantiles(trade_off_ut, probs = 0.975)
    )

    if (isFALSE(return_posterior)) {
        return(df)
    } else {
        return(list(
            df = df,
            posterior = list(
                ntn = nb_mced_ut,
                treat_none = treat_none_ut,
                delta = delta_ut,
                treat_all = 0
            )
        ))
    }
}

get_population_scaling_factor <- function(as_string = FALSE) {
    if (isTRUE(as_string)) {
        return("100,000")
    } else {
        return(100000)
    }
}


compute_trade_off_positive <- function(fit) {
    nb <- fit$fit$distributions$net_benefit$mced_test
    treat_all <- fit$fit$distributions$treat_all
    treat_none <- 0.0
    trade_off_positive_posterior <- matrix(
        nrow = nrow(nb), ncol = ncol(nb)
    )
    for (i in seq_along(fit$thresholds)) {
        .t <- fit$thresholds[i]
        .nb <- nb[, i]
        .default_nb <- pmax(treat_all[, i], treat_none)
        .delta_nb <- .nb - .default_nb
        trade_off_positive_posterior[, i] <- 1 / .delta_nb
    }

    trade_off_positive_summary <- data.frame(
        thr = fit$thresholds,
        estimate_trade_off_positive = matrixStats::colQuantiles(trade_off_positive_posterior, probs = 0.5),
        lower_trade_off_positive = matrixStats::colQuantiles(trade_off_positive_posterior, probs = 0.025),
        upper_trade_off_positive = matrixStats::colQuantiles(trade_off_positive_posterior, probs = 0.975)
    )

    return(trade_off_positive_summary)
}
