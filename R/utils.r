process_extracted_data <- function(n, d, se, sp, name = NA_character_) {
    tp <- round(se * d)
    tn <- round(sp * (n - d))
    fn <- d - tp
    fp <- (n - d) - tn
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
        ppv_hat = round(ppv_hat$mean * 100, 1),
        ppv_lower = round(ppv_hat$lower * 100, 1),
        ppv_upper = round(ppv_hat$upper * 100, 1),
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

compute_nb_untreated <- function(bdca_fit, decision_strategy = "mced_test") {
    nb_mced <- bdca_fit$fit$distributions$net_benefit[[decision_strategy]]
    nb_ta <- bdca_fit$fit$distributions$treat_all
    prev <- bdca_fit$fit$distributions$prevalence

    delta_ut <- treat_none_ut <- nb_mced_ut <- matrix(
        nrow = nrow(nb_mced),
        ncol = ncol(nb_mced)
    )
    for (i in seq_along(bdca_fit$thresholds)) {
        .t <- bdca_fit$thresholds[i]
        w_t <- .t / (1 - .t)
        nb_mced_ut[, i] <- 1000 * ((nb_mced[, i] - nb_ta[, i]) / w_t)
        treat_none_ut[, i] <- 1000 * ((1 - prev) - (prev * (1 / w_t)))
        delta_ut[, i] <- nb_mced_ut[, i] - pmax(treat_none_ut[, i], 0.0)
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
        p_useful = colMeans(delta_ut > 0)
    )
    return(df)
}

plot_net_benefit_treated <- function(bdca_fit, .color, .label) {
    bayesDCA:::plot.BayesDCA(
        bdca_fit,
        colors = list(mced_test = .color),
        labels = list(mced_test = .label)
    ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(legend.position = c(.85, .85))
}

plot_net_benefit_untreated <- function(bdca_fit, .color, .label) {
    plot_df <- compute_nb_untreated(bdca_fit)

    .colors <- c(.color, "gray40", "black")
    names(.colors) <- c(.label, "Refer all", "Refer none")

    plot_df %>%
        ggplot2::ggplot(ggplot2::aes(thr, estimate, ymin = lower, ymax = upper)) +
        # Plot curve for model/test
        ggplot2::geom_ribbon(alpha = 0.5, fill = .colors[.label]) +
        ggplot2::geom_line(linewidth = 1.5, ggplot2::aes(color = .label)) +
        # Plot curve for Refer none
        ggplot2::geom_ribbon(
            alpha = 0.5,
            ggplot2::aes(
                x = thr,
                ymin = lower_none,
                ymax = upper_none
            ),
            inherit.aes = FALSE,
            fill = .colors["Refer none"]
        ) +
        ggplot2::geom_line(
            linewidth = 1.5,
            ggplot2::aes(
                x = thr,
                y = estimate_none,
                color = "Refer none"
            ),
            inherit.aes = FALSE
        ) +
        # Plot curve for Refer all (i.e., zero net true negatives)
        ggplot2::geom_hline(
            linetype = 2,
            alpha = 0.7,
            linewidth = 1.5,
            ggplot2::aes(color = "Refer all", yintercept = 0)
        ) +
        ggplot2::coord_cartesian(ylim = c(-0.01, NA)) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Net true negatives\nper 1000 patients"
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            name = NULL
        )
}

plot_prob_useful <- function(bdca_fit, .color, .label) {
    bayesDCA:::plot_superiority_prob(
        bdca_fit,
        type = "useful",
        colors = list(mced_test = .color),
        labels = list(mced_test = .label)
    ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::labs(subtitle = NULL) +
        ggplot2::guides(color = "none")
}

plot_evpi <- function(bdca_fit, .color, .label) {
    suppressMessages({
        bayesDCA:::plot_evpi(
            bdca_fit,
            type = "useful", # same as type = "best" in this case, just using for plotting functionality
            colors = list(mced_test = .color),
            labels = list(mced_test = .label)
        ) +
            ggplot2::theme_bw(base_size = 18) +
            ggplot2::theme(legend.position = c(.75, .85)) +
            ggplot2::labs(subtitle = NULL, y = "Net true positives") +
            ggplot2::scale_y_continuous(
                breaks = scales::pretty_breaks(10),
                labels = function(x) 1000 * x
            ) +
            ggplot2::guides(color = "none")
    })
}
