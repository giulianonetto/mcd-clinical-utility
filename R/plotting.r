get_color_values <- function(.color, .label) {
    # your model, refer all, refer none
    marrom <- "#7e3f1a"
    cinza_escuro <- "gray20"
    cols <- c(rep(.color, 2), rep(cinza_escuro, 2), rep("gray60", 2))
    names(cols) <- c(.label, "mced_test", "Refer all", "Treat all", "Refer none", "Treat none")
    return(cols)
}
get_labels_values <- function(.label) {
    list(
        mced_test = .label,
        "Treat all" = "Refer all",
        "Treat none" = "Refer none"
    )
}
plot_net_benefit_treated <- function(fit, .color, .label, .pathway) {
    color_values <- get_color_values(.color = .color, .label = .label)
    label_values <- get_labels_values(.label = .label)

    if (stringr::str_to_lower(.pathway) != "overall") {
        ix <- stringr::str_detect(
            str_to_lower(names(color_values)),
            "all|none",
            negate = TRUE
        )
        .breaks <- names(color_values)[ix]
        if (stringr::str_to_lower(.pathway) == "lung") {
            legend_position <- c(.485, .2)
        } else {
            legend_position <- c(.485, .9)
        }
    } else {
        .breaks <- names(color_values)
        legend_position <- c(.485, .84)
    }
    p <- bayesDCA:::plot.BayesDCA(fit) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(
            legend.position = legend_position,
            legend.background = ggplot2::element_rect(fill = "transparent"),
            legend.key = ggplot2::element_rect(fill = "transparent"),
            legend.key.width = ggplot2::unit(2.5, "lines"),
            legend.spacing.y = ggplot2::unit(0.1, "cm")
        ) +
        ggplot2::scale_color_manual(
            values = color_values, labels = label_values,
            breaks = .breaks
        ) +
        ggplot2::scale_fill_manual(
            values = color_values, labels = label_values,
            breaks = .breaks
        ) +
        ggplot2::coord_cartesian(ylim = c(-0.001, NA)) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(
                byrow = TRUE,
                override.aes = list(linewidth = 2)
            ),
            fill = "none"
        )
    p$layers[[1]]$aes_params <- list(alpha = 0.5)
    return(p)
}

plot_delta_treated <- function(fit, .color, .label) {
    color_values <- get_color_values(.color = .color, .label = .label)
    label_values <- get_labels_values(.label = .label)
    p <- bayesDCA:::plot_delta(fit) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(legend.position = c(.25, .9)) +
        ggplot2::scale_color_manual(values = color_values, labels = label_values) +
        ggplot2::scale_fill_manual(values = color_values, labels = label_values) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        )
    return(p)
}

plot_net_benefit_untreated <- function(fit, .color, .label) {
    if (is.null(fit$nb_untreated)) {
        fit$nb_untreated <- compute_nb_untreated(fit)
    }

    if (is.data.frame(fit$nb_untreated)) {
        plot_df <- fit$nb_untreated
    } else {
        plot_df <- fit$nb_untreated$df
    }

    .colors <- get_color_values(.color = .color, .label = .label)
    label_function <- function(x) {
        x_lab <- paste0(x / 1000, "K")
        x_lab[x == 0L] <- "0"
        return(x_lab)
    }
    plot_df %>%
        dplyr::mutate_at(
            dplyr::vars(dplyr::matches("estimate|lower|upper")),
            ~ . * get_population_scaling_factor()
        ) %>%
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
            alpha = 1,
            linewidth = 1.5,
            ggplot2::aes(color = "Refer all", yintercept = 0)
        ) +
        ggplot2::coord_cartesian(ylim = c(-0.01, 9e4)) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        ) +
        ggplot2::scale_y_continuous(
            breaks = scales::pretty_breaks(10),
            labels = scales::label_number(big.mark = ",")
        ) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = stringr::str_glue(
                "Net true negatives per {get_population_scaling_factor(as_string = TRUE)} patients"
            )
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            name = NULL
        )
}

plot_prob_useful <- function(fit, .color, .label) {
    bayesDCA:::plot_superiority_prob(
        fit,
        type = "useful",
        colors = list(mced_test = .color),
        labels = list(mced_test = .label)
    ) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::labs(subtitle = NULL) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        )
}

plot_evpi <- function(fit, .color, .label) {
    suppressMessages({
        bayesDCA:::plot_evpi(
            fit,
            type = "useful", # same as type = "best" in this case, just using for plotting functionality
            colors = list(mced_test = .color),
            labels = list(mced_test = .label)
        ) +
            ggplot2::theme_bw(base_size = 24) +
            ggplot2::theme(legend.position = c(.75, .85)) +
            ggplot2::labs(subtitle = NULL, y = "Net true positives") +
            ggplot2::scale_y_continuous(
                breaks = scales::pretty_breaks(10, min.n = 6),
                labels = function(x) get_population_scaling_factor() * x
            ) +
            ggplot2::scale_x_continuous(
                labels = scales::label_percent(accuracy = 1),
                breaks = seq(0, 0.07, by = 0.01)
            ) +
            ggplot2::guides(color = "none")
    })
}

plot_trade_off_negative <- function(fit, .color, .label, max_abs_value = NULL, max_abs_value_tol = 10) {
    if (is.null(fit$nb_untreated)) {
        fit$nb_untreated <- compute_nb_untreated(fit)
    }

    if (is.data.frame(fit$nb_untreated)) {
        plot_df <- fit$nb_untreated
    } else {
        plot_df <- fit$nb_untreated$df
    }

    .colors <- get_color_values(.color = .color, .label = .label)
    label_function <- function(x) {
        x_lab <- paste0(x / 1000, "K")
        x_lab[x == 0L] <- "0"
        return(x_lab)
    }
    if (is.null(max_abs_value)) {
        max_abs_value <- pmin(
            max(abs(plot_df$upper_trade_off_negative)),
            max_abs_value_tol
        )
    }
    x_breaks <- seq(0, max_abs_value, by = 1)
    plot_df %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = thr,
                y = estimate_trade_off_negative,
                ymin = lower_trade_off_negative,
                ymax = upper_trade_off_negative
            )
        ) +
        ggplot2::geom_ribbon(alpha = 0.5, fill = .colors[.label]) +
        ggplot2::geom_line(linewidth = 1.5, ggplot2::aes(color = .label)) +
        ggplot2::coord_cartesian(ylim = c(-0.01, max_abs_value)) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "# tests",
            title = "MCED test trade-off",
            subtitle = "Number of tests for each additional net true negative"
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        ) +
        ggplot2::scale_y_continuous(
            breaks = x_breaks
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            name = NULL
        )
}

combine_trade_off_negative_plots <- function(
    pathways_dca_results,
    output_dir,
    .colors,
    legend_order,
    pathways = NULL,
    min_p_useful = 0.8) {
    trade_off_plots <- purrr::map(pathways_dca_results, ~ .x$trade_off_negative)
    if (!is.null(pathways)) {
        legend_order <- legend_order[legend_order %in% pathways]
    }
    trade_off_plots <- trade_off_plots[legend_order]
    df <- purrr::map_df(
        trade_off_plots,
        ~ .x$data %>%
            dplyr::filter(thr %in% c(0.02, 0.03, 0.04, 0.05), p_useful > min_p_useful) %>%
            dplyr::select(thr, contains("trade_off")),
        .id = "pathway"
    )
    df %>%
        dplyr::mutate(
            x = factor(
                paste0(thr * 100, "%")
            ),
            pathway = factor(
                as.character(pathway),
                levels = legend_order
            )
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = x,
                y = estimate_trade_off_negative,
                ymin = lower_trade_off_negative,
                ymax = upper_trade_off_negative
            )
        ) +
        ggplot2::geom_linerange(
            aes(color = pathway),
            position = ggplot2::position_dodge(width = .65),
            linewidth = 2
        ) +
        ggplot2::geom_point(
            ggplot2::aes(color = pathway),
            position = ggplot2::position_dodge(width = .65),
            size = 5
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = 2,
            linewidth = 1,
            alpha = 0.5
        ) +
        ggplot2::theme_bw(base_size = 24) +
        # ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::coord_cartesian(ylim = c(-0.01, NA)) +
        ggplot2::scale_y_continuous(
            breaks = scales::pretty_breaks(10)
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            breaks = legend_order
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "# tests per net true negative",
            title = "MCED test trade-off",
            subtitle = "Number of tests for each correctly avoided referral",
            color = NULL
        )
}

plot_trade_off_positive <- function(fit, .color, .label, max_abs_value = NULL, max_abs_value_tol = 200) {
    plot_df <- compute_trade_off_positive(fit = fit)

    .colors <- get_color_values(.color = .color, .label = .label)
    label_function <- function(x) {
        x_lab <- paste0(x / 1000, "K")
        x_lab[x == 0L] <- "0"
        return(x_lab)
    }
    if (is.null(max_abs_value)) {
        max_abs_value <- min(
            max(abs(plot_df$upper_trade_off_positive)),
            max_abs_value_tol
        )
    }
    x_breaks <- seq(0, max_abs_value, by = 1)
    plot_df %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = thr,
                y = estimate_trade_off_positive,
                ymin = lower_trade_off_positive,
                ymax = upper_trade_off_positive
            )
        ) +
        ggplot2::geom_ribbon(alpha = 0.5, fill = .colors[.label]) +
        ggplot2::geom_line(linewidth = 1.5, ggplot2::aes(color = .label)) +
        ggplot2::coord_cartesian(ylim = c(-0.01, max_abs_value)) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "# tests per net true positive",
            title = "MCED test trade-off",
            subtitle = "Number of tests to refer one true cancer"
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        ) +
        ggplot2::scale_y_continuous(
            breaks = scales::pretty_breaks(10)
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            name = NULL
        )
}

combine_trade_off_positive_plots <- function(
    pathways_dca_results,
    output_dir,
    .colors,
    legend_order,
    pathways = NULL,
    min_p_useful = 0.8) {
    trade_off_plots <- purrr::map(pathways_dca_results, ~ .x$trade_off_positive)
    if (!is.null(pathways)) {
        legend_order <- legend_order[legend_order %in% pathways]
    }
    trade_off_plots <- trade_off_plots[legend_order]
    df <- purrr::map_df(
        trade_off_plots,
        ~ .x$data %>%
            dplyr::filter(thr %in% c(0.02, 0.03, 0.04, 0.05)) %>%
            dplyr::select(thr, contains("trade_off")),
        .id = "pathway"
    )
    p_useful_data <- purrr::map_df(pathways_dca_results, ~ .x$p_useful$data, .id = "pathway") %>%
        dplyr::rename(thr := threshold, p_useful := prob)
    dplyr::left_join(
        df,
        p_useful_data,
        by = c("pathway", "thr")
    ) %>%
        dplyr::filter(p_useful > min_p_useful) %>%
        dplyr::mutate(
            x = factor(
                paste0(thr * 100, "%")
            ),
            pathway = factor(
                as.character(pathway),
                levels = legend_order
            )
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = x,
                y = estimate_trade_off_positive,
                ymin = lower_trade_off_positive,
                ymax = upper_trade_off_positive
            )
        ) +
        ggplot2::geom_linerange(
            aes(color = pathway),
            position = ggplot2::position_dodge(width = .65),
            linewidth = 2
        ) +
        ggplot2::geom_point(
            ggplot2::aes(color = pathway),
            position = ggplot2::position_dodge(width = .65),
            size = 5
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = 2,
            linewidth = 1,
            alpha = 0.5
        ) +
        ggplot2::theme_bw(base_size = 24) +
        # ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::coord_cartesian(ylim = c(-0.01, NA)) +
        ggplot2::scale_y_continuous(
            breaks = scales::pretty_breaks(10)
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            breaks = legend_order
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "# tests per net true positive",
            title = "MCED test trade-off",
            subtitle = "Number of tests to refer one true cancer",
            color = NULL
        )
}
