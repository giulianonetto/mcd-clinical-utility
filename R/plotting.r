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
plot_net_benefit_treated <- function(fit, .color, .label) {
    color_values <- get_color_values(.color = .color, .label = .label)
    label_values <- get_labels_values(.label = .label)
    p <- bayesDCA:::plot.BayesDCA(fit) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(legend.position = c(.85, .85)) +
        ggplot2::scale_color_manual(values = color_values, labels = label_values) +
        ggplot2::scale_fill_manual(values = color_values, labels = label_values) +
        ggplot2::coord_cartesian(ylim = c(-0.001, NA))
    p$layers[[1]]$aes_params <- list(alpha = 0.5)
    return(p)
}

plot_net_benefit_untreated <- function(fit, .color, .label) {
    plot_df <- compute_nb_untreated(fit)
    .colors <- get_color_values(.color = .color, .label = .label)
    label_function <- function(x) {
        x_lab <- paste0(x / 1000, "K")
        x_lab[x == 0L] <- "0"
        return(x_lab)
    }
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
            alpha = 1,
            linewidth = 1.5,
            ggplot2::aes(color = "Refer all", yintercept = 0)
        ) +
        ggplot2::coord_cartesian(ylim = c(-0.01, NA)) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
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
        ggplot2::guides(color = "none")
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
            ggplot2::guides(color = "none")
    })
}
