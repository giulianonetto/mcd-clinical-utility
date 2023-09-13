plot_net_benefit_treated <- function(fit, .color, .label) {
    bayesDCA:::plot.BayesDCA(
        fit,
        colors = list(mced_test = .color),
        labels = list(mced_test = .label)
    ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(legend.position = c(.85, .85))
}

plot_net_benefit_untreated <- function(fit, .color, .label) {
    plot_df <- compute_nb_untreated(fit)

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

plot_prob_useful <- function(fit, .color, .label) {
    bayesDCA:::plot_superiority_prob(
        fit,
        type = "useful",
        colors = list(mced_test = .color),
        labels = list(mced_test = .label)
    ) +
        ggplot2::theme_bw(base_size = 18) +
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


create_pathways_figures <- function(pathways_dca_results) {
    combined_plots <- purrr::imap(
        pathways_dca_results,
        ~ {
            p1 <- .x$dca_treated +
                ggplot2::labs(y = .y, subtitle = NULL) +
                ggplot2::theme(axis.title.y = ggplot2::element_text(size = 25))
            p2 <- .x$p_useful + ggplot2::labs(y = NULL, subtitle = NULL)
            p3 <- .x$dca_untreated + ggplot2::labs(y = NULL, subtitle = NULL)
            p1 | p2 | p3
        }
    )

    combined_plots$Overall[[1]] <- combined_plots$Overall[[1]] +
        ggplot2::labs(
            title = "Net benefit"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 35)
        )
    combined_plots$Overall[[2]] <- combined_plots$Overall[[2]] +
        ggplot2::labs(
            title = "P(useful)"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 35)
        )
    combined_plots$Overall[[3]] <- combined_plots$Overall[[3]] +
        ggplot2::labs(
            title = "Net true negatives",
            subtitle = "per 1000 patients"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 35),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 30)
        )



    fig01 <- (
        combined_plots$Overall /
            combined_plots$Gynaecology /
            combined_plots$`Upper GI` /
            combined_plots$`Lower GI` /
            combined_plots$Lung /
            combined_plots$RDC
    )
    ggplot2::ggsave(
        here::here("output/testing/fig01.png"),
        fig01,
        width = 20, height = 20
    )
}
