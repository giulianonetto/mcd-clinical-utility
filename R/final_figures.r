create_final_figures <- function(
    pathways_dca_results,
    .colors,
    output_dir = "output/testing") {
    # figures config
    figure_width <- 24
    figure_height <- figure_width * 1.4
    column_titles_font_size <- figure_width * 2.25
    column_subtitles_font_size <- figure_width * 1.3
    row_titles_font_size <- figure_width * 1.5
    x_axis_font_size <- figure_width * 0.85
    y_axis_font_size <- figure_width * 0.7
    legend_font_size <- figure_width * 1.25
    lines_order <- c(
        "RDC",
        "Gynaecology",
        "Upper GI",
        "Lower GI",
        "Overall",
        "Lung"
    )
    legend_order <- c(
        "Overall",
        "Gynaecology",
        "Upper GI",
        "Lower GI",
        "Lung",
        "RDC"
    )
    # "ground truth" Net TN for testing (using point estimates only)
    ground_truth_0.03 <- purrr::map_df(
        list(
            Overall = c(se = 0.663, sp = 0.984, p = 0.067),
            Gynaecology = c(se = 0.491, sp = 0.983, p = 0.037),
            "Lower GI" = c(se = 0.685, sp = 0.991, p = 0.065),
            "Upper GI" = c(se = 0.804, sp = 0.981, p = 0.045),
            "Lung" = c(se = 0.708, sp = 0.962, p = 0.298),
            "RDC" = c(se = 0.556, sp = 0.982, p = 0.073)
        ), ~ {
            w <- (1 - 0.03) / 0.03
            ntn <- .x["sp"] * (1 - .x["p"]) - w * (1 - .x["se"]) * .x["p"]
            c(ntn = unname(ntn))
        },
        .id = "pathway"
    ) %>%
        dplyr::mutate(
            pathway = factor(
                as.character(pathway),
                levels = legend_order
            )
        )

    # Figure 1
    combined_plots <- purrr::imap(
        pathways_dca_results,
        ~ {
            p1 <- .x$dca_treated +
                ggplot2::labs(y = .y, subtitle = NULL) +
                ggplot2::theme(
                    axis.title.y = ggplot2::element_text(
                        size = row_titles_font_size,
                        face = "bold",
                        vjust = 1.5
                    )
                )
            p2 <- .x$p_useful + ggplot2::labs(y = NULL, subtitle = NULL)
            p3 <- .x$dca_untreated +
                ggplot2::labs(y = NULL, subtitle = NULL) +
                ggplot2::guides(color = "none")
            (p1 | p2 | p3) +
                patchwork::plot_layout(
                    guides = "collect",
                    widths = c(3.5, 3, 3.5)
                ) &
                ggplot2::theme(
                    legend.position = "top",
                    legend.text = ggplot2::element_text(size = legend_font_size),
                    legend.justification = "left",
                    axis.text.x = ggplot2::element_text(size = x_axis_font_size),
                    axis.text.y = ggplot2::element_text(size = y_axis_font_size)
                )
        }
    )

    column_title_config <- ggplot2::element_text(
        hjust = 0.5,
        vjust = -3.5,
        size = column_titles_font_size,
        face = "bold"
    )
    combined_plots$Overall[[1]] <- combined_plots$Overall[[1]] +
        ggplot2::labs(
            title = "Net benefit"
        ) +
        ggplot2::theme(
            plot.title = column_title_config
        )
    combined_plots$Overall[[2]] <- combined_plots$Overall[[2]] +
        ggplot2::labs(
            title = "P(useful)"
        ) +
        ggplot2::theme(
            plot.title = column_title_config
        )
    combined_plots$Overall[[3]] <- combined_plots$Overall[[3]] +
        ggplot2::labs(
            title = "Unnecessary referrals\navoided",
            subtitle = stringr::str_glue(
                "Net true negatives\nper {get_population_scaling_factor(as_string = TRUE)} patients"
            )
        ) +
        ggplot2::theme(
            plot.title = column_title_config,
            plot.subtitle = ggplot2::element_text(hjust = 0.5, vjust = -6, size = column_subtitles_font_size)
        )

    fig01 <- (
        combined_plots$Overall /
            combined_plots$Gynaecology /
            combined_plots$`Upper GI` /
            combined_plots$`Lower GI` /
            combined_plots$Lung /
            combined_plots$RDC
    ) & theme(plot.margin = margin(t = 0.5, r = 2, b = 0.5, l = 1, unit = "cm"))
    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/fig01.png")),
        fig01,
        width = figure_width,
        height = figure_height,
        dpi = 600
    )

    evpi_data <- purrr::map_df(
        pathways_dca_results,
        ~ .x$evpi$data,
        .id = "pathway"
    ) %>%
        dplyr::mutate(
            pathway = factor(
                as.character(pathway),
                levels = lines_order
            )
        )

    column_title_config <- ggplot2::element_text(
        hjust = 0,
        size = figure_width * 1.35,
        face = "bold"
    )
    supp_fig01 <- evpi_data %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = threshold,
                y = .evpi,
                color = pathway
            )
        ) +
        ggplot2::geom_vline(xintercept = 0.03, linetype = 3, linewidth = 1.5, alpha = 1) +
        ggplot2::geom_line(linewidth = 2) +
        ggplot2::theme_bw(base_size = 24) +
        ggplot2::theme(
            plot.title = column_title_config,
            plot.subtitle = ggplot2::element_text(hjust = 0, size = column_title_config$size * 0.8),
            legend.position = c(.8, .75),
            panel.grid = ggplot2::element_blank(),
            legend.key.height = ggplot2::unit(0.75, "cm")
        ) +
        ggplot2::labs(
            title = "Expected Value of Perfect Information",
            subtitle = stringr::str_glue(
                "per {get_population_scaling_factor(as_string = TRUE)} patients"
            ),
            y = "Net true positives",
            x = "Decision threshold",
            color = NULL
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            breaks = legend_order
        ) +
        ggplot2::scale_y_continuous(
            breaks = scales::pretty_breaks(10, min.n = 6),
            labels = function(x) get_population_scaling_factor() * x
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::label_percent(accuracy = 1),
            breaks = seq(0, 0.07, by = 0.01)
        )

    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/evpi.png")),
        supp_fig01,
        width = 10.5,
        height = 7.5
    )

    fig02 <- combine_trade_off_negative_plots(
        pathways_dca_results = pathways_dca_results,
        output_dir = here::here(output_dir),
        legend_order = legend_order,
        pathways = legend_order[legend_order != "Overall"],
        .colors = .colors
    )

    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/tradeoff.png")),
        fig02,
        width = 10.5,
        height = 7.5
    )

    tradeoff_positive <- combine_trade_off_positive_plots(
        pathways_dca_results = pathways_dca_results,
        output_dir = here::here(output_dir),
        legend_order = legend_order,
        pathways = legend_order[legend_order != "Overall"],
        .colors = .colors
    )

    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/tradeoff_positive.png")),
        tradeoff_positive,
        width = 10.5,
        height = 7.5
    )

    a <- fig02 +
        ggplot2::ggtitle(NULL) +
        ggplot2::theme(
            legend.position = c(0.16, 0.74),
            legend.key.height = ggplot2::unit(0.8, "cm")
        )
    b <- tradeoff_positive + ggplot2::ggtitle(NULL) + ggplot2::guides(color = "none")
    both_tradeoffs <- (a | b) + patchwork::plot_annotation(tag_levels = c("A", "B"))
    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/both_tradeoffs.png")),
        both_tradeoffs,
        width = 20,
        height = 7
    )
}
