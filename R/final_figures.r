create_pathways_figures <- function(
    pathways_dca_results,
    .colors,
    output_dir = "output/testing") {
    figure_width <- 24
    figure_height <- figure_width * 1.4
    column_titles_font_size <- figure_width * 2.25
    column_subtitles_font_size <- figure_width * 1.3
    row_titles_font_size <- figure_width * 1.5
    x_axis_font_size <- figure_width * 0.85
    y_axis_font_size <- figure_width * 0.7
    legend_font_size <- figure_width * 1.25
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
            p4 <- .x$evpi +
                ggplot2::labs(y = NULL, subtitle = NULL) +
                ggplot2::guides(color = "none") +
                ggplot2::theme(
                    plot.margin = margin(t = 0.3, r = 4, b = 0.3, l = 2, "cm")
                )
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
            title = "Net true negatives",
            subtitle = stringr::str_glue(
                "Unnecessary referrals avoided\nper {get_population_scaling_factor(as_string = TRUE)} patients"
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
    )
    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/fig01.png")),
        fig01,
        width = figure_width,
        height = figure_height
    )


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
    fig02 <- evpi_data %>%
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
        here::here(stringr::str_glue("{output_dir}/fig02.png")),
        fig02,
        width = 10.5,
        height = 7.5
    )
}

create_trade_off_figure <- function(pathways_dca_results, output_dir) {
    trade_off_plots <- purrr::map(pathways_dca_results, ~ .x$trade_off)
    df <- purrr::map_df(
        trade_off_plots,
        ~ .x$data %>%
            dplyr::filter(thr %in% c(0.02, 0.03)) %>%
            dplyr::select(thr, contains("trade_off")),
        .id = "pathway"
    )
    df %>%
        dplyr::mutate(
            x = factor(
                paste0(thr * 100, "%")
            )
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = x,
                y = estimate_trade_off,
                ymin = lower_trade_off,
                ymax = upper_trade_off
            )
        ) +
        ggplot2::geom_pointrange(
            aes(color = pathway),
            position = ggplot2::position_dodge(width = .4)
        ) +
        ggplot2::theme_bw(base_size = 24) +
        # ggplot2::theme(legend.position = c(.15, .85)) +
        ggplot2::coord_cartesian(ylim = c(-0.5, NA)) +
        ggplot2::scale_color_discrete(name = NULL) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "# tests",
            title = "MCED test trade-off",
            subtitle = "Number of tests for each additional net true negative"
        )
}
