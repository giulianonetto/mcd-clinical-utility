create_pathways_figures <- function(pathways_dca_results) {
    figure_width <- 24
    figure_height <- figure_width * 1.3
    column_titles_font_size <- figure_width * 1.3
    row_titles_font_size <- figure_width
    x_axis_font_size <- figure_width * 0.85
    y_axis_font_size <- figure_width * 0.7
    legend_font_size <- figure_width * 0.9
    net_neg_label_function <- function(x) {
        x_lab <- paste0(x / 1000, "K")
        x_lab[x == 0L] <- "0"
        return(x_lab)
    }
    combined_plots <- purrr::imap(
        pathways_dca_results,
        ~ {
            p1 <- .x$dca_treated +
                ggplot2::labs(y = .y, subtitle = NULL) +
                ggplot2::theme(
                    axis.title.y = ggplot2::element_text(
                        size = row_titles_font_size,
                        face = "bold"
                    )
                )
            p2 <- .x$p_useful + ggplot2::labs(y = NULL, subtitle = NULL)
            p3 <- .x$dca_untreated +
                ggplot2::scale_y_continuous(
                    labels = net_neg_label_function
                ) +
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
                    axis.text.x = ggplot2::element_text(size = x_axis_font_size),
                    axis.text.y = ggplot2::element_text(size = y_axis_font_size)
                )
        }
    )

    column_title_config <- ggplot2::element_text(
        hjust = 0.5,
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
                "per {get_population_scaling_factor(as_string = TRUE)} patients"
            )
        ) +
        ggplot2::theme(
            plot.title = column_title_config,
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = column_titles_font_size * 0.85)
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


    evpi_plots <- purrr::imap(
        pathways_dca_results,
        ~ {
            .x$evpi +
                ggplot2::labs(y = NULL, subtitle = NULL, title = .y) +
                ggplot2::theme(
                    plot.title = column_title_config
                )
        }
    )

    column_title_config$hjust <- 0

    fig02 <- (
        evpi_plots$Overall |
            evpi_plots$Gynaecology |
            evpi_plots$`Upper GI`
    ) /
        (
            evpi_plots$`Lower GI` |
                evpi_plots$Lung |
                evpi_plots$RDC
        ) +
        patchwork::plot_annotation(
            title = "Expected Value of Perfect Information",
            subtitle = stringr::str_glue(
                "per {get_population_scaling_factor(as_string = TRUE)} patients"
            ),
            theme = ggplot2::theme(
                plot.title = column_title_config,
                plot.subtitle = ggplot2::element_text(hjust = 0, size = column_titles_font_size * 0.85)
            )
        )
    ggplot2::ggsave(
        here::here(stringr::str_glue("{output_dir}/fig02.png")),
        fig02,
        width = figure_width * 0.8,
        height = figure_height * 0.4
    )
}
