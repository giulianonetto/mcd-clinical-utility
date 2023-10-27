run_symplify_pathways_dca <- function(symplify_pathways_data, output_dir, n_draws = 4000) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    extracted_data <- readr::read_tsv(symplify_pathways_data, show_col_types = FALSE)
    extracted_data <- purrr::map_df(
        seq_len(nrow(extracted_data)),
        function(i) {
            x <- extracted_data[i, ]
            process_extracted_data(
                n = x$N, d = x$D,
                se = x$se, sp = x$sp,
                name = x$pathway
            )
        }
    ) %>%
        dplyr::rename(pathway := name)

    .colors <- c(
        "Overall" = "#1B9E77",
        "Gynaecology" = "#8d36f0",
        "Lower GI" = "#0e4674",
        "Lung" = "#1c87df",
        "RDC" = "#E7298A",
        "Upper GI" = "#D95F02"
    )
    pathways <- names(.colors)
    # decision_thresholds <- seq(0, ceiling(extracted_data_pathway$p * 1.2) / 100, length = 100)
    decision_thresholds <- seq(0, 0.075, by = 0.001)
    pathways_dca_results <-
        dca_treated_tables <-
        dca_untreated_tables <-
        p_useful_tables <-
        trade_off_positive_tables <-
        trade_off_negative_tables <-
        setNames(vector("list", length = length(pathways)), pathways)
    for (.pathway in pathways) {
        msg <- stringr::str_glue(
            "Analyzing the {.pathway} pathway..."
        )
        logger::log_info(cli::col_br_blue(msg))
        extracted_data_pathway <- extracted_data[extracted_data$pathway == .pathway, ]
        dca_data <- generate_data_from_counts(
            n = extracted_data_pathway$n,
            d = extracted_data_pathway$d,
            tp = extracted_data_pathway$tp,
            tn = extracted_data_pathway$tn,
            fp = extracted_data_pathway$fp,
            fn = extracted_data_pathway$fn
        )
        .fit <- bayesDCA::dca(
            dca_data,
            thresholds = decision_thresholds,
            n_draws = n_draws
        )
        .fit$nb_untreated <- compute_nb_untreated(.fit)

        .label <- ifelse(
            .pathway == "Overall",
            stringr::str_glue("MCED test ({.pathway})"),
            stringr::str_glue("MCED test ({.pathway} pathway)")
        )
        .color <- .colors[[.pathway]]
        dca_treated <- plot_net_benefit_treated(fit = .fit, .color = .color, .label = .label)
        dca_untreated <- plot_net_benefit_untreated(fit = .fit, .color = .color, .label = .label)
        delta_treated <- plot_delta_treated(fit = .fit, .color = .color, .label = .label)
        p_useful <- plot_prob_useful(fit = .fit, .color = .color, .label = .label)
        evpi <- plot_evpi(fit = .fit, .color = .color, .label = .label)
        trade_off_negative <- plot_trade_off_negative(fit = .fit, .color = .color, .label = .label)
        trade_off_positive <- plot_trade_off_positive(fit = .fit, .color = .color, .label = .label)

        pathways_dca_results[[.pathway]] <- list(
            dca_treated = dca_treated,
            dca_untreated = dca_untreated,
            p_useful = p_useful,
            delta_treated = delta_treated,
            evpi = evpi,
            trade_off_negative = trade_off_negative,
            trade_off_positive = trade_off_positive
        )
        # add vertical line at t=3%
        pathways_dca_results[[.pathway]] <- map(
            pathways_dca_results[[.pathway]],
            ~ {
                .x$layers <- c(
                    ggplot2::geom_vline(xintercept = 0.03, linetype = 3, linewidth = 1.5, alpha = 1),
                    .x$layers
                )
                .x + ggplot2::theme(panel.grid = ggplot2::element_blank())
            }
        )
        dca_treated_tables[[.pathway]] <- .fit$summary$net_benefit
        dca_untreated_tables[[.pathway]] <- .fit$nb_untreated
        p_useful_tables[[.pathway]] <- p_useful$data
        trade_off_negative_tables[[.pathway]] <- trade_off_negative$data
        trade_off_positive_tables[[.pathway]] <- trade_off_positive$data
    }

    logger::log_info("Saving all results tables")
    dca_treated_all <- dplyr::bind_rows(dca_treated_tables, .id = "pathway")
    dca_untreated_all <- dplyr::bind_rows(dca_untreated_tables, .id = "pathway")
    p_useful_all <- dplyr::bind_rows(p_useful_tables, .id = "pathway")
    trade_off_negative_all <- dplyr::bind_rows(trade_off_negative_tables, .id = "pathway")
    trade_off_positive_all <- dplyr::bind_rows(trade_off_positive_tables, .id = "pathway")
    readr::write_tsv(
        dca_untreated_all,
        here::here(stringr::str_glue("{output_dir}/dca_untreated.tsv"))
    )
    readr::write_tsv(
        dca_treated_all,
        here::here(stringr::str_glue("{output_dir}/dca_treated.tsv"))
    )
    readr::write_tsv(
        p_useful_all,
        here::here(stringr::str_glue("{output_dir}/p_useful.tsv"))
    )
    readr::write_tsv(
        trade_off_negative_all,
        here::here(stringr::str_glue("{output_dir}/trade_off_negative.tsv"))
    )
    readr::write_tsv(
        trade_off_positive_all,
        here::here(stringr::str_glue("{output_dir}/trade_off_positive.tsv"))
    )
    logger::log_info("Done!")

    logger::log_info("Creating final pathways figures")
    create_final_figures(
        pathways_dca_results = pathways_dca_results,
        .colors = .colors,
        output_dir = output_dir
    )
    logger::log_info("Figures successfully created!")
    return("done")
}
