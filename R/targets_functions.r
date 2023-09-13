run_symplify_pathways_dca <- function(symplify_pathways_data, output_dir, n_draws = 4000) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    extracted_data <- readr::read_tsv(symplify_pathways_data, show_col_types = FALSE)
    extracted_data <- map_df(
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
    pathways_dca_results <- vector("list", length = length(pathways))
    names(pathways_dca_results) <- pathways
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
            thresholds = seq(0, ceiling(extracted_data_pathway$p * 1.2) / 100, length = 100),
            n_draws = n_draws
        )

        .label <- stringr::str_glue("{.pathway} pathway")
        .color <- .colors[[.pathway]]
        # TODO: prevalence vertical lines
        dca_treated <- plot_net_benefit_treated(fit = .fit, .color = .color, .label = .label)
        dca_untreated <- plot_net_benefit_untreated(fit = .fit, .color = .color, .label = .label)
        p_useful <- plot_prob_useful(fit = .fit, .color = .color, .label = .label)
        evpi <- plot_evpi(fit = .fit, .color = .color, .label = .label)

        pathways_dca_results[[.pathway]] <- list(
            dca_treated = dca_treated,
            dca_untreated = dca_untreated,
            p_useful = p_useful,
            evpi = evpi
        )
    }

    logger::log_info("Creating final pathways figures")
    create_pathways_figures(
        pathways_dca_results = pathways_dca_results
    )
    logger::log_info("Figures successfully created!")
    return(pathways_dca_results[[1]][[1]])
}
