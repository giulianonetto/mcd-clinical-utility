run_symplify_pathways_dca <- function(symplify_pathways_data, output_dir) {
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

    pathways <- unique(extracted_data$pathway)
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
        fit <- bayesDCA::dca(
            dca_data,
            thresholds = seq(0, ceiling(extracted_data_pathway$p * 1.2) / 100, length = 100),
            n_draws = 2e4
        )
        .colors <- list(
            mced_test = dplyr::case_match(
                .pathway,
                "Overall" ~ "#1B9E77",
                "Gynaecology" ~ "#8d36f0",
                "Lower GI" ~ "#0e4674",
                "Lung" ~ "#1c87df",
                "RDC" ~ "#E7298A",
                "Upper GI" ~ "#D95F02"
            )
        )

        dca_untreated <- plot_net_benefit_opt_out(fit, .title = NULL, .colors = .colors$mced_test)

        .labels <- list(mced_test = stringr::str_glue("Gallery test ({.pathway} pathway)"))
        p1 <- bayesDCA::compare_dca(
            fit,
            .evpi = TRUE,
            colors = .colors,
            labels = .labels
        ) # TODO: add prevalence line
        .pathway <- stringr::str_replace_all(.pathway, " ", "-")
        file_name_p1 <- stringr::str_glue(
            "pathway-dca-{.pathway}.png"
        )
        ggplot2::ggsave(
            file.path(output_dir, file_name_p1),
            width = 12,
            height = 7.5
        )
    }
}
