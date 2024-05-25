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
            "Analyzing the {.pathway} pathway... "
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

run_optimizing_mced_test <- function(symplify_pathways_data, output_dir, l = 201) {
    .colors <- c(
        "Overall" = "#1B9E77",
        "Gynaecology" = "#8d36f0",
        "Upper GI" = "#D95F02",
        "Lower GI" = "#0e4674",
        "Lung" = "#1c87df",
        "RDC" = "#E7298A"
    )

    get_ntn <- function(.se, .sp, .p, .t) .sp * (1 - .p) - .p * (1 - .se) * ((1 - .t) / .t)
    get_refer_none <- function(.p, .t) 1 * (1 - .p) - .p * (1 - 0) * ((1 - .t) / .t)
    get_best_competitor <- function(.p, .t) {
        x <- get_refer_none(.p = .p, .t = .t)
        x[x < 0] <- 0
        return(x)
    }

    compute_clinical_utility <- function(l = 201,
                                         potential_sens_spec = seq(0, 1, length = l),
                                         potential_prev = c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.3),
                                         potential_thresholds = c(0.02, 0.03, 0.05)) {
        expand.grid(
            se = potential_sens_spec,
            sp = potential_sens_spec,
            p = potential_prev,
            thr = potential_thresholds
        ) %>%
            mutate(
                ntn = get_ntn(se, sp, p, thr) * 1e5,
                ntn_default = get_best_competitor(p, thr) * 1e5,
                deltaNB = ntn - ntn_default,
                useful = factor(ifelse(deltaNB > 0, "Useful", "Not useful"), levels = c("Useful", "Not useful")),
                tradeoff = 1 / deltaNB,
                p = fct_relabel(ordered(round(p * 100, 1)), \(x) paste0("prev. ", x, "%")),
                thr = fct_relabel(ordered(round(thr * 100)), \(x) paste0(x, "%\nthreshold"))
            )
    }

    plot_nb_gain_regions <- function(df) {
        limits <- c(-0.01, NA)
        color_breaks <- c(10, 25, 50, 75, 90)
        df %>%
            mutate(
                tradeoff2 = case_when(
                    tradeoff < 0 ~ NA,
                    tradeoff == 0 ~ "0",
                    tradeoff <= 10 ~ "10",
                    tradeoff <= 50 ~ "50",
                    tradeoff <= 100 ~ "100",
                    tradeoff <= 500 ~ "500",
                    tradeoff <= 1000 ~ "1000",
                    tradeoff <= 2000 ~ "5000",
                    TRUE ~ ">5000",
                ),
                tradeoff2 = factor(
                    tradeoff,
                    levels = c("0", "10", "50", "100", "500", "1000", "5000", ">5000"),
                    ordered = TRUE
                )
            ) %>%
            ggplot(aes(sp, se)) +
            geom_raster(
                data = . %>% dplyr::filter(useful == "Useful"),
                ggplot2::aes(fill = deltaNB)
            ) +
            scale_fill_viridis_b(
                name = "Unnecessary referrals avoided\n(Net TN per 100K)",
                limits = c(0, 100) * 1e3,
                breaks = color_breaks * 1e3,
                labels = paste0(color_breaks, "K"),
                direction = -1
            ) +
            ggnewscale::new_scale_fill() +
            # geom_raster(
            #     data = . %>% dplyr::filter(useful == "Not useful"),
            #     ggplot2::aes(fill = "Not useful"),
            #     show.legend = FALSE
            # ) +
            # scale_fill_manual(values = "#c0c0c0") +
            scale_fill_manual(values = "#f0f0f0") +
            facet_wrap(~thr) +
            scale_color_manual(values = .colors) +
            scale_y_continuous(labels = \(x) scales::percent(x, suffix = NULL), breaks = scales::pretty_breaks()) +
            scale_x_continuous(labels = \(x) scales::percent(x, suffix = NULL), breaks = scales::pretty_breaks()) +
            theme_minimal(base_size = 20) +
            theme(
                axis.text.y = element_text(size = 12),
                legend.position = "top",
                legend.title = element_text(hjust = 0.5),
                axis.text.x = element_text(size = 12),
                panel.grid.major = element_line(linewidth = 0.25),
                panel.grid.minor = element_line(linewidth = 0.1),
                legend.key.width = unit(2, "cm")
            ) +
            labs(
                fill = NULL, color = NULL,
                x = "Specificity (%)", y = "Sensitivity (%)"
            ) +
            # geom_vline(xintercept = c(0.8, 0.9), linetype = 2, linewidth = 1, alpha = 0.7) +
            facet_grid(cols = vars(p), rows = vars(thr))
    }

    df_clinical_utility <- compute_clinical_utility(l = l)
    extracted_data <- readr::read_tsv(symplify_pathways_data, show_col_types = FALSE) %>%
        dplyr::mutate(
            p = ordered(paste0("prev. ", round(100 * D / N), "%"), levels = levels(df_clinical_utility$p)),
            pathway = factor(pathway, levels = names(.colors))
        )

    extracted_data <- purrr::map_df(
        c(2, 3, 5), ~ {
            extracted_data %>%
                dplyr::mutate(
                    thr = paste0(.x, "%\nthreshold")
                )
        }
    )

    extracted_data$thr <- ordered(extracted_data$thr, levels = levels(df_clinical_utility$thr))

    p <- df_clinical_utility %>%
        plot_nb_gain_regions() +
        ggplot2::geom_point(
            data = extracted_data,
            ggplot2::aes(x = sp, y = se, color = pathway),
            inherit.aes = FALSE,
            size = 3.5
        ) +
        ggplot2::geom_point(
            data = extracted_data,
            ggplot2::aes(x = sp, y = se),
            inherit.aes = FALSE,
            size = 3.5, pch = 21, color = "gray20"
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4.5)))

    ggplot2::ggsave(
        here::here(file.path(output_dir, "fig02.png")),
        p,
        width = 16, height = 8,
        bg = "white",
        dpi = 600
    )

    df_clinical_utility %>%
        dplyr::mutate(thr = stringr::str_replace_all(thr, "\n", " ")) %>%
        dplyr::filter(deltaNB > 0) %>%
        readr::write_tsv(
            here::here(file.path(output_dir, "clinical-utility-regions.tsv"))
        )
}

run_estimating_optimal_cutoff <- function(output_dir, seed = 1234567) {
    withr::with_seed(
        seed = seed,
        code = {
            n <- 100
            preds <- seq(-4, 4, length = n)
            noise <- rpois(n, 3e2) * sample(c(-1, 1), n, replace = TRUE)
            ntn <- (
                8e4 + 1e3 * preds - 7e3 * preds^2 - 1e2 * preds^3 + noise
            )
            error_margin <- 1.5e4
            p <- data.frame(
                preds = preds,
                ntn = ntn
            ) %>%
                ggplot2::ggplot(ggplot2::aes(preds, ntn)) +
                ggplot2::geom_ribbon(
                    ggplot2::aes(ymin = ntn - error_margin, ymax = ntn + error_margin),
                    fill = "gray80", alpha = 0.8
                ) +
                ggplot2::geom_line(linewidth = 1.75) +
                ggplot2::geom_segment(
                    x = preds[which.max(ntn)],
                    xend = preds[which.max(ntn)],
                    y = 0,
                    yend = max(ntn) + 1,
                    linetype = 2,
                    color = "red",
                    linewidth = 1
                ) +
                ggplot2::geom_point(
                    ggplot2::aes(x = preds[which.max(ntn)], y = max(ntn)),
                    color = "red",
                    size = 6
                ) +
                ggplot2::geom_segment(
                    x = min(preds),
                    xend = preds[which.max(ntn)],
                    y = max(ntn) + 1,
                    yend = max(ntn) + 1,
                    linetype = 2,
                    color = "red",
                    linewidth = 1
                ) +
                ggplot2::scale_y_continuous(labels = scales::comma) +
                ggplot2::geom_hline(yintercept = 0, linetype = 2) +
                ggplot2::theme_minimal(base_size = 20) +
                ggplot2::labs(
                    x = "Predicted score",
                    y = "Net true negatives per 100K",
                    title = "Hypothetical example of optimal cutoff estimation",
                    subtitle = "Unnecessary referrals avoided by predicted score cutoff",
                ) +
                ggplot2::coord_cartesian(ylim = c(1000, 1e5))

            ggplot2::ggsave(
                here::here(file.path(output_dir, "fig03.png")),
                p,
                width = 10, height = 6,
                bg = "white",
                dpi = 600
            )
        }
    )
}
