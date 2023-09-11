library(targets)

tar_option_set(
  packages = c(
    "tidyverse",
    "patchwork"
  )
)
# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Target list:
list(
  tar_target(
    name = pathways_dca_results,
    command = run_symplify_pathways_dca(
      symplify_pathways_data = here::here("data/symplify-pathways.tsv")
    )
  ),
  tar_target(
    name = pathways_figures,
    command = create_pathways_figures(
      pathways_dca_results = pathways_dca_results,
      output_dir = here::here("output/pathways")
    )
  )
)
