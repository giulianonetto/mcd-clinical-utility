library(targets)
library(here)

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
      symplify_pathways_data = here("data/symplify-pathways.tsv"),
      output_dir = here::here("output/"),
      n_draws = 5e4
    )
  ),
  tar_target(
    name = optimizing_mced_test,
    command = run_optimizing_mced_test(
      symplify_pathways_data = here("data/symplify-pathways.tsv"),
      output_dir = here::here("output/")
    )
  ),
  tar_target(
    name = estimating_optimal_cutoff,
    command = run_estimating_optimal_cutoff(
      output_dir = here::here("output/")
    )
  )
)
