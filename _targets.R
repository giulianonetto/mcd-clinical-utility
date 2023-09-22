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
      output_dir = here::here("output/testing"),
      n_draws = 2e4
    )
  )
)
