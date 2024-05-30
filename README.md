# Clinical utility of multi-cancer detection in symptomatic patients: a decision-making perspective

This repository assesses the clinical utility of the Galleri test for multi-cancer detection in symptomatic patients referred for cancer investigation, using data from the [SYMPLIFY cohort study](https://doi.org/10.1016/S1470-2045(23)00277-2).

To reproduce all results in the `output` folder, please run:

```
docker build -t cruzandkorthauer2024img .
docker run -it -v ${PWD}:/home/rstudio/ \
    cruzandkorthauer2024img \
    R -e "targets::tar_destroy(ask = FALSE); targets::tar_make()"
```
