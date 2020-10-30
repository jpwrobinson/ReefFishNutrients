# ReefFishNutrients

Code and data for Robinson et al.

```models.../``` directories contain R code (.Rmd) used for Bayesian model analyses of micronutrient concentrations in reef fish, for the trait model (```models_inter_species```), trait-habitat model (```models_intra_habitat```), and reef benthic gradient model (```models_uvc_benthic```). 

All models run in [Stan](https://mc-stan.org/) using [rethinking](https://github.com/rmcelreath) in [R](https://cran.r-project.org/).

Data provided are raw nutrient measurements (```nutrient_data.csv```), species-level posterior predictions (```nutrient_posteriors.csv```), and UVC datasets for fish (`uvc_fish.csv`) and benthos (`uvc_benthic.csv`).

