## Inference

This directory contains all of the code necessary to estimate the parameters of each model, simulate the latent subtype-level infection dynamics, and perform model validation. Please refer to the [Introductory page](https://github.com/cobeylab/Influenza-immune-dynamics) for a step-by-step overview of the manuscript analysis. Breifly, the worfklow for the model inference is as follows:

1. **Fit the sub-model of the short-term titer dynamics** Navigate to the [Short_term_sub_model](./Short_term_sub_model) directory. Infer the parameters that govern the short-term titer boosting for children and adults with H1N1pdm09 and H3N2. Then, fix these parameters in the full longitudinal models (step 2).

2. **Fit the full single-subtype models of the immune dynamics** Navigate to the [Single_strain_model](./Single_strain_model) directory. First, fit the transmission rate for each subtype to the full data. Then, infer the parameters that govern the dynamics of protection. 

3. **Fit the full multi-subtype model to infer the duration of heterosubtypic protection** Navigate to the [Mulit_strain_model](./Multi_strain_model) directory. Fix the parameters that govern the single-strain dynamics of each subtype based on the output of steps 1 and 2. Then infer the rate of waning of heterosubtypic protection. 

## Simulations
The `Simulations` subdirectories contain the code to simulate the best-fit models using the inferred parameters and to produce the figures that appear in the text from the model simulations.

## References
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.

2. Ionides EL, Breto C, Park J, Smith RA, King AA (2017) Monte Carlo profile confidence
 intervals for dynamic systems. Journal of The Royal Society Interface 14(132).
 


