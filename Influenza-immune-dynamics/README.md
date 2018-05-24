# Summary
This model tests hypotheses about the dynamics of protective immunity against influenza A viruses by fitting age- and subtype-specific models of influenza immune dynamics to longitudinal serology. 
The model is implemented in R, using an interface with C++ for the dynamic model within the "pomp" [package](http://kingaa.github.io/pomp/install.html)<sup>1</sup>.

# Requirements and Setup 
The inference code was built and run using R version 3.3.2. R can be downloaded [here](https://www.r-project.org).
The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the `installation.R` script. To run this script from the command line, simply navigate to the directory and execute:
```
R CMD BATCH ./installation.R 
```
The packages should be installed in your R library. 

Before running any of the models, users should familiarize themselves with the *pomp* statistical inference software <sup>1</sup>. A helpful introduction can be found [here](https://kingaa.github.io/pomp/vignettes/getting_started.html).

# Raw data
The data that was used in the modeling analysis is available in a Sqlite file in [Data](./Data). 

# Running the Model 
Each model runs separately from a self-contained directory. Navigate to the directory corresponding to the model that you wish to run. The underlying dynamic model and the observation model are specified in `rprocess.R`. The code was written to run each MIF search, or "chain", as a separate process, such that the exploration of the likelihood surface from different starting conditions can be parallelized across computing cores. The `example` subdirectory of each model directory contains an `example_job_submission.sbatch` script to run parallel MIF searches using a high performance computing cluster. The inference proceeds as follows for any candidate model:

2. **Global exploration of the likelihood surface** Maximize the likelihood via a global exploration of the likelihood space from random starting parameter sets.  The`example_global_search_experiment.R` in the `example` folder of each model directory contains example code to set up and run one MIF search from a set of random model parameters. Here you can specify the subtype, host age group, and the parameters of the MIF search. You can run any number of MIF searches in parallel. Each MIF search, or "chain" has an associated chain Id. This process will generate several output files in the `results` folder:

* A `.rda` file that stores the entire MIF chain.
* A `.csv` file containing the output from the search 

After the initial set of MIF iterations for each search, or chain, run the `continue_global_search.R` script to continue chains that have not converged. The parameters of the MIF search, including the random walk of the standard deviation for each parameter, can be updated within this script.

3. **Likelihood profiles** Each inference subdirectory contains the code to construct likelihood profiles to calculate maximum likelihood parameter estiamtes and 95% confidence intervals. To profile over a parameter of interest first update the `example_profile_likelihood.R` script to sweep over the desired range and parameter. Next, construct the profile by generating MIF searches from a series of starting parameter sets that sweep over the desired (fixed) range of the focal parameter. As with the global likelihood search, this script generates one MIF chain for one profile point, and multiple profile points can be run in parallel by specifying a series of "chainId" variables. The output, generated in the `results` folder, consists of:
an associated chain Id. This process will generate several output files in the `results` folder:
* A `.rda` file that stores the entire MIF chain for the profile point.
* A `.csv` file containing the output from the search (just in case there are overwrite issues with the SQLITE database due to many chains being run in parallel.)

## Calculating Confidence Intervals from Likelihood Profiles
Once the profile likelihood search has been completed, select the point of maximum likelihood for each value of the profile parameter to represent the inferred parameter. Then, use the Monte Carlo Adjusted Profile (MCAP) method<sup>2</sup> to calculate a smoothed estimate of the profile and the corresponding 95% confidence interval. A function containing the MCAP algorithm is given in the `model_functions.R` script within the `Utility_scripts` folder of each model subdirectory. The `calculate_MLE_and_CI_for_param.R` script in the `analysis` folder calculates the smoothed profile, MLE, and 95% CI for each parameter inferred.

## Figures from model results 

## Simulations 

## References
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.

2. Ionides EL, Breto C, Park J, Smith RA, King AA (2017) Monte Carlo profile confidence
 intervals for dynamic systems. Journal of The Royal Society Interface 14(132).
 


