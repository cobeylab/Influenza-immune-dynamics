# Summary
This model tests hypotheses about the dynamics of protective immunity against influenza A viruses by fitting age- and subtype-specific models of influenza immune dynamics to longitudinal serology. 
The model is implemented in R, using an interface with C++ for the dynamic model within the "pomp" [package](http://kingaa.github.io/pomp/install.html)<sup>1</sup>.

# Requirements and Setup 
The inference code was built and run using R version 3.3.2. R can be downloaded [here](https://www.r-project.org).
The approximate installation time is one hour.
The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the `installation.R` script. To run this script from the command line, simply navigate to the directory and execute:
```
R CMD BATCH ./installation.R 
```
The packages should be installed in your R library. 

Before running any of the models, users should familiarize themselves with the *pomp* statistical inference software <sup>1</sup>. A helpful introduction can be found [here](https://kingaa.github.io/pomp/vignettes/getting_started.html).

# Raw data
The data that was used in the modeling analysis is available in a Sqlite file in [Data](./Data). 

# Workflow for Subtype-level Analysis 
This section outlines the workflow of the analysis used in the manuscript to estimate the model parameters and generate the results and figures. Follow these steps separately for each subtype (H3N2 and H1N1pdm09).  

1. **Fit the sub-model of the short-term post-infection titer dynamics**  The first step of the analysis is to fit the sub-model of the short-term post-infection titer dynamics to the subset of PCR-confirmed infections. The data include the time of PCR-confirmed infection and the immediate pre- and post-infection titers for a subset of individuals. The magnitude and standard deviation of the short-term titer boost and the dependence of the titer boost on the pre-infection titer are estimated at this stage for children and adults. Navigate to the [Short_term_sub_model](./Models/Short_term_sub_model) directory and follow the instructions. Parameters are estimated using likelihood profiles to estimate the MLE and the 95% confidence interval. Details of the parameter estimation are provided in the next section ("Running the Model to Estimate Parameters"). 

2. **Fit the parameters of the full model to titer data from the entire cohort** The next step is to fit the parameters of the full model to the longitudinal titers from the entire cohort, fixing the parameters associated with short-term titer dynamics from step 1. The contribution of HI-correlated and non-HI-correlated protection, the titer waning rate, the 50% protective titer, and the long-term boost after infection are estimated at this stage. Navigate to the [Single_strain_model](./Models/Single_strain_model) directory and follow the instructions. A detailed overview of the parameter estimation is provided in the next section ("Running the Model to Estimate Parameters").

3. **Simulate the latent infection dynamics of infection** After all model parameters have been estimated for a given subtype (steps 1 and 2), the next step is to simulate synthetic data using the best-fit parameters. Navigate to the [Simulations](./Models/Single_strain_model/Simulations/Latent_states) directory and follow the directions to simulate the latent dynamics.  

4. **Generate figures from the synthetic data** Using the synthetic data from step 3, the next step is to plot the individual-level susceptibility and the epidemic behavior. Navigate to the [Simulations](./Models/Single_strain_model/Simulations/Latent_states) directory and execute the `make_output_plots.R` script. Output plots will appear as `.pdf` files in the `Output_plots` subdirectory.

5. **Perform model validation** The next step is to generate plots to compare the simulated data for each individual to the corresponding observed data. Navigate to the  [Single_strain_model](./Models/Single_strain_model/Simulations/Model_Validation and follow the instructions. Output plots will appear as `.pdf` files in the `Output_plots` subdirectory.

6. **Estimate the rate of heterosubtypic protection** The final step is to estimate the rate of waning of heterosubtypic protection using a multi-strain model, fixing the subtype-specific parameters. Therefore, this step must be performed after steps 1-2 have been completed for each subtype. Navigate to the [Multi_strain_model](./Models/Multi_strain_model) directory and follow the instructions to estimate the rate of heterosubtypic immune waning. A detailed overview of the parameter estimation is provided in the next section ("Running the Model to Estimate Parameters").

# Running the Model to Estimate Parameters 
The [Models](./Models) directory contains the code to run the inference for the sub-model of the short-term titer dynamics, the full model of the longitudinal titer dynamics for either H1N1pdm09 or H3N2, and the two-strain model of heterosubytpic immunity. Each model runs separately from a self-contained sub-directory. Please consult the `README` files in the `Models` directory to guide the workflow for that particular inference. Navigate to the directory corresponding to the model that you wish to run. The underlying dynamic model and the observation model are specified in the `rprocess` R files. The code was written to run each MIF search, or "chain", as a separate process, such that the exploration of the likelihood surface from different starting conditions can be parallelized across computing cores.  Each model directory contains an `example_job_submission.sbatch` script to run parallel MIF searches using a high performance computing cluster. The components of the inference for any particular model are as follows:

1. **Global exploration of the likelihood surface** Maximize the likelihood via a global exploration of the likelihood space from random starting parameter sets.  The`example_global_likelihood.R` script in the `Inference` folder of the short-term and single-strain model directoires contains example code to set up and run one MIF search from a set of random model parameters. Here you can specify the subtype, host age group, and the parameters of the MIF search. You can run any number of MIF searches in parallel. Each MIF search, or "chain" has an associated chain Id. This process will generate several output files:

* A `.rda` file that stores the entire MIF chain.
* A `.csv` file containing the output from the search.


2. **Likelihood profiles** Each model subdirectory contains the code to construct likelihood profiles to calculate maximum likelihood parameter estiamtes and 95% confidence intervals. To profile over a parameter of interest first update the `example_profile.R` script to sweep over the desired range and parameter. Next, construct the profile by generating MIF searches from a series of starting parameter sets that sweep over the desired (fixed) range of the focal parameter. As with the global likelihood search, this script generates one MIF chain for one profile point, and multiple profile points can be run in parallel by specifying a series of "chainId" variables. The output consists of:
* A `.rda` file that stores the entire MIF chain for the profile point.
* A `.csv` file containing the output from the search.

3. **Calculating Confidence Intervals from Likelihood Profiles**
Once the profile likelihood search has been completed, select the point of maximum likelihood for each value of the profile parameter to represent the inferred parameter. Then, use the Monte Carlo Adjusted Profile (MCAP) method<sup>2</sup> to calculate a smoothed estimate of the profile and the corresponding 95% confidence interval. A function containing the MCAP algorithm is given in the `model_functions.R` script within the `Utility_scripts` folder of each model subdirectory. 


# References
1. King AA, Nguyen D and Ionides EL (2015) Statistical inference for partially observed Markov processes via the R package pomp. arXiv preprint arXiv:1509.00503.

2. Ionides EL, Breto C, Park J, Smith RA, King AA (2017) Monte Carlo profile confidence
 intervals for dynamic systems. Journal of The Royal Society Interface 14(132).
 


