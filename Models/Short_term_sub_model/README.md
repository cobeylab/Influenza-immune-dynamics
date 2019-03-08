This directory contains all of the code necessary to estimate the parameters of the sub-model of the short-term titer dynamics, to simulate the dynamics from the best-fit sub-model parameters, and to perform model validation. Please refer to the [Introductory page](../) for a step-by-step overview of the manuscript analysis. Briefly, the worfklow for this directory is as follows:

1. Estimate the parameters that guide the short-term post-infection titer dynamics (using the code in the `Inference` subdirectory). 

2. Analyze the inference results and calculate the parameter MLEs/95% CIs (using the code in the `Results_analysis` subdirectory). After you have estimated all of the fitted parameters, create a `.rda` file containing a named vector of all parameters for the model simulations (step 3).

3. Simulate the dynamics from the best-fit model parameters and generate the output plots (using code in the `Simulations` subdirectory). The output plots will save as `.pdf` files in the `Output_plots` subdirectory. 