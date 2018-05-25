
The worfklow of this directory is as follows:

1. Perform the model inference (using the code in the `Inference` subdirectory). 

2. Analyze the inference results and calculate the parameter MLEs/95% CIs (using the code in the `Results_analysis` subdirectory). After you have estimated all of the fitted parameters, create a `.rda` file containing a named
vector of all parameters for the model simulations (step 3).

3. Simulate the dynamics from the best-fit model parameters and generate the output plots (using code in the `Simulations` subdirectory).