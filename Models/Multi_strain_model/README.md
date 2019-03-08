
Use this directory to erform the model inference (using the code in the `Inference` subdirectory). Specifically, estimate the rate of waning of heterosubtypic protection by constructing a likelihood profile. Please refer to the [Introductory page](../) for a step-by-step overview of the manuscript analysis. Briefly, the overview of the inference is as follows: 

1. Execute the `make_host_data_list.R` script to generate a list of host data objects

2. Convert the host data list from step 1 to a POMP object using the `data_to_pomp_object_structure.R` file.

3. Fixing the subtype-specific parameters from the single-strain models for H1N1pdm09 and H3N2, execute the `experiment_profile_w_hetero.R` to generate the likelihood profile for the rate of waning of heterosubtypic immunity. 
