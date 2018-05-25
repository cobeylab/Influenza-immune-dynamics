
The worfklow of this directory is as follows:

1. Make a list of host objects using the `make_host_data_structure.R` file. This code will preserve the defining features of individuals in the dataset (the initial titer, minimum observed titer, age, birth-date, first visit date, ending visit date, etc.), but it will allow you to simulate many observations for each host. Therefore, the model simulations will sample hosts closely and provide a better approximation of the latent dynamics.

2. Convert the host data list from step 1 to a POMP object using the `data_structure_to_pomp_object_structure.R` file. 

3. Simulate using the `simulate_data.R` file.

4. Make figures using the `make_output_plots.R` file. Figures will appear in the `Output_plots` subdirectory.