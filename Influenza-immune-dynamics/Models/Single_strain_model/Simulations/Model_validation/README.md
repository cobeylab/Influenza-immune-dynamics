
The worfklow of this directory is as follows:

1. Simulate using the `simulate_data.R` file. The POMP object file for any simulation should be the corresponding POMP object that you created in the `Inference` subdirectory. The goal here is to simulate titers for each individual and compare them to the observed titers at the same visit dates.

2. Make figures using the `make_output_plots.R` file. Figures will appear in the `Output_plots` subdirectory.