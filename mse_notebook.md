# MSE notebook

*For taking any and all notes regarding MSE for BSB. Posts arranged from latest to oldest*
** Notes by and (mostly for the benefit) of RMWJ Bandara

## `bsb_mse_example_01.R`

This is a brief explanation of the example MSE used for black sea bass.
The model sets up a 3 year long MSE with assessments every 3 years. This is 
preceded by 35 years of historical data (1989 to 2023) after which the MSE 
takes over. The population is structured as two meta-stocks that have movement 
between each other. 

## 08/22/2025 log

Issues related to `bsb_mse_example_01.R` have been fixed. There was an error when 
running the MSE loop (`loop_through_fn`) resulting in non-convergence. This was due to 
not using all the years in the 