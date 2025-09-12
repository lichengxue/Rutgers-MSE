# MSE Notebook

- Notes by and (mostly for the benefit) of RMWJ Bandara

## `bsb_mse_example_01.R`

This is a brief explanation of the example MSE used for black sea bass.
The model sets up a 3 year long MSE with assessments every 3 years. This is 
preceded by 35 years of historical data (1989 to 2023) after which the MSE 
takes over for 10 years. The population is structured as two meta-stocks that
have movement between each other. The following outlines the important
sections of the code.

### Biology of black sea bass



### Time

The historical period of the model is for 35 years (1989 to 2023; `years`).
The MSE is set to run for 10 years (`MSE_years`). Each year is divided to 11 
seasons (`n_seasons`). Then each of these seasons are assigned weighted fractions 
that add up to 1 (`fracyr_seasons`). They spawn in the middle of the year which is 
denoted by `fracyr_spawn`. 

## Movement

This is the most interesting part of how the model is set up. Movement between the 
two subunits is defined by several parameters. `basic_info$NAA_Where` decides which ages 
can be in which subunit on January 1st of each year. Here we have set it up so 
that Age 1 in Stock 1 cannot be in Stock 2 on January 1st and that No individuals 
from Stock 2 move into region 1.

```{r}
basic_info$NAA_where <- array(1, dim = c(n_stocks, n_regions, n_ages))
basic_info$NAA_where[1, 2, 1] <- 0 # Stock 1, age 1 cannot be in region 2 on Jan 1
basic_info$NAA_where[2, 1, ] <- 0  # Stock 2 does not move into region 1
```

Then we define the movement parameters. We say

## `bsb_mse_example_03.R`

