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

### Spatial considerations

It's important to remember that *stock* and *region* have different meanings in 

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

```r
basic_info$NAA_where <- array(1, dim = c(n_stocks, n_regions, n_ages))
basic_info$NAA_where[1, 2, 1] <- 0 # Stock 1, age 1 cannot be in region 2 on Jan 1
basic_info$NAA_where[2, 1, ] <- 0  # Stock 2 does not move into region 1
```

Movement parameters are defined using the `move` list object. It has the following 
items.
- `move$stock_move` - Says whether the individuals can move to and from that subunit
- `move$separable` - Says whether the subunits are separate from one another
- `move$must_move` - An array of dimensions `n_stocks` x `n_seasons` x `n_regions`. 
This defines whether there is mandatory movement of individuals belonging to a 
certain stock from one region to the other. 
- `move$can_move` - An array of dimensions  `n_stocks` x `n_seasons` x `n_regions` x `n_regions`. 
This defines whether there can be movement of individuals belong to a stock from 
one region to another and vice versa. 

> All individuals belonging to stock 1 who are in region 2, should move
back to region 1 in season 5 (May).
> We also allow for movement 

## `bsb_mse_example_03.R`

