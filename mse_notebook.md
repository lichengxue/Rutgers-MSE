# MSE Notebook

- Notes by and (mostly for the benefit) of RMWJ Bandara

Background readings to fully understand WHAM (Woods Hall Assessment Model) and 
it's derivatives.


## `bsb_mse_example_01.R`

This is a brief explanation of the example MSE used for black sea bass.
The model sets up a 3 year long MSE with assessments every 3 years. This is 
preceded by 35 years of historical data (1989 to 2023) after which the MSE 
takes over for 10 years. The population is structured as two meta-stocks that
have movement between each other. The following outlines the important
sections of the code.


### Biology of black sea bass

### Spatial considerations

It's important to remember that *stock* and *region* have different meanings in the model.
The `stock` refers to the biological unit and it's individuals. 

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
- `move$separable` - If `separable` is equal to TRUE, fish will only move at the 
end of the time box.
- `move$must_move` - An array of dimensions `n_stocks` x `n_seasons` x `n_regions`. 
This defines whether there is mandatory movement of individuals belonging to a 
certain stock from one region to the other. 
- `move$can_move` - An array of dimensions  `n_stocks` x `n_seasons` x `n_regions` x `n_regions`. 
This defines whether there can be movement of individuals belong to a stock from 
one region to another. 

Only individuals belonging to stock 1 can move between regions.
All individuals belonging to stock 1 who are in region 2, should move
back to region 1 in season 5 (May).
We also allow for movement of stock 1 individuals from region 2 to region 1 
from January to April. We allow for movement of stock 1 individuals from 
region 1 to region 2 from August to December. We allow the movement of all 
individuals belonging to stock 1 in region 2 to move to region 1.
Then we set  movement model for each stock (`mean_vals`). Then we set the 
movement model (*stock_constant*). Reference for the movement models can be found 
[here](https://timjmiller.github.io/wham/reference/set_move.html).

> "stock_constant" - estimate a movement rate for each stock to each region shared across all 
seasons, ages, years

```r
move <- list(stock_move = c(TRUE, FALSE), separable = TRUE)
move$must_move <- array(0, dim = c(2, 11, 2))
move$must_move[1, 5, 2] <- 1
move$can_move <- array(0, dim = c(2, 11, 2, 2))
move$can_move[1, 1:4, 2, 1] <- 1
move$can_move[1, 7:11, 1, 2] <- 1
move$can_move[1, 5, 2, ] <- 1
mean_vals <- array(0, dim = c(2, length(fracyr_seasons), 2, 1))
mean_vals[1, 1:11, 1, 1] <- 0.02214863
mean_vals[1, 1:11, 2, 1] <- 0.3130358
move$mean_vals <- mean_vals
move$mean_model <- matrix("stock_constant", 2, 1)
```

## Numbers at age dynamics

Recruitment is random around a mean (`recruit_model = 2`). `sigma`($\sigma$) determines 
the type of error associated with recruitment and how numbers at age relate to
the previous age. The options are either "rec" or "rec+1" where "rec+1" estimates 
two sigmas ($\sigma_a, \sigma$). One for recruitment and one shared among other 
ages. <br>
Correlation between age classes is determined by `re_cor`. We have set it to "iid", 
where numbers at age are uncorrelated with each other (The errors are drawn from 
Independent and Identical distribution). <br>
Initial numbers at age ($NAA_1$) is determined by `ini.opt`. We have set it to 
"age_specific_fe" which is age-and-region specific fixed effects parameters.
The initial standard deviation values for NAA is defined by `sigma_vals`. 
Mean recruitment ($\mu$) for each stock is determined by `recruit_pars`.

```r
sigma <- "rec+1" # Full state-space model where all numbers at age are random effects
re_cor <- "iid" # iid - Independent and identically distributed covariate
ini.opt <- "age-specific-fe" # Need to figure out what this means....
sigma_vals <- array(0.2, dim = c(n_stocks, n_regions, n_ages)) # NAA survival sigma - Another potential fix is increasing this value
sigma_vals[, , 1] <- 0.75 # Recruitment sigma/error

NAA_re <- list(
  N1_model = rep(ini.opt, n_stocks),
  sigma = rep(sigma, n_stocks),
  cor = rep(re_cor, n_stocks),
  recruit_model = 2, # 2 = density-independent. Random about mean
  recruit_pars = list(16000, 21000), # Mean recruitment for each stock
  sigma_vals = sigma_vals
)
```

### Natural mortality

Natural mortality ($M$) for black sea bass is set as a time-invariant, 
age constant value (0.4).

### Gear selectivity

Gear selectivity is defined for all four different fleets (commercial and recreational) 
and for the two trawl surveys.

### Building the Operating Model (OM) object

Then all the above information is gathered into a single list for wham_input (function is `prepare_wham_input`). 
Then this input list is edited because it has some nice defaults which makes it 
possible to make minor changes to get the model that you want.  
Weights at age for this input list is updated using `waa_info`. 
Then initial numbers at age ($NAA_{a,s}$) are set using numbers from the stock assessment.

Then the operating model (OM) is generated using function 
[`fit_wham`](https://timjmiller.github.io/wham/reference/fit_wham.html) without 
actually fitting it to data using the following options.

- `do.fit` - Fit the model using TMB
- `do.brps` - Calculate biological reference points
- `MakeADFun.silent` - If to fit log likelihoods in TMB to the fitted model.

We are just generating an operating model with the data and doing no fitting.

```r
om <- fit_wham(input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
```

### Setting up the MSE and the estimation model

We set up the MSE to run every 3 years starting from the year where the historical 
data ends. 

Reference for specifying numbers at age in wham is 
[here](https://timjmiller.github.io/wham/reference/set_NAA.html)

### Specify MSE timeline and HCR

A single realization of the OM is generated with `update_om_fn()`. 

#### Specify Harvest Control Rule

HCRs are specified as a list. WhamMSE supports three types of HCRs. These include 

1. Fixed Percentage of SPR - `hcr.type <- 1`
2. Constant catch - `hcr.type <- 2`
3. Hockey stick / Sliding - `hcr.type <- 3`

Each of these must be accompanied by `hcr.opts`.

```r
  hcr <- list()
  hcr$hcr.type <- 1 # FXSPR - Fishing pressure to keep the SPR at a certain percentage
  hcr$hcr.opts <- list(use_FXSPR = TRUE, percentFXSPR = 75) # Apply F at 75% unfished SPR
```



### Run the MSE loop

MSE loop is run through `loop_through_fn`. 

```r
  mod <- loop_through_fn(
    om = om_with_data,
    em_info = info,
    random = random,
    sel_em = sel,
    M_em = M,
    NAA_re_em = NAA_re,
    move_em = move_em,
    em.opt = list(
      separate.em = FALSE,
      separate.em.type = 3,
      do.move = TRUE,
      est.move = TRUE
    ),
    assess_years = assess.years,
    assess_interval = assess.interval,
    base_years = base.years,
    year.use = 35,
    add.years = TRUE,
    seed = 123,
    hcr = hcr,
    save.last.em = TRUE
  )
```

Here both the om and em information must be provided along with information on which 
years the assessment happens. 

### Known issues

The model can run into convergence issues. Potential fixes are shown below.

- Changing the years-at-age for Year 1 in the estimation model to _equilibrium_
  - ```r
      NAA_re$N1_model[] = "equilibrium"
    ```
- Increasing `sigma_vals` for the operating model to 0.75 to reflect the real stock assessment
- Increase `prior_sigma` of the estimation model's movement 
- Change whether we want to estimate movement with the MSE function (`loop_through_fn`)
- Change the operating model's population dynamics model from a state-space model with random effects to 
just a state-space model. This needs to be reflected in the estimation model as well. 
    - ```r
      sigma <- "rec"
      re_cor <- "iid"
      ini.opt <- "age-specific-fe"
      sigma_vals <- array(0.2, dim = c(n_stocks, n_regions, n_ages)) # NAA survival sigma
      sigma_vals[, , 1] <- 0.75 # Recruitment sigma
      # For the estimation model
      NAA_re$sigma = "rec"
      ```

