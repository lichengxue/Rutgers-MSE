# reusable_functions.R #
# Set of reusable functions 


# Reconcile index timing
# This function is used to convert when 
# You don't need to understand this function, it's used to back calculate fracyr_indices...
reconcile_index_timing <- function(fracyr_seasons, fracyr_indices, index_seasons,
                                   boundary_rule = c("next_season", "this_season")) {
  boundary_rule <- match.arg(boundary_rule)
  n_seasons <- length(fracyr_seasons)
  stopifnot(length(index_seasons) == length(fracyr_indices))
  if (abs(sum(fracyr_seasons) - 1) > 1e-8) stop("fracyr_seasons must sum to 1.")
  if (any(index_seasons < 1 | index_seasons > n_seasons)) stop("index_seasons out of range.")
  
  starts <- c(0, cumsum(fracyr_seasons))     # season starts, length n_seasons+1
  eps <- sqrt(.Machine$double.eps)
  
  # Clamp absolute times into [0, 1)
  x <- pmin(pmax(fracyr_indices, 0), 1 - eps)
  
  # Implied seasons from the absolute times (for diagnostics)
  x_nudge <- if (boundary_rule == "next_season") {
    y <- x + eps
    y[x == 0] <- eps
    y
  } else {
    pmax(x - eps, 0)
  }
  implied_season <- findInterval(x_nudge, starts, rightmost.closed = FALSE)
  implied_season[implied_season < 1] <- 1
  implied_season[implied_season > n_seasons] <- n_seasons
  
  # Use the provided index_seasons but report mismatches
  mismatch <- which(implied_season != index_seasons)
  if (length(mismatch)) {
    msg <- paste0(
      "Note: ", length(mismatch), " index/indices have implied_season != index_seasons: [",
      paste(mismatch, collapse = ", "), "]. Using index_seasons as authoritative."
    )
    message(msg)
  }
  
  # Offsets within the chosen seasons
  season_len <- fracyr_seasons[index_seasons]
  offset <- x - starts[index_seasons]   # may be ~0 at boundaries
  # Clean tiny negatives/positives; clamp to [0, season_len)
  offset[abs(offset) < eps] <- 0
  offset <- pmin(pmax(offset, 0), season_len - eps)
  
  # Proportion within season (0..1)
  prop_within <- ifelse(season_len > 0, offset / season_len, 0)
  
  # Rebuild absolute times from the chosen decomposition
  fracyr_rebuilt <- starts[index_seasons] + offset
  fracyr_rebuilt <- pmin(pmax(fracyr_rebuilt, 0), 1 - eps)
  
  list(
    starts = starts,
    implied_season = implied_season,
    used_season = index_seasons,
    within_season_offset = offset,        # absolute fraction of year inside season
    within_season_prop = prop_within,     # 0..1 within-season proportion
    fracyr_rebuilt = fracyr_rebuilt       # should match input up to ~1e-8
  )
}

#' @title 
#' Generate inverse logit
#' 
#' @description
#' `sum` returns the sum of all the values present in its arguments.
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
# inverse logit
gen.invlogit <- function(eta, low, upp, s = 1) {
  low + (upp - low) * plogis(s * eta)
}


#' @title 
#' Gaussian link
#' @description
#' Return the adjustment factor for recruitment based on a gaussian link
#' between recruitment and temperature
#' @details
#' This function should be supplied both a temperature, an optimal temperature
#' and the width of the relationship
gauss_rec <- function(tx, topt, twidth) {
  exp(-0.5*((tx-topt)/twidth)^2)
}

