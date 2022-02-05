# hubbert

This is an R package that contains tools I've developed for analyzing oil-and-gas production and depletion.

### Functions

* `cumsum_year`: a function for stacking the production curves of individual wells.

* `depletion_stats`: a function for calculating the average depletion curve across all wells.

* `lotka_volterra`: a fast implementation of the oil-and-gas version of the Lotka-Volterra model proposed by [Ugo Bardi and Alessandro Lavacchi](https://www.mdpi.com/1996-1073/2/3/646). I use it fit model parameters to empirical data.

* `prod_curve`: a function for estimating production curves for individual wells. You can read more details [here](https://github.com/blairfix/prod_curve).

* `prod_extrap`: a function for extrapolating the production curves of active wells into the future.

* `prod_sim`: a function for simulating the production curves for future wells.
