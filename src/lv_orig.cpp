#include <Rcpp.h>
#include <iostream>
#include <math.h>

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]

Rcpp::DataFrame  lv_orig (
	double prey_initial,
	double predator_initial,
	double prey_growth_rate,
	double predator_death_rate,
	double prey_effect_on_predator,
	double predator_effect_on_prey,
	double time,
	double step_size
	)

{

    int n_steps = std::round(time / step_size);

    // output vectors
    Rcpp::NumericVector time_out(n_steps);
    Rcpp::NumericVector prey_pop_out(n_steps);
    Rcpp::NumericVector predator_pop_out(n_steps);
    Rcpp::NumericVector delta_prey_out(n_steps);
    Rcpp::NumericVector delta_predator_out(n_steps);


    // initial conditions
    double predator_pop = predator_initial;
    double prey_pop = prey_initial;

    // rescale parameters for step size
    prey_growth_rate = prey_growth_rate * step_size;
    predator_death_rate = predator_death_rate * step_size;
    predator_effect_on_prey = predator_effect_on_prey * step_size;
    prey_effect_on_predator = prey_effect_on_predator * step_size;

    // run model
    for( int i = 0; i < n_steps; i++ ){

        //  prey population change
        double delta_prey = ( prey_growth_rate * prey_pop  ) - ( predator_effect_on_prey * prey_pop * predator_pop );

        double delta_predator =  (- predator_death_rate * predator_pop ) + ( prey_effect_on_predator * prey_pop * predator_pop );

        // new populations
        prey_pop  = prey_pop +  delta_prey;
        predator_pop  = predator_pop + delta_predator;

	// check for negatives
        if(prey_pop < 0){ prey_pop = 0; }
        if(predator_pop < 0){ predator_pop = 0; }

        // export
        time_out[i] = i*step_size;
        prey_pop_out[i] = prey_pop;
	predator_pop_out[i] = predator_pop;
        delta_prey_out[i] = delta_prey;
        delta_predator_out[i] = delta_predator;

    }

    Rcpp::DataFrame output = Rcpp::DataFrame::create(
	    Rcpp::Named("time") = time_out,
	    Rcpp::Named("prey_population") = prey_pop_out,
	    Rcpp::Named("predator_population") = predator_pop_out,
	    Rcpp::Named("delta_prey") = delta_prey_out,
	    Rcpp::Named("delta_predator") = delta_predator_out
	    );

    return output;

}
