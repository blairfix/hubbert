#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat prod_sim (
		     arma::vec & well_id_vec,
		     arma::vec & date_start_vec,
		     double year_start,
		     double reserve_rate,
		     double reserve_intercept,
		     double half_rate,
		     double half_intercept,
		     double n_year_extrap,
		     double stop_threshold,
		     double reserve_sd_log
		    )

{

    // make time vector 
    double t_start = 0;
    double n_month = n_year_extrap * 12 + 1;

    arma::vec t_vec = arma::linspace( 
				      t_start,
				      t_start + n_year_extrap ,
				      n_month 
				     );

    // estimate size of output container
    double n_wells = well_id_vec.size();
    int nrows = n_month * n_wells;

    // output matrix (empty vals set to -999)
    arma::mat output(nrows, 4, arma::fill::value(-999) );
    int id_out = 0; // output index counter


    // -----------------------------------------------------
    // loop over wells

    for( int i = 0; i < n_wells; i++ ){
	
	// well data
	double well_id = well_id_vec[i];
	double date_start = date_start_vec[i];

	// get half life from well start date
	double half_life = std::exp( half_rate * date_start + half_intercept );

	// convert half life to exponential decline rate
	double decline_rate = -log(2.0) / half_life; 
	
	// get reserve size from start date
	double reserve = std::exp ( reserve_rate * date_start + reserve_intercept );

	// add lognormal noise to reserve size
	// lognormal distribution with mean = 0
	double sigma = reserve_sd_log;
	double mu = log(1) - 0.5 * pow( sigma, 2);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::lognormal_distribution<> d(mu, sigma);

	reserve = reserve *  d(gen);


	// loop over t_vec
	arma::vec reserve_extrap( t_vec.size() );

	for( int t = 0; t < t_vec.size(); t++ ){
	    
	    // extrapolate future remaining reserve using half life function
	    reserve_extrap[t] =  reserve * std::exp( decline_rate * t_vec[t] );

	    // if t > 0, get monthly production data 
	    if( t > 0 ){

		// monthly production = difference between consecutive reserve values
		double month_prod =  reserve_extrap[t-1] - reserve_extrap[t];
		
		// test if monthly production greater than stop threshold
		if( month_prod > stop_threshold){

		    // make current year 
		    double year = t_vec[t] + year_start;
		    
		    // write to output containers
		    output(id_out, 0) = well_id;
		    output(id_out, 1) = year;
		    output(id_out, 2) = month_prod; 
		    output(id_out, 3) = date_start;
		    
		    // update output iterator
		    id_out++;	

		}
	    }
	}
    }


    // remove empty values from output matrix
    arma::uvec id_keep = find(output.col(0) != -999 );
    output = output.rows( id_keep ) ;

    // output matrix contains following columns:
    // 1. well id
    // 2. year
    // 3. monthly production
    // 4. well start date

    return output;
}



