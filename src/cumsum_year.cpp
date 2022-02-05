#include <RcppArmadillo.h>
#include <list>
#include <math.h>
#include <vector>

#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::vec cumsum_by_date (
    const arma::vec &x,
    const arma::vec &year,
    double steps_per_year
    )

{


    // year min
    double y_min = arma::min(year);
    y_min = floor(y_min * steps_per_year ) /  steps_per_year;

    // year max
    double y_max = arma::max(year);
    y_max = ceil(y_max * steps_per_year ) /  steps_per_year;


    // make bin vector
    double year_span = y_max - y_min;
    int step = steps_per_year + 1;
    int n_bins = (y_max - y_min) * step + 1;


    arma::vec bin_vector(n_bins);
    bin_vector.zeros();


    // get cumulative sum by year
    arma::vec output( x.size() );

    for(int i = 0; i < x.size(); ++i){

        // bin id
        int id = ( year[i] - y_min ) * step ;

//        if(id >= bin_vector.size() | id < 0 ){
//            std::cout << id << std::endl;
//            std::cout << bin_vector.size() << std::endl;
//            throw std::invalid_argument( "bin overeach" );
//        }

        // sum production
        bin_vector[id] = bin_vector[id] + x[i];

        // output
        output[i] = bin_vector[id];

    }

    return output;

}


