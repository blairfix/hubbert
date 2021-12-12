#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat depletion_stats (
    const arma::vec &depletion,
    const arma::vec &well_id,
    const arma::vec &month
    )

{

    // output matrix
    arma::vec wells_unique = arma::unique(well_id);
    int n_wells = wells_unique.size();
    int max_t = arma::max(month) + 2;

    arma::mat output(n_wells, max_t);
    output.zeros();

    // initial value
    output(0, 0) = depletion[0];
    int t_counter = 1;
    int well_counter = 0;

    // loop over depletion data
    for(int i = 1; i < depletion.size(); i++){

        // if more well data, iterate row number
        // else iterate column number, set column to 0
        if( well_id[i] == well_id[i-1] ){

            output(well_counter, t_counter) = depletion[i];
            t_counter++;

        } else{

            well_counter++;
            t_counter = 0;

            output(well_counter, t_counter) = depletion[i];

            t_counter++;

        }

    }

    // get stats by column
    arma::mat stats(2, max_t);

    stats.row(0) = arma::mean( output, 0 );
    stats.row(1) = arma::stddev( output, 0 );


    // transpose to long matrix
    arma::mat stats_long = stats.t();

    return stats_long;

}


