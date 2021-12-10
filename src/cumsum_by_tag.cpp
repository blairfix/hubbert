#include <RcppArmadillo.h>
#include <list>
#include <map>
#include <math.h>
#include <string>
#include <vector>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


std::list<double> cumsum_by_tag (
    const arma::vec &x,
    const arma::vec &tag
    )

{


    std::list<double> output;
    std::map< std::string, double> counter;

    // initialize sums to 0
    arma::vec tag_unique = arma::unique(tag);

    for(int i = 0; i < tag_unique.size(); ++i){

        // tag to string
        std::string tag_string = std::to_string( tag_unique[i] );

        // set init to 0
        counter[tag_string] = 0;


    }



    // loop over data
    // get cumulative sum by tag
    for(int i = 0; i < x.size(); ++i){

        // tag to string
        std::string tag_string = std::to_string( tag[i] );

        // cumulative sum
        counter[tag_string] = counter[tag_string] + x[i];

        // append to output
        output.push_back( counter[tag_string] ) ;

    }


    return output;

}


