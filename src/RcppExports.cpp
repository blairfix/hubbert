// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cumsum_by_tag
std::list<double> cumsum_by_tag(const arma::vec& x, const arma::vec& tag);
RcppExport SEXP _hubbert_cumsum_by_tag(SEXP xSEXP, SEXP tagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tag(tagSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum_by_tag(x, tag));
    return rcpp_result_gen;
END_RCPP
}
// cumsum_by_date
arma::vec cumsum_by_date(const arma::vec& x, const arma::vec& year, double steps_per_year);
RcppExport SEXP _hubbert_cumsum_by_date(SEXP xSEXP, SEXP yearSEXP, SEXP steps_per_yearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type year(yearSEXP);
    Rcpp::traits::input_parameter< double >::type steps_per_year(steps_per_yearSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum_by_date(x, year, steps_per_year));
    return rcpp_result_gen;
END_RCPP
}
// depletion_stats
arma::mat depletion_stats(const arma::vec& depletion, const arma::vec& well_id, const arma::vec& month);
RcppExport SEXP _hubbert_depletion_stats(SEXP depletionSEXP, SEXP well_idSEXP, SEXP monthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type depletion(depletionSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type well_id(well_idSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type month(monthSEXP);
    rcpp_result_gen = Rcpp::wrap(depletion_stats(depletion, well_id, month));
    return rcpp_result_gen;
END_RCPP
}
// lotka_volt
arma::mat lotka_volt(double r_0, double c_0, double k_1, double k_2, double k_3, double time, double step_size);
RcppExport SEXP _hubbert_lotka_volt(SEXP r_0SEXP, SEXP c_0SEXP, SEXP k_1SEXP, SEXP k_2SEXP, SEXP k_3SEXP, SEXP timeSEXP, SEXP step_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r_0(r_0SEXP);
    Rcpp::traits::input_parameter< double >::type c_0(c_0SEXP);
    Rcpp::traits::input_parameter< double >::type k_1(k_1SEXP);
    Rcpp::traits::input_parameter< double >::type k_2(k_2SEXP);
    Rcpp::traits::input_parameter< double >::type k_3(k_3SEXP);
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(lotka_volt(r_0, c_0, k_1, k_2, k_3, time, step_size));
    return rcpp_result_gen;
END_RCPP
}
// lotka_volt_cumulative
arma::mat lotka_volt_cumulative(double r_0, double c_0, double k_1, double k_2, double k_3_init, double k_3_slope, double time, double step_size);
RcppExport SEXP _hubbert_lotka_volt_cumulative(SEXP r_0SEXP, SEXP c_0SEXP, SEXP k_1SEXP, SEXP k_2SEXP, SEXP k_3_initSEXP, SEXP k_3_slopeSEXP, SEXP timeSEXP, SEXP step_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r_0(r_0SEXP);
    Rcpp::traits::input_parameter< double >::type c_0(c_0SEXP);
    Rcpp::traits::input_parameter< double >::type k_1(k_1SEXP);
    Rcpp::traits::input_parameter< double >::type k_2(k_2SEXP);
    Rcpp::traits::input_parameter< double >::type k_3_init(k_3_initSEXP);
    Rcpp::traits::input_parameter< double >::type k_3_slope(k_3_slopeSEXP);
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(lotka_volt_cumulative(r_0, c_0, k_1, k_2, k_3_init, k_3_slope, time, step_size));
    return rcpp_result_gen;
END_RCPP
}
// prod_curve
arma::mat prod_curve(const arma::uvec& well_id, const arma::vec& prod_first_6, const arma::vec& prod_first_12, const arma::vec& prod_first_24, const arma::vec& prod_first_60, const arma::vec& prod_2nd_last_year, const arma::vec& prod_last_year, const arma::vec& prod_daily_last_year, const arma::vec& prod_cumulative, const arma::vec& prod_peak_daily, const arma::vec& date_start, const arma::vec& date_end, const arma::vec& peak_month, const arma::vec& decline_3, const arma::vec& decline_12, const arma::vec& decline_24, const arma::vec& decline_60, double time_step, double step_threshold);
RcppExport SEXP _hubbert_prod_curve(SEXP well_idSEXP, SEXP prod_first_6SEXP, SEXP prod_first_12SEXP, SEXP prod_first_24SEXP, SEXP prod_first_60SEXP, SEXP prod_2nd_last_yearSEXP, SEXP prod_last_yearSEXP, SEXP prod_daily_last_yearSEXP, SEXP prod_cumulativeSEXP, SEXP prod_peak_dailySEXP, SEXP date_startSEXP, SEXP date_endSEXP, SEXP peak_monthSEXP, SEXP decline_3SEXP, SEXP decline_12SEXP, SEXP decline_24SEXP, SEXP decline_60SEXP, SEXP time_stepSEXP, SEXP step_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type well_id(well_idSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_first_6(prod_first_6SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_first_12(prod_first_12SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_first_24(prod_first_24SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_first_60(prod_first_60SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_2nd_last_year(prod_2nd_last_yearSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_last_year(prod_last_yearSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_daily_last_year(prod_daily_last_yearSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_cumulative(prod_cumulativeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prod_peak_daily(prod_peak_dailySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type date_start(date_startSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type date_end(date_endSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type peak_month(peak_monthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type decline_3(decline_3SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type decline_12(decline_12SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type decline_24(decline_24SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type decline_60(decline_60SEXP);
    Rcpp::traits::input_parameter< double >::type time_step(time_stepSEXP);
    Rcpp::traits::input_parameter< double >::type step_threshold(step_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(prod_curve(well_id, prod_first_6, prod_first_12, prod_first_24, prod_first_60, prod_2nd_last_year, prod_last_year, prod_daily_last_year, prod_cumulative, prod_peak_daily, date_start, date_end, peak_month, decline_3, decline_12, decline_24, decline_60, time_step, step_threshold));
    return rcpp_result_gen;
END_RCPP
}
// prod_extrap
arma::mat prod_extrap(arma::vec& well_id_vec, arma::vec& cum_prod_vec, arma::vec& date_start_vec, arma::vec& month_prod_vec, double year_last, double reserve_rate, double reserve_intercept, double half_rate, double half_intercept, double n_year_extrap, double stop_threshold);
RcppExport SEXP _hubbert_prod_extrap(SEXP well_id_vecSEXP, SEXP cum_prod_vecSEXP, SEXP date_start_vecSEXP, SEXP month_prod_vecSEXP, SEXP year_lastSEXP, SEXP reserve_rateSEXP, SEXP reserve_interceptSEXP, SEXP half_rateSEXP, SEXP half_interceptSEXP, SEXP n_year_extrapSEXP, SEXP stop_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type well_id_vec(well_id_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type cum_prod_vec(cum_prod_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type date_start_vec(date_start_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type month_prod_vec(month_prod_vecSEXP);
    Rcpp::traits::input_parameter< double >::type year_last(year_lastSEXP);
    Rcpp::traits::input_parameter< double >::type reserve_rate(reserve_rateSEXP);
    Rcpp::traits::input_parameter< double >::type reserve_intercept(reserve_interceptSEXP);
    Rcpp::traits::input_parameter< double >::type half_rate(half_rateSEXP);
    Rcpp::traits::input_parameter< double >::type half_intercept(half_interceptSEXP);
    Rcpp::traits::input_parameter< double >::type n_year_extrap(n_year_extrapSEXP);
    Rcpp::traits::input_parameter< double >::type stop_threshold(stop_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(prod_extrap(well_id_vec, cum_prod_vec, date_start_vec, month_prod_vec, year_last, reserve_rate, reserve_intercept, half_rate, half_intercept, n_year_extrap, stop_threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hubbert_cumsum_by_tag", (DL_FUNC) &_hubbert_cumsum_by_tag, 2},
    {"_hubbert_cumsum_by_date", (DL_FUNC) &_hubbert_cumsum_by_date, 3},
    {"_hubbert_depletion_stats", (DL_FUNC) &_hubbert_depletion_stats, 3},
    {"_hubbert_lotka_volt", (DL_FUNC) &_hubbert_lotka_volt, 7},
    {"_hubbert_lotka_volt_cumulative", (DL_FUNC) &_hubbert_lotka_volt_cumulative, 8},
    {"_hubbert_prod_curve", (DL_FUNC) &_hubbert_prod_curve, 19},
    {"_hubbert_prod_extrap", (DL_FUNC) &_hubbert_prod_extrap, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_hubbert(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
