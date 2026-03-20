#include <Rcpp.h>
#include "robust_core.h"
#include "robscale_config.h"

// [[Rcpp::export]]
Rcpp::NumericVector cpp_get_are_constants() {
  Rcpp::NumericVector are = Rcpp::NumericVector::create(
    Rcpp::_["sd_c4"]      = robscale::ARE_SD_C4,
    Rcpp::_["gmd"]        = robscale::ARE_GMD,
    Rcpp::_["mad_scaled"] = robscale::ARE_MAD,
    Rcpp::_["iqr_scaled"] = robscale::ARE_IQR,
    Rcpp::_["sn"]         = robscale::ARE_SN,
    Rcpp::_["qn"]         = robscale::ARE_QN,
    Rcpp::_["robScale"]   = robscale::ARE_ROBSCALE
  );
  return are;
}
