#include "../inst/include/binom_beta.h"

// --------------------
// Binomial and Beta-Binomial Log Likelihood

double Rcpp_log_binom(const double& xx,const double &nn,
	const double& pi){
	return xx * std::log(pi) + (nn - xx) * std::log(1.0 - pi);
}

double Rcpp_log_BB_2(const double& xx,const double& nn,
	const double& aa,const double& bb){
	
	return std::lgamma(xx + aa) +
		std::lgamma(nn - xx + bb) +
		std::lgamma(aa + bb) -
		std::lgamma(nn + aa + bb) -
		std::lgamma(aa) - std::lgamma(bb);
}

