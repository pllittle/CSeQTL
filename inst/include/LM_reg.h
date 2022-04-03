#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _LMREG_
#define _LMREG_

void Rcpp_LM_calc_sigma2(const arma::uword& N,
	const arma::uword& P,const arma::vec& RES,
	double& sigma2);

void Rcpp_LM_residuals(const arma::vec& Y,
	const arma::mat& M,const arma::mat& EYE,
	arma::vec& RES);

Rcpp::List Rcpp_LM(const arma::mat& YY,
	const arma::mat& XX,const arma::uvec& var_groups,
	const arma::uword& test_type0,
	const int& ncores);


#endif