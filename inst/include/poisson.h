#include <RcppArmadillo.h>
#include "minor.h"

#ifndef _POIS_
#define _POIS_

double Rcpp_POIS_reg_LL(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA,const arma::vec& LGX1);

arma::vec Rcpp_POIS_reg_GRAD(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA);

arma::mat Rcpp_POIS_reg_HESS(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA);

void Rcpp_POIS_reg(const arma::vec& Y,const arma::mat& X,
	const arma::vec& O,arma::vec& BETA,
	const arma::uword& max_iter,
	const double& eps,const bool& show);

#endif
