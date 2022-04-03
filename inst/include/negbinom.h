#include <RcppArmadillo.h>
#include "minor.h"
#include "poisson.h"
#include "binom_beta.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _NEGBINOM_
#define _NEGBINOM_

double Rcpp_log_NB(const double& yy,const double& mu,
	const double& phi,const double& lgy1);

double Rcpp_NB_reg_LL(const arma::vec& YY,const arma::mat& XX,
	const arma::vec& OO,const arma::vec& LGY1,const arma::vec& BETA,
	const double& PHI);

arma::vec Rcpp_NB_reg_GRAD(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& BETA,const double& PHI);

arma::mat Rcpp_NB_reg_HESS(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& BETA,const double& PHI);

void Rcpp_NB_reg_NR(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& LGY1,arma::vec& PARS,
	const arma::uword& max_iter,
	const double& eps,const bool& show);

void Rcpp_NB_iPARS(const arma::vec& iBETA,
	const double& iPHI,const arma::umat& GI,
	arma::vec& PARS);

arma::vec Rcpp_NB_reg_one(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::uword& max_iter,
	const double& eps,const bool& show);

void Rcpp_CSeQTL_BETA_PHI(const arma::mat& XX,
	const arma::mat& TREC,arma::mat& iBETA,
	arma::vec& iPHI,const bool& show,
	const int& ncores);

#endif

