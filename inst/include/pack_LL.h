#include <RcppArmadillo.h>
#include "binom_beta.h"

#ifndef _PACK_LL_
#define _PACK_LL_

double Rcpp_CSeQTL_LL(const arma::vec& TREC,
	const arma::vec& LGX1,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,
	const arma::vec& BETA,const double& PHI,
	const double& PSI,const arma::vec& KAPPA2Q,
	const arma::vec& ETA_TREC,const arma::vec& ALPHA);

double Rcpp_CSeQTL_calc_LL(const arma::vec& TREC,
	const arma::vec& LGX1,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,
	const arma::umat& GI,const arma::vec& PARS);

#endif

