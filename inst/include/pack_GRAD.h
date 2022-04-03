#include <RcppArmadillo.h>

#ifndef _PACK_GG_
#define _PACK_GG_

arma::vec Rcpp_CSeQTL_GRAD(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& BETA,
	const double& PHI,const double& PSI,const arma::vec& KAPPA2Q,
	const arma::vec& ETA_TREC,const arma::vec& ALPHA,const arma::umat& GI,
	const arma::vec& upPARS);

arma::vec Rcpp_CSeQTL_calc_GRAD(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::umat& GI,
	const arma::vec& upPARS,const arma::vec& PARS);

void Rcpp_CSeQTL_control_PAR(arma::vec& PARS,const arma::uword& QQ,
	const arma::umat& GI,const arma::vec& upPARS);

#endif
