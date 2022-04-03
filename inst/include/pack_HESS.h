#include <RcppArmadillo.h>

#ifndef _PACK_HH_
#define _PACK_HH_

arma::mat Rcpp_CSeQTL_HESS(const arma::vec& TREC,
	const arma::vec& hap2,const arma::vec& ASREC,
	const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,
	const arma::vec& BETA,const double& PHI,
	const double& PSI,const arma::vec& KAPPA2Q,
	const arma::vec& ETA_TREC,const arma::vec& ALPHA,
	const arma::umat& GI,const arma::vec& upPARS,
	const arma::mat& I_np,const double& shift,
	const bool& show);

arma::mat Rcpp_CSeQTL_calc_HESS(const arma::vec& TREC,
	const arma::vec& hap2,const arma::vec& ASREC,
	const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,
	const arma::umat& GI,const arma::vec& PARS,
	const arma::vec& upPARS,const arma::mat& I_np,
	const double& shift);

void Rcpp_CSeQTL_hessBR(const arma::mat& hess,
	const arma::umat& GI,const arma::vec& upPARS,
	bool& rcond_nz,const bool& show);

#endif

