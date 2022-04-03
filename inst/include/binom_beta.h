#include <RcppArmadillo.h>

#ifndef _BINBETA_
#define _BINBETA_

double Rcpp_log_binom(const double& xx,
	const double &nn,const double& pi);

double Rcpp_log_BB_2(const double& xx,
	const double& nn,const double& aa,
	const double& bb);

#endif


