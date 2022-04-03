#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _MINOR_
#define _MINOR_

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}

double Rcpp_norm(const arma::vec& a);

arma::umat Rcpp_calc_GI(const arma::uword& PP,
	const arma::uword& QQ);

arma::mat Rcpp_CSeQTL_MU(const arma::umat& GI,
	const arma::vec& PARS);

double Rcpp_calc_MAF(const arma::vec& SNP,
	const bool& phasing,const bool& show);

arma::vec Rcpp_calc_MAF_all(const arma::mat& SNP,
	const bool& phasing,const bool& show,
	const int& ncores);

#endif

