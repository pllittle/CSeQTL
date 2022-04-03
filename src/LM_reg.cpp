#include "../inst/include/LM_reg.h"

// --------------------
// Linear Regression Functions

void Rcpp_LM_calc_sigma2(const arma::uword& N,const arma::uword& P,
	const arma::vec& RES,double& sigma2){
	sigma2 = arma::dot(RES,RES) / (1.0 * (N - P) );
}

void Rcpp_LM_residuals(const arma::vec& Y,const arma::mat& M,
	const arma::mat& EYE,arma::vec& RES){
	RES = (EYE - M) * Y;
}

// [[Rcpp::export]]
Rcpp::List Rcpp_LM(const arma::mat& YY,const arma::mat& XX,
	const arma::uvec& var_groups,const arma::uword& test_type0 = 1,
	const int& ncores = 1){
	
	arma::uword test_type = test_type0;
	if( test_type != 1 ){
		test_type = 3;
	}
	
	arma::uvec uniq_vars = arma::unique(var_groups);
	arma::uword N = XX.n_rows, K = YY.n_rows, mm, 
		P = uniq_vars.n_elem - 1;
	arma::mat RESID = arma::zeros<arma::mat>(K,N);
	arma::mat PVALUE = arma::zeros<arma::mat>(K,P);
	arma::mat EYE = arma::eye<arma::mat>(N,N);
	
	arma::mat MM_full = XX * arma::inv(XX.t() * XX) * XX.t();
	Rcpp::Rcout << "Store differences in projection matrices ...\n";
	arma::cube MM_diff = arma::zeros<arma::cube>(N,N,P);
	arma::uvec var_df = arma::zeros<arma::uvec>(P);
	for(mm = 0; mm < P; mm++){
		if( test_type == 1 ){
			// For Type I tests: sequencial order of covariates
			arma::uvec sub_cols = arma::find( var_groups <= mm );
			// Calculate projection matrix
			MM_diff.slice(mm) = XX.cols(sub_cols) * arma::inv(XX.cols(sub_cols).t()
				* XX.cols(sub_cols)) * XX.cols(sub_cols).t();
			if( mm > 0 ){
				MM_diff.slice(mm - 1) = MM_diff.slice(mm) - MM_diff.slice(mm - 1);
				if( mm == P - 1 ) MM_diff.slice(mm) = MM_full - MM_diff.slice(mm);
			}
		} else {
			// For Type III tests: added last tests
			arma::uvec sub_cols = arma::find(var_groups != mm + 1);
			// Calculate projection matrix
			MM_diff.slice(mm) = MM_full - ( XX.cols(sub_cols) * arma::inv(XX.cols(sub_cols).t()
				* XX.cols(sub_cols)) * XX.cols(sub_cols).t() );
		}
		
		arma::uvec index_cols = arma::find( var_groups == mm + 1 );
		var_df.at(mm) = index_cols.n_elem;
	}
	
	Rcpp::Rcout << "Loop thru all " << K 
		<< " outcomes to calculate residuals and Type " 
		<< test_type << " pvalues ...\n";
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(ncores,YY,MM_full,EYE,N,RESID,K,P, \
			XX,MM_diff,var_df,PVALUE)
	#endif
	for(arma::uword kk = 0; kk < K; kk++){
		if( ncores == 1 ){
			if( (kk+1) % 100 == 0 ) Rcpp::Rcout << ".";
			if( (kk+1) % 5000 == 0 ) Rcpp::Rcout << kk + 1 << "\n";
		}
		
		double sigma2_full = 0;
		arma::vec tmp_resid = arma::zeros<arma::vec>(N);
		
		// Calculate full model: residuals, sigma2, likelihood
		Rcpp_LM_residuals(YY.row(kk).t(),MM_full,EYE,tmp_resid);
		RESID.row(kk) = tmp_resid.t();
		Rcpp_LM_calc_sigma2(N,XX.n_cols,RESID.row(kk).t(),sigma2_full);
		
		// Calculate reduced models stats
		for(arma::uword mm2 = 0; mm2 < P; mm2++){
			// difference in sum_square_error Y' (M_Larger - M_Smaller) Y 
			double ss_red = arma::dot(YY.row(kk).t(),MM_diff.slice(mm2) * YY.row(kk).t())
				/ var_df.at(mm2);
			double F_stat = ss_red / sigma2_full;
			PVALUE.at(kk,mm2) = 1.0 - R::pf(F_stat,var_df.at(mm2),N - XX.n_cols,1,0);
		}
	
	}
	
	Rcpp::Rcout << "\n";
	return Rcpp::List::create(
    Rcpp::Named("RESID",RESID),
    Rcpp::Named("PVALUE",PVALUE));
}

