#include "../inst/include/poisson.h"

// --------------------
// Poisson related functions

double Rcpp_POIS_reg_LL(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA,const arma::vec& LGX1){
	
	arma::uword ii, N = Y.n_elem;
	double LL = 0.0, lambda_i;
	
	for(ii = 0; ii < N; ii++){
		lambda_i = std::exp(arma::dot(X.row(ii).t(),BETA) + O.at(ii));
		LL += -1.0 * lambda_i + Y.at(ii) * std::log(lambda_i) - LGX1.at(ii);
	}
	
	return LL;
}

arma::vec Rcpp_POIS_reg_GRAD(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA){
	
	arma::uword ii, N = Y.n_elem;
	arma::vec GRAD = arma::zeros<arma::vec>(BETA.n_elem);
	double lambda_i;
	
	for(ii = 0; ii < N; ii++){
		lambda_i = std::exp(arma::dot(X.row(ii).t(),BETA) + O.at(ii));
		GRAD += (Y.at(ii) - lambda_i) * X.row(ii).t();
	}
	
	return GRAD;
}

arma::mat Rcpp_POIS_reg_HESS(const arma::vec& Y,
	const arma::mat& X,const arma::vec& O,
	const arma::vec& BETA){
	
	arma::uword ii, N = Y.n_elem;
	arma::mat HESS = arma::zeros<arma::mat>(BETA.n_elem,BETA.n_elem);
	double lambda_i;
	
	for(ii = 0; ii < N; ii++){
		lambda_i = std::exp(arma::dot(X.row(ii).t(),BETA) + O.at(ii));
		HESS += -1.0 * lambda_i * X.row(ii).t() * X.row(ii);
	}
	
	return HESS;
}

void Rcpp_POIS_reg(const arma::vec& Y,const arma::mat& X,
	const arma::vec& O,arma::vec& BETA,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-7,const bool& show = true){
	
	arma::uword iter = 0, jj, uu, num_params = X.n_cols, NN = X.n_rows;
	double curr_LL = 0.0, old_LL, new_LL;
	arma::vec curr_BETA = BETA, new_BETA = BETA,
		old_grad = arma::zeros<arma::vec>(num_params),
		old_hess_grad = old_grad,
		LGX1 = arma::zeros<arma::vec>(NN);
	for(jj = 0; jj < NN; jj++){
		LGX1.at(jj) = std::lgamma( Y.at(jj) + 1.0 );
	}
	
	while(iter < max_iter){
		old_LL = Rcpp_POIS_reg_LL(Y,X,O,BETA,LGX1);
		old_grad = Rcpp_POIS_reg_GRAD(Y,X,O,BETA);
		old_hess_grad = -1.0 * arma::inv(Rcpp_POIS_reg_HESS(Y,X,O,BETA)) * old_grad;
		old_grad /= std::max(1.0,Rcpp_norm(old_grad));
		old_hess_grad /= std::max(1.0,Rcpp_norm(old_hess_grad));
		uu = 0;
		for(jj = 0; jj <= 30; jj++){
			new_BETA = BETA + old_hess_grad / std::pow(4.0,jj);
			new_LL = Rcpp_POIS_reg_LL(Y,X,O,new_BETA,LGX1);
			if( new_LL > old_LL ){
				BETA = new_BETA;
				old_LL = new_LL;
				uu = 1;
				break;
			} else {
				new_BETA = BETA + old_grad / std::pow(4.0,jj);
				new_LL = Rcpp_POIS_reg_LL(Y,X,O,new_BETA,LGX1);
				if( new_LL > old_LL ){
					BETA = new_BETA;
					old_LL = new_LL;
					uu = 2;
					break;
				}
			}
		}
		
		if( show ){
			// if(uu == 0) printR_obj("\tNo more update");
			if(uu == 0) Rcpp::Rcout << "\tNo more update\n";
		}
		
		if(uu == 0) break;
		
		if(iter > 0){
			if( std::abs(old_LL - curr_LL) < eps && Rcpp_norm(BETA - curr_BETA) < eps ){
				old_grad = Rcpp_POIS_reg_GRAD(Y,X,O,BETA);
				old_hess_grad = -1.0 * arma::inv(Rcpp_POIS_reg_HESS(Y,X,O,BETA)) * old_grad;
				if( Rcpp_norm(old_grad) < eps && Rcpp_norm(old_hess_grad) < eps ){
					break;
				}
			}
		}
		
		curr_LL = old_LL;
		curr_BETA = BETA;
		iter++;
	}
	
	if( show ){
		old_LL = Rcpp_POIS_reg_LL(Y,X,O,BETA,LGX1);
		old_grad = Rcpp_POIS_reg_GRAD(Y,X,O,BETA);
		arma::mat hess = Rcpp_POIS_reg_HESS(Y,X,O,BETA), covar = hess;
			covar = arma::inv(-1.0 * hess);
		Rcpp::Rcout << "\tIter = " << iter
			<< "; LL = " << old_LL << "\n";
		Rcpp::Rcout << "\tNorm_Grad = " << Rcpp_norm(old_grad) 
			<< "; Norm_iHG = " << Rcpp_norm(covar * old_grad) << "\n";
		Rcpp::Rcout << "\tBeta = " << BETA.t();
		Rcpp::Rcout << "\tVar = " << arma::diagvec(covar).t();
	}
	
}

