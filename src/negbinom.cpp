#include "../inst/include/negbinom.h"

// --------------------
// Negative Binomial related functions

double Rcpp_log_NB(const double& yy,const double& mu,
	const double& phi,const double& lgy1){
	
	double loglike = 0.0, nn = 1.0 / phi;
	// lgy1 = std::lgamma(yy + 1.0)
	
	if( yy > 0 ){
		loglike += std::lgamma(yy + nn) -
			std::lgamma(nn) - lgy1;
		loglike += yy * std::log(mu);
	}
	
	loglike += nn * std::log(nn) -
		(nn + yy) * std::log(nn + mu);
	
	return loglike;
}

double Rcpp_NB_reg_LL(const arma::vec& YY,const arma::mat& XX,
	const arma::vec& OO,const arma::vec& LGY1,const arma::vec& BETA,
	const double& PHI){
	
	arma::uword ii;
	double mu, LL = 0.0;
	
	for(ii = 0; ii < YY.n_elem; ii++){
		mu = std::exp( arma::dot(XX.row(ii).t(),BETA) + OO.at(ii) );
		LL += Rcpp_log_NB(YY.at(ii),mu,PHI,LGY1.at(ii));
	}
	
	return LL;
}

arma::vec Rcpp_NB_reg_GRAD(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& BETA,const double& PHI){
	
	arma::uword ii, pp = XX.n_cols;
	double mu, nn = 1.0 / PHI, digam_nn = R::digamma(nn),
		log_nn = std::log(nn);
	arma::vec GRAD = arma::zeros<arma::vec>(pp + 1);
	
	for(ii = 0; ii < YY.n_elem; ii++){
		mu = std::exp( arma::dot(XX.row(ii).t(),BETA) + OO.at(ii) );
		// Part BETA
		GRAD.subvec(0,pp - 1) += (YY.at(ii) / mu - 
			(nn + YY.at(ii)) / (nn + mu)) * mu * XX.row(ii).t();
		// Part PHI*
		GRAD.at(pp) += (R::digamma(YY.at(ii) + nn) - digam_nn + 
			1.0 + log_nn - (nn + YY.at(ii)) / (nn + mu) - 
			std::log(nn + mu)) * (-1.0 * nn);
	}

	return GRAD;
}

arma::mat Rcpp_NB_reg_HESS(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& BETA,const double& PHI){
	
	arma::uword ii, pp = XX.n_cols;
	double mu, nn = 1.0 / PHI, trigam_nn = R::trigamma(nn),
		digam_nn = R::digamma(nn), log_nn = std::log(nn);
	arma::mat HESS = arma::zeros<arma::mat>(pp + 1,pp + 1);
	arma::vec hess_beta_phi = arma::zeros<arma::vec>(pp);
	
	for(ii = 0; ii < YY.n_elem; ii++){
		mu = std::exp( arma::dot(XX.row(ii).t(),BETA) + OO.at(ii) );
		// Part2 BETA
		HESS.submat(0,0,pp - 1,pp - 1) += 
			-1.0 * (nn + YY.at(ii)) / std::pow(nn + mu,2.0) * 
			nn * mu * XX.row(ii).t() * XX.row(ii);
		
		// Part2 BETA*PHI*
		hess_beta_phi = -1.0 * (YY.at(ii) - mu) / 
			std::pow(nn + mu,2.0) * nn * mu * XX.row(ii).t();
		HESS(arma::span(0,pp - 1),pp) += hess_beta_phi;
		HESS(pp,arma::span(0,pp - 1)) += hess_beta_phi.t();
		
		// Part2 PHI*
		HESS.at(pp,pp) += nn * ( nn * (R::trigamma(YY.at(ii) + nn) - trigam_nn + PHI +
			(YY.at(ii) - mu)/std::pow(nn + mu,2.0) - 1.0/(nn + mu) ) + 
			R::digamma(YY.at(ii) + nn) - digam_nn + 1.0 + 
			log_nn - (nn + YY.at(ii)) / (nn + mu) - std::log(nn + mu)
			);
	}

	return HESS;
}

void Rcpp_NB_reg_NR(const arma::vec& YY,
	const arma::mat& XX,const arma::vec& OO,
	const arma::vec& LGY1,arma::vec& PARS,
	const arma::uword& max_iter = 4e3,
	const double& eps = 1e-7,const bool& show = true){
	
	arma::uword iter = 0, jj, uu, pp = XX.n_cols;
	double curr_LL = 0.0, old_LL, new_LL;
	arma::vec curr_PARS = PARS, new_PARS = PARS,
		GRAD = arma::zeros<arma::vec>(pp + 1),
		iHESS_GRAD = GRAD;

	while(iter < max_iter){
		old_LL = Rcpp_NB_reg_LL(YY,XX,OO,LGY1,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp)));
		GRAD = Rcpp_NB_reg_GRAD(YY,XX,OO,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp)));
		iHESS_GRAD = -1.0 * arma::inv(Rcpp_NB_reg_HESS(YY,XX,OO,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp)))) * GRAD;
		GRAD /= std::max(1.0,Rcpp_norm(GRAD));
		iHESS_GRAD /= std::max(1.0,Rcpp_norm(iHESS_GRAD));
		uu = 0;
		for(jj = 0; jj <= 30; jj++){
			new_PARS = PARS + iHESS_GRAD / std::pow(4.0,jj);
			new_LL = Rcpp_NB_reg_LL(YY,XX,OO,LGY1,
				new_PARS.subvec(0,pp-1),std::exp(new_PARS.at(pp)));
			if( new_LL > old_LL ){
				PARS = new_PARS;
				old_LL = new_LL;
				uu = 1;
				break;
			} else {
				new_PARS = PARS + GRAD / std::pow(4.0,jj);
				new_LL = Rcpp_NB_reg_LL(YY,XX,OO,LGY1,
					new_PARS.subvec(0,pp-1),std::exp(new_PARS.at(pp)));
				if( new_LL > old_LL ){
					PARS = new_PARS;
					old_LL = new_LL;
					uu = 2;
					break;
				}
			}
		}
		
		if(show){
			if(uu == 0){
				printR_obj("\tNo more update");
			}
		}
		
		if( uu == 0 ) break;
		
		if(iter > 0){
			if( std::abs(curr_LL - old_LL) < eps && Rcpp_norm(curr_PARS - PARS) < eps ){
				GRAD = Rcpp_NB_reg_GRAD(YY,XX,OO,
					PARS.subvec(0,pp-1),std::exp(PARS.at(pp)));
				iHESS_GRAD = -1.0 * arma::inv(Rcpp_NB_reg_HESS(YY,XX,OO,
					PARS.subvec(0,pp-1),std::exp(PARS.at(pp)))) * GRAD;
				if( Rcpp_norm(GRAD) < eps && Rcpp_norm(iHESS_GRAD) < eps ){
					break;
				}
			}
		}
		
		curr_PARS = PARS;
		curr_LL = old_LL;
		iter++;
	}
	
	if( show ){
		old_LL = Rcpp_NB_reg_LL(YY,XX,OO,LGY1,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp)));
		GRAD = Rcpp_NB_reg_GRAD(YY,XX,OO,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp)));
		arma::mat hess = Rcpp_NB_reg_HESS(YY,XX,OO,
			PARS.subvec(0,pp-1),std::exp(PARS.at(pp))),
			covar = hess;
		covar = arma::inv(-1.0 * hess);
		Rcpp::Rcout << "\tIter = " << iter 
			<< "; LL = " << old_LL << "\n";
		Rcpp::Rcout << "\tNorm_Grad = " << Rcpp_norm(GRAD) 
			<< "; Norm_iHG = " << Rcpp_norm(covar * GRAD) << "\n";
		Rcpp::Rcout << "\tBeta = " << PARS.subvec(0,pp-1).t();
		Rcpp::Rcout << "\tPhi = " << std::exp(PARS.at(pp)) << "\n";
		Rcpp::Rcout << "\tVar = " << arma::diagvec(covar).t();
	}
	
}

void Rcpp_NB_iPARS(const arma::vec& iBETA,const double& iPHI,
	const arma::umat& GI,arma::vec& PARS){
	
	PARS.subvec(GI.at(0,0),GI.at(0,1)) = iBETA;
	PARS.at(GI.at(1,0)) = std::log(iPHI);
}


// [[Rcpp::export]]
arma::vec Rcpp_NB_reg_one(const arma::vec& YY,const arma::mat& XX,
	const arma::vec& OO,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-5,const bool& show = true){
	
	arma::uword PP = XX.n_cols;
	arma::vec LGY1 = arma::lgamma(YY + 1.0),
		PARS = arma::zeros<arma::vec>(PP + 1),
		tmp_BETA = arma::zeros<arma::vec>(PP);
	
	if(show) Rcpp::Rcout << "### Run Poisson regression ...\n";
	if( PP == 1 ){
		PARS.at(0) = std::log(arma::mean(YY));
	} else {
		// Poisson Reg to initialize beta
		Rcpp_POIS_reg(YY,XX,OO,tmp_BETA,4e3,1e-7,show);
		PARS.subvec(0,PP-1) = tmp_BETA;
	}
	
	// Negative Binomial Regression
	if(show) Rcpp::Rcout << "### Run Negative Binomial regression ...\n";
	Rcpp_NB_reg_NR(YY,XX,OO,LGY1,PARS,4e3,1e-7,show);
	
	return PARS;
}

void Rcpp_CSeQTL_BETA_PHI(const arma::mat& XX,const arma::mat& TREC,
	arma::mat& iBETA,arma::vec& iPHI,const bool& show = true,
	const int& ncores = 1){
		
	// This function will run regular negative binomial regression, 
	//	with offset = 0, aka ignoring snps/kappa/eta
	
	arma::uword PP = XX.n_cols, thres = 0,
		NN = XX.n_rows, GG = TREC.n_rows;
	arma::vec OO = arma::zeros<arma::vec>(NN);
	
	if( show ){
		if(PP == 1){
			Rcpp::Rcout << "Storing intercept-only NBreg BETA and PHI parameters";
		} else {
			Rcpp::Rcout << "Storing initial NBreg BETA and PHI parameters";
		}
	}
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(ncores,PP,GG,OO,iBETA,iPHI,TREC,XX)
	#endif
	for(arma::uword gg = 0; gg < GG; gg++){
		if( (ncores == 1) && show && (100.00 * (gg + 1.0) / GG > thres) ){
			Rcpp::Rcout << ".";
			thres += 10;
		}
		
		// Initialize NB parameters
		arma::vec PARS = Rcpp_NB_reg_one(TREC.row(gg).t(),XX,
			OO,4e3,1e-7,false);
		
		// Store estimates
		iBETA.row(gg) = PARS.subvec(0,PP-1).t();
		iPHI.at(gg) = std::exp(PARS.at(PP));
	}
	
	if( show ){
		if( ncores == 1 ){
			Rcpp::Rcout << ";\n";
		} else {
			Rcpp::Rcout << "..........;\n";
		}
	}
	
}


