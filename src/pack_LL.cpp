#include "../inst/include/pack_LL.h"

// --------------------
// Log Likelihood Functions

double Rcpp_CSeQTL_LL(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::vec& BETA,const double& PHI,const double& PSI,
	const arma::vec& KAPPA2Q,const arma::vec& ETA_TREC,const arma::vec& ALPHA){
	
	arma::uword ii, NN = SNP.n_elem, QQ = RHO.n_cols;
	arma::vec KAPPA = arma::ones<arma::vec>(QQ);
	if(QQ > 1) KAPPA.subvec(1,QQ-1) = KAPPA2Q;
	double LL = 0.0, log_mu_AA, log_mu, mu, xi_TREC,
		xi_ASE, rhoK, LL_TREC, LL_ASREC,
		iPHI = 0, iPHI_log_iPHI = 0, lg_iPHI = 0,
		pi, aa, bb;
	arma::vec ETA_KAP_TREC = ETA_TREC % KAPPA,
		ETA_KAP_ASE = (ETA_TREC % ALPHA) % KAPPA;
	
	if( PHI > 0.0 ){
		iPHI = 1.0 / PHI;
		iPHI_log_iPHI = iPHI * std::log(iPHI);
		lg_iPHI = std::lgamma(iPHI);
	}
	
	for(ii = 0; ii < NN; ii++){
		rhoK = arma::dot(RHO.row(ii).t(),KAPPA);
		log_mu_AA = arma::dot(BETA,XX.row(ii).t()) + std::log(rhoK);
		xi_TREC = arma::dot(RHO.row(ii).t(),ETA_KAP_TREC) / rhoK;
		if( SNP.at(ii) == 0 ){ // AA
			pi = 0.5;
			log_mu = log_mu_AA;
		} else if( SNP.at(ii) == 1 || SNP.at(ii) == 2 ){ // AB or BA
			log_mu = log_mu_AA + std::log(( 1.0 + xi_TREC ) / 2.0 );
			xi_ASE = arma::dot(RHO.row(ii).t(),ETA_KAP_ASE) / rhoK;
			pi = xi_ASE / (1.0 + xi_ASE);
		} else if( SNP.at(ii) == 3 ){ // BB
			log_mu = log_mu_AA + std::log(xi_TREC);
			pi = 0.5;
		} else {
			continue;
		}
		
		// TReC
		mu = std::exp(log_mu);
		if( PHI > 0.0 ){ // Negative Binomial
			LL_TREC = iPHI_log_iPHI - (TREC.at(ii) + iPHI) * std::log(iPHI + mu) - LGX1.at(ii);
			if( TREC.at(ii) > 0 ){
				LL_TREC += std::lgamma(TREC.at(ii) + iPHI) - lg_iPHI + TREC.at(ii) * log_mu;
			}
		} else { // Poisson
			LL_TREC = TREC.at(ii) * log_mu - mu - LGX1.at(ii);
		}
		
		// ASReC
		if( PSI > 0.0 ){ // Beta-binomial
			aa = pi / PSI; bb = (1.0 - pi) / PSI;
			
			if( PHASE.at(ii) == 1 ){
				if( SNP.at(ii) == 1 ){ 				// AB
					LL_ASREC = Rcpp_log_BB_2(hap2.at(ii),
						ASREC.at(ii),aa,bb);
				} else if( SNP.at(ii) == 2 ){ 	// BA
					LL_ASREC = Rcpp_log_BB_2(ASREC.at(ii) - hap2.at(ii),
						ASREC.at(ii),aa,bb);
				} else { 											// AA or BB
					LL_ASREC = Rcpp_log_BB_2(hap2.at(ii),
						ASREC.at(ii),aa,bb);
				}
				LL_ASREC += LBC.at(ii);
			} else {
				LL_ASREC = 0.0; // aka no phased reads
			}
		
		} else { // Binomial
			
			if( PHASE.at(ii) == 1 ){
				if( SNP.at(ii) == 1 ){ // AB
					LL_ASREC = Rcpp_log_binom(hap2.at(ii),ASREC.at(ii),pi);
				} else if( SNP.at(ii) == 2 ){ // BA
					LL_ASREC = Rcpp_log_binom(ASREC.at(ii) - hap2.at(ii),ASREC.at(ii),pi);
				} else { // AA or BB
					LL_ASREC = Rcpp_log_binom(hap2.at(ii),ASREC.at(ii),pi);
				}
				LL_ASREC += LBC.at(ii);
			} else {
				LL_ASREC = 0.0;
			}
			
		}
		
		// TReC + ASReC
		LL += LL_TREC + LL_ASREC;
	}
	
	return LL;
}

double Rcpp_CSeQTL_calc_LL(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,const arma::vec& PARS){
	
	return Rcpp_CSeQTL_LL(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
		PARS.subvec(GI.at(0,0),GI.at(0,1)),
		std::exp(PARS.at(GI.at(1,0))),
		std::exp(PARS.at(GI.at(4,0))),
		arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1))),
		arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))),
		arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1))));
	
}

