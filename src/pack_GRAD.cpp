#include "../inst/include/pack_GRAD.h"

// --------------------
// Gradient Functions

arma::vec Rcpp_CSeQTL_GRAD(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& BETA,
	const double& PHI,const double& PSI,const arma::vec& KAPPA2Q,
	const arma::vec& ETA_TREC,const arma::vec& ALPHA,const arma::umat& GI,
	const arma::vec& upPARS){
	
	// GRAD = (BETA,log_PHI,log_KAPPA2Q,log_ETA,log_PSI,ALPHA)
	arma::uword ii, jj, N = SNP.n_elem, QQ = RHO.n_cols;
	// check_ETA, check_ALPHA are meant to check if the gradient should be calculated
	bool check_KAPPA = true,check_ETA = arma::any(upPARS.subvec(GI.at(3,0),GI.at(3,1)) == 1.0),
		check_ALPHA = arma::any(upPARS.subvec(GI.at(5,0),GI.at(5,1)) == 1.0);
	if( QQ > 1 ){
		check_KAPPA = arma::any(upPARS.subvec(GI.at(2,0),GI.at(2,1)) == 1.0);
	}
	arma::vec KAPPA = arma::ones<arma::vec>(QQ);
	arma::mat jacob_k_log_k2q = arma::zeros<arma::mat>(QQ,QQ),
		jacob_alpha = arma::diagmat(ALPHA);
	if( QQ > 1 ){
		KAPPA.subvec(1,QQ-1) = KAPPA2Q;
		for(jj = 0; jj < QQ - 1; jj++){
			jacob_k_log_k2q.at(jj+1,jj) = KAPPA2Q.at(jj);
		}
	}
	arma::vec GRAD = arma::zeros<arma::vec>(GI.at(5,1)+1),
		ETA_ASE = ETA_TREC % ALPHA,ETA_KAP_TREC = ETA_TREC % KAPPA,
		ETA_KAP_ASE = ETA_ASE % KAPPA,ETA_KAP_RHO_TREC = arma::zeros<arma::vec>(QQ),
		ETA_KAP_RHO_ASE = ETA_KAP_RHO_TREC;
	
	double log_mu_AA, log_mu, mu, aa = 0.0, bb = 0.0,
		dig_aa = 0.0, dig_bb = 0.0,xi_TREC, xi_ASE, pi,
		rhoK, iPHI = 0.0,dig_iPHI = 0.0, log_iPHI = 0.0, 
		iPSI = 0.0, dig_iPSI = 0.0, part_li_mu_i,
		part_offset_xi_TREC, dig_hB_a = 0.0,
		dig_hA_b = 0.0, dig_hA_a = 0.0, dig_hB_b = 0.0,
		part_pi_xiA_rK, part_li2_pi = 0.0;
	arma::vec part_mu_kap = arma::zeros<arma::vec>(QQ),
		eta_xi_TREC = part_mu_kap,eta_xi_ASE = part_mu_kap,
		part_pi_kap_star = part_mu_kap,part_pi_eta_star = part_mu_kap,
		part_pi_alpha_star = part_mu_kap;
	if( PHI > 0.0 ){
		iPHI = 1.0 / PHI;
		dig_iPHI = R::digamma(iPHI);
		log_iPHI = std::log(iPHI);
	}
	if( PSI > 0.0 ){
		iPSI = 1.0 / PSI;
		dig_iPSI = R::digamma(iPSI);
	}
	
	for(ii = 0; ii < N; ii++){
		rhoK = arma::dot(RHO.row(ii).t(),KAPPA);
		log_mu_AA = arma::dot(BETA,XX.row(ii).t()) + std::log( rhoK );
		ETA_KAP_RHO_TREC = ETA_KAP_TREC % RHO.row(ii).t();
		xi_TREC = arma::sum(ETA_KAP_RHO_TREC) / rhoK;
		ETA_KAP_RHO_ASE = ETA_KAP_ASE % RHO.row(ii).t();
		xi_ASE = arma::sum(ETA_KAP_RHO_ASE) / rhoK;
		
		if( SNP.at(ii) == 0 ){
			// Genotype AA
			log_mu = log_mu_AA;
			pi = 0.5;
		} else if( SNP.at(ii) == 1 || SNP.at(ii) == 2 ){
			// Genotype AB or BA
			log_mu = log_mu_AA + std::log(( 1.0 + xi_TREC ) / 2.0 );
			pi = xi_ASE / (1.0 + xi_ASE);
		} else if( SNP.at(ii) == 3 ){
			// Genotype BB
			log_mu = log_mu_AA + std::log(xi_TREC);
			pi = 0.5;
		} else {
			continue;
		}
		mu = std::exp(log_mu);
		eta_xi_TREC = ETA_TREC - xi_TREC;
		eta_xi_ASE = ETA_ASE - xi_ASE;
		part_pi_xiA_rK = std::pow(1.0 - pi,2.0) / rhoK;
		part_li2_pi = 0.0;
		if( PHASE.at(ii) == 1 ){
			if( PSI > 0.0 ){
				aa = pi / PSI; bb = (1.0 - pi) / PSI;
				dig_aa = R::digamma(aa);
				dig_bb = R::digamma(bb);
				dig_hB_a = R::digamma(hap2.at(ii) + aa);
				dig_hA_b = R::digamma(ASREC.at(ii) - hap2.at(ii) + bb);
				dig_hA_a = R::digamma(ASREC.at(ii) - hap2.at(ii) + aa);
				dig_hB_b = R::digamma(hap2.at(ii) + bb);
				if( SNP.at(ii) == 1 ){
					part_li2_pi = iPSI * ( dig_hB_a - dig_hA_b - dig_aa + dig_bb );
				} else if( SNP.at(ii) == 2 ){
					part_li2_pi = iPSI * ( dig_hA_a - dig_hB_b - dig_aa + dig_bb );
				}
			} else {
				if( SNP.at(ii) == 1 ){
					part_li2_pi = hap2.at(ii) / pi - (ASREC.at(ii) - hap2.at(ii)) / (1.0 - pi);
				} else if( SNP.at(ii) == 2 ){
					part_li2_pi = (ASREC.at(ii) - hap2.at(ii)) / pi - hap2.at(ii) / (1.0 - pi);
				}
			}
		}
		
		part_offset_xi_TREC = ( 1.0*(SNP.at(ii) == 1 || SNP.at(ii) == 2) 
			/ (1.0 + xi_TREC) + 1.0*(SNP.at(ii) == 3) / xi_TREC );
		
		// part_beta
		if( PHI > 0.0 ){ // Negative Binomial
			part_li_mu_i = TREC.at(ii) / mu - (TREC.at(ii) + iPHI) / (mu + iPHI);
		} else { // Poisson
			part_li_mu_i = TREC.at(ii) / mu - 1.0;
		}
		GRAD.subvec(GI.at(0,0),GI.at(0,1)) += part_li_mu_i * mu * XX.row(ii).t();
		
		// part_log_phi
		if( PHI > 0.0 ){
			GRAD.at(GI.at(1,0)) += 
				-1.0 * iPHI * ( R::digamma(TREC.at(ii) + iPHI) - dig_iPHI + 
				1.0 + log_iPHI - (TREC.at(ii) + iPHI) / (mu + iPHI) - 
				std::log(mu + iPHI) );
		}
		
		// part_log_kappa2Q
		if( QQ > 1 && check_KAPPA ){
			
			part_mu_kap = mu * (1.0 + part_offset_xi_TREC * eta_xi_TREC) % RHO.row(ii).t() / rhoK;
			part_pi_kap_star = part_pi_xiA_rK * jacob_k_log_k2q.cols(0,QQ-2).t() * (eta_xi_ASE % RHO.row(ii).t());
			
			// TReC piece
			GRAD.subvec(GI.at(2,0),GI.at(2,1)) +=
				jacob_k_log_k2q.cols(0,QQ-2).t() * part_li_mu_i * part_mu_kap;
			
			// ASReC piece
			if( PHASE.at(ii) == 1 ){
				GRAD.subvec(GI.at(2,0),GI.at(2,1)) += part_li2_pi * part_pi_kap_star;
			}
		
		}
		
		// part_log_eta
		if( check_ETA ){
			// part_pi_eta_star = part_pi/part_eta * part_eta/part_eta_star
			part_pi_eta_star = part_pi_xiA_rK * ETA_KAP_RHO_ASE;
			
			// TReC piece
			GRAD.subvec(GI.at(3,0),GI.at(3,1)) += part_li_mu_i * 
				part_offset_xi_TREC * mu * ETA_KAP_RHO_TREC / rhoK;

			// ASReC piece
			if( PHASE.at(ii) == 1 ){
				GRAD.subvec(GI.at(3,0),GI.at(3,1)) += part_li2_pi * part_pi_eta_star;
			}
		
		}
		
		// part_log_PSI
		if( PSI > 0.0 && PHASE.at(ii) == 1 ){
			if( SNP.at(ii) == 2 ){ // Genotype BA
				GRAD.at(GI.at(4,0)) += 
					-iPSI * (pi * (dig_hA_a - dig_aa) + (1.0 - pi) * (dig_hB_b - dig_bb) +
					dig_iPSI - R::digamma(ASREC.at(ii) + iPSI));
			} else { // Genotypes AA, AB, BB
				GRAD.at(GI.at(4,0)) += 
					-iPSI * (pi * (dig_hB_a - dig_aa) + (1.0 - pi) * (dig_hA_b - dig_bb) +
					dig_iPSI - R::digamma(ASREC.at(ii) + iPSI));
			}
		}
		
		// part_ALPHA
		if( check_ALPHA && PHASE.at(ii) == 1 ){
			part_pi_alpha_star = part_pi_xiA_rK * ETA_KAP_RHO_ASE;
			GRAD.subvec(GI.at(5,0),GI.at(5,1)) += part_li2_pi * part_pi_alpha_star;
		}
		
	}
	
	// Forcing gradient terms to 0, by checking if and upPARAM elements are not 1.0
	if( upPARS.at(0) == 0.0 ){
		GRAD.at(0) = 0.0;
	}
	if( QQ > 1 ){
		if( arma::any(upPARS.subvec(GI.at(2,0),GI.at(2,1)) == 0.0) ){
			GRAD.subvec(GI.at(2,0),GI.at(2,1)) %= upPARS.subvec(GI.at(2,0),GI.at(2,1));
		}
	}
	if( arma::any(upPARS.subvec(GI.at(3,0),GI.at(3,1)) == 0.0) ){
		GRAD.subvec(GI.at(3,0),GI.at(3,1)) %= upPARS.subvec(GI.at(3,0),GI.at(3,1));
	}
	if( arma::any(upPARS.subvec(GI.at(5,0),GI.at(5,1)) == 0.0) ){
		GRAD.subvec(GI.at(5,0),GI.at(5,1)) %= upPARS.subvec(GI.at(5,0),GI.at(5,1));
	}
	
	return GRAD;	
}

arma::vec Rcpp_CSeQTL_calc_GRAD(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::umat& GI,
	const arma::vec& upPARS,const arma::vec& PARS){
	
	return Rcpp_CSeQTL_GRAD(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,
		PARS.subvec(GI.at(0,0),GI.at(0,1)),std::exp(PARS.at(GI.at(1,0))),
		std::exp(PARS.at(GI.at(4,0))),arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1))),
		arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))),
		arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1))),
		GI,upPARS);
}

void Rcpp_CSeQTL_control_PAR(arma::vec& PARS,const arma::uword& QQ,
	const arma::umat& GI,const arma::vec& upPARS){
	
	arma::uword ct;
	double neg_inf = -arma::datum::inf, small_log_mu = -3.0;
	arma::vec upKAPPA = arma::zeros<arma::vec>(QQ),
		upETA = upKAPPA,upALPHA = upKAPPA;
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ){
		upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	}
	upETA = upPARS.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upPARS.subvec(GI.at(5,0),GI.at(5,1));
	
	for(ct = 0; ct < QQ; ct++){
		if( upKAPPA.at(ct) == 0.0 ){
			if( ct == 0 ){
				PARS.at(0) = small_log_mu;
			} else {
				PARS.at(GI.at(2,0) + ct - 1) = neg_inf;
			}
		}
		
		if( upETA.at(ct) == 0.0 ){
			PARS.at(GI.at(3,0) + ct) = 0.0;
		}
		if( upALPHA.at(ct) == 0.0 ){
			PARS.at(GI.at(5,0) + ct) = 0.0;
		}
	}
	
	PARS.at(GI.at(1,0)) = 0.0;
	PARS.at(GI.at(4,0)) = 0.0;
	if( upPARS.at(GI.at(1,0)) == 0.0 ) PARS.at(GI.at(1,0)) = neg_inf;
	if( upPARS.at(GI.at(4,0)) == 0.0 ) PARS.at(GI.at(4,0)) = neg_inf;
	
}
