#include <RcppArmadillo.h>
#include "../inst/include/minor.h"
#include "../inst/include/negbinom.h"
#include "../inst/include/poisson.h"
#include "../inst/include/LM_reg.h"
#include "../inst/include/binom_beta.h"
#include "../inst/include/pack_LL.h"
#include "../inst/include/pack_GRAD.h"
#include "../inst/include/pack_HESS.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]



// --------------------
// CSeQTL Fundamental Functions: BFGS, LRT, fullOptTest

void Rcpp_CSeQTL_upPARS(const arma::uword& QQ,const arma::umat& GI,
	const arma::vec& upPARS){
	
	arma::vec upKAPPA = arma::zeros<arma::vec>(QQ);
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	
	Rcpp::Rcout << "   upPHI = " << upPARS.at(GI.at(1,0)) << "\n";
	Rcpp::Rcout << "   upKAPPA = " << upKAPPA.t();
	Rcpp::Rcout << "   upETA = " << upPARS.subvec(GI.at(3,0),GI.at(3,1)).t();
	Rcpp::Rcout << "   upPSI = " << upPARS.at(GI.at(4,0)) << "\n";
	Rcpp::Rcout << "   upALPHA = " << upPARS.subvec(GI.at(5,0),GI.at(5,1)).t();
	
}

void Rcpp_CSeQTL_BFGS(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,arma::vec& PARS,
	arma::uword& converge,const arma::mat& I_np,const arma::vec& upPARS,
	const arma::uword& max_iter,const double& eps,const double& gr_eps,
	const double& conv_eps,const bool& show){
	/*
	max_iter = 4e3; eps = 1e-10; gr_eps = 1e-2; conv_eps = 5e-5
	*/
	
	converge = 0;
	arma::uword iter = 0, jj, uu, reset_Bk = 0,
		np = GI.at(5,1)+1, QQ = GI.at(5,1) - GI.at(5,0) + 1;
	
	// Initialize parameters
	Rcpp_CSeQTL_control_PAR(PARS,QQ,GI,upPARS);
	// if( show ) Rcpp::Rcout << "\t##\n\tiPARS = " << PARS.t();
	if( show ) Rcpp::Rcout << "\t### Begin BFGS ...\n";
	
	arma::mat inv_Bk = I_np, ISYT = inv_Bk;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk, gr_k = curr_xk,
		p_k = curr_xk, s_k = curr_xk, y_k = curr_xk;
	double old_LL, new_LL, inv_norm_p_k, tmp_step, ys,
		fnscale = -1.0, curr_LL = 0.0, hess_shift = 1e-5;
	double orig_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
		LBC,PHASE,SNP,RHO,XX,GI,xk);
	
	while(iter < max_iter){
		// Calculate Direction p_k
		gr_k = fnscale * Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,
			PHASE,SNP,RHO,XX,GI,upPARS,xk);
		
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));

		// Line search for new xk
		uu = 0;
		old_LL = fnscale * Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
			LBC,PHASE,SNP,RHO,XX,GI,xk);
		
		// Less than 0 b/c fnscale = -1
		if( old_LL < 0.0 ){
			break; 
		}
		
		for(jj = 0; jj <= 20; jj++){
			tmp_step 	= inv_norm_p_k / std::pow(4.0,jj);
			new_xk 		= xk + tmp_step * p_k;
			new_LL = fnscale * Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
				LBC,PHASE,SNP,RHO,XX,GI,new_xk);
			if( new_LL < old_LL ){ // minimizing
				s_k = tmp_step * p_k;
				y_k = fnscale * Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,
					PHASE,SNP,RHO,XX,GI,upPARS,new_xk) - gr_k;
				if( y_k.has_nan() ) continue;
				ys = arma::dot(y_k,s_k);
				if( ys > 0.0 ){
					ISYT = I_np - (s_k * y_k.t()) / ys;
					inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
				}
				
				xk = new_xk;
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if( show ) printR_obj("\tReset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( show ) printR_obj("\tNo more update");
				break;
			}
		}
		
		if( reset_Bk > 5 ) break;
		
		// Check Convergence
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				Rcpp_norm(curr_xk - xk) < eps ){
				gr_k = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,
					PHASE,SNP,RHO,XX,GI,upPARS,xk);
				if( Rcpp_norm(gr_k) < eps ){
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	
	}
	
	// Update parameters
	PARS = xk;
	
	// Calculate LL, GRAD
	old_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,GI,PARS);
	gr_k = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,GI,upPARS,PARS);
	
	// Assessing convergence
	arma::mat hess = arma::zeros<arma::mat>(np,np),covar = hess;
	bool rcond_nz = false;
	double nG = Rcpp_norm(gr_k),nIHG = 0;
	if( show ) Rcpp::Rcout << "\tgr_eps = " << gr_eps << ";norm_GRAD = " << nG << "\n";
	if( nG < gr_eps ){
		hess = Rcpp_CSeQTL_calc_HESS(TREC,hap2,ASREC,PHASE,SNP,
			RHO,XX,GI,PARS,upPARS,I_np,hess_shift);
		Rcpp_CSeQTL_hessBR(hess,GI,upPARS,rcond_nz,false);
		arma::uvec nz = arma::find( hess.diag() != 0.0 );
		if( rcond_nz ){
			covar.submat(nz,nz) = arma::inv(-1.0 * hess.submat(nz,nz));
			nIHG = Rcpp_norm(covar * gr_k);
			if( arma::any( arma::diagvec(covar.submat(nz,nz)) < 0.0 ) ){
				converge = 0;
				if( show ){
					Rcpp::Rcout << "\t***Converge fail: negative variance!\n";
					arma::vec vars = arma::diagvec(covar);
					if( arma::any( vars.subvec(GI.at(0,0),GI.at(0,1)) < 0.0 ) ) 
						Rcpp::Rcout << "\tvar_beta = " << vars.subvec(GI.at(0,0),GI.at(0,1)).t();
					if( vars.at(GI.at(0,1)) < 0.0 ) 
						Rcpp::Rcout << "\tvar_logPhi = " << vars.at(GI.at(1,0)) << "\n";
					if( QQ > 1 && arma::any( vars.subvec(GI.at(2,0),GI.at(2,1)) < 0.0 ) ) 
						Rcpp::Rcout << "\tvar_logK = " << vars.subvec(GI.at(2,0),GI.at(2,1)).t();
					if( arma::any( vars.subvec(GI.at(3,0),GI.at(3,1)) < 0.0 ) ) 
						Rcpp::Rcout << "\tvar_logE = " << vars.subvec(GI.at(3,0),GI.at(3,1)).t();
					if( vars.at(GI.at(4,0)) < 0.0 ) 
						Rcpp::Rcout << "\tvar_logPsi = " << vars.at(GI.at(4,0)) << "\n";
					if( arma::any( vars.subvec(GI.at(5,0),GI.at(5,1)) < 0.0 ) ) 
						Rcpp::Rcout << "\tvar_logA = " << vars.subvec(GI.at(5,0),GI.at(5,1)).t();
				}
			} else if( nIHG < conv_eps ){
				converge = 1;
			} else {
				if( show ) Rcpp::Rcout << "\t***Converge fail: norm_iHG > threshold!\n";
			}
		} else {
			if( show ) Rcpp::Rcout << "\t***Converge fail: rcond issue!\n";
			converge = 0;
		}
	} else {
		if( show ) Rcpp::Rcout << "\t***Converge fail: norm_GRAD > threshold!\n";
	}
	
	if( show ){
		Rcpp::Rcout << "\tLL_old = " << orig_LL << "; Iter = " << iter 
			<< "; LL_new = " << old_LL << "\n";
		// Rcpp::Rcout << "\tGR = " << gr_k.t();
		nG = Rcpp_norm(gr_k);
		Rcpp::Rcout << "\tConvergence Indicators: \n"
			<< "\t   NormGrad = " << nG << "; NormIHessGrad = " << nIHG << "\n";
		arma::vec upKAPPA = arma::zeros<arma::vec>(QQ), upETA = upKAPPA, upALPHA = upKAPPA;
			upKAPPA.at(0) = upPARS.at(0);
			if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
			upETA = upPARS.subvec(GI.at(3,0),GI.at(3,1));
			upALPHA = upPARS.subvec(GI.at(5,0),GI.at(5,1));
		Rcpp::Rcout << "\tupKAP = " << upKAPPA.t();
		Rcpp::Rcout << "\tupETA = " << upETA.t();
		Rcpp::Rcout << "\tupALP = " << upALPHA.t();
		Rcpp::Rcout << "\t   BETA = " << PARS.subvec(GI.at(0,0),GI.at(0,1)).t();
		Rcpp::Rcout << "\t   PHI = " << std::exp(PARS.at(GI.at(1,0))) << "\n";
		arma::vec kappa = arma::ones<arma::vec>(QQ);
		if(QQ > 1) kappa.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
		Rcpp::Rcout << "\t   KAP = " << kappa.t();
		Rcpp::Rcout << "\t   ETA = " << arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))).t();
		Rcpp::Rcout << "\t   PSI = " << std::exp(PARS.at(GI.at(4,0))) << "\n";
		Rcpp::Rcout << "\t   ALP = " << arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1)).t());
		Rcpp::Rcout << "\tuseASREC = " << arma::any(PHASE == 1) << "\n";
		Rcpp::Rcout << "\t   Convergence Status = ";
		if( converge == 1 ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
		arma::mat MU = Rcpp_CSeQTL_MU(GI,PARS);
		Rcpp::Rcout << "\tMUa = " << MU.row(0);
		Rcpp::Rcout << "\tMUb = " << MU.row(1);
		
		// Calculate variances
		hess = Rcpp_CSeQTL_calc_HESS(TREC,hap2,ASREC,PHASE,SNP,
			RHO,XX,GI,PARS,upPARS,I_np,hess_shift);
		arma::uvec nz = arma::find( hess.diag() != 0.0 );
		Rcpp_CSeQTL_hessBR(hess,GI,upPARS,rcond_nz,false);
		if( rcond_nz ){
			covar.submat(nz,nz) = arma::inv(-1.0 * hess.submat(nz,nz));
			// Rcpp::Rcout << "\tVAR_BETA = " << arma::diagvec(covar.submat(GI.at(0,0),GI.at(0,0),GI.at(0,1),GI.at(0,1))).t();
			Rcpp::Rcout << "\tVAR_logPHI = " << covar.at(GI.at(1,0),GI.at(1,0)) << ";\n";
			if( covar.at(GI.at(1,0),GI.at(1,0)) < 0.0 ){
				Rcpp::Rcout << "\t****Negative VAR_logPHI!\n";
			}
			arma::mat eqtl_vars = arma::zeros<arma::mat>(3,QQ);
			// Kappa, Eta, Alpha variances
			eqtl_vars.at(0,0) = covar.at(0,0);
			if(QQ > 1) eqtl_vars(0,arma::span(1,QQ-1)) = arma::diagvec(covar.submat(GI.at(2,0),
				GI.at(2,0),GI.at(2,1),GI.at(2,1))).t();
			eqtl_vars.row(1) = arma::diagvec(covar.submat(GI.at(3,0),GI.at(3,0),GI.at(3,1),GI.at(3,1))).t();
			eqtl_vars.row(2) = arma::diagvec(covar.submat(GI.at(5,0),GI.at(5,0),GI.at(5,1),GI.at(5,1))).t();
			eqtl_vars.print("\tVAR_logEQTLS = ");
			Rcpp::Rcout << "\tVAR_logPSI = " << covar.at(GI.at(4,0),GI.at(4,0));
			if( covar.at(GI.at(4,0),GI.at(4,0)) < 0.0 ){
				Rcpp::Rcout << "; ****Negative VAR_logPSI!\n";
			}
		}
		
		Rcpp::Rcout << "\n";
	}
	
}

void Rcpp_CSeQTL_optPARAMS(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,arma::uword& converge,const arma::uword& QQ,
	const arma::umat& GI,const arma::mat& I_np,const arma::vec& iBETA,const double& iPHI,
	const arma::vec& upPARS,arma::vec& PARS,const bool& useASREC,const arma::uword& max_iter,
	const double& eps,const double& gr_eps,const double& conv_eps,const bool& show){
	
	// Want to adaptively optimize, to avoid re-optimizing 
	// If QQ > 1 and KAPPA2Q parameters all == 1, start with KAPPA. Otherwise,
	// If all ETA == 1, start with ETA optimization w/o ASREC. Otherwise,
	// If Psi == 1, start at ETA optimization w/ ASREC
	
	converge = 0;
	bool any_PHASE = arma::any(PHASE == 1);
	arma::vec upPARS_2 = upPARS,
		upKAPPA = arma::zeros<arma::vec>(QQ),
		upETA = upKAPPA, upALPHA = upETA;
	
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upETA = upPARS.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upPARS.subvec(GI.at(5,0),GI.at(5,1));
	
	if( QQ > 1 ){
		if( arma::any(PARS.subvec(GI.at(2,0),GI.at(2,1)) == 0.0 && upKAPPA.subvec(1,QQ-1) == 1.0) ){
			if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize BETA,PHI,KAPPA ...\n";
			PARS.zeros();
			Rcpp_NB_iPARS(iBETA,iPHI,GI,PARS);
			
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)).zeros();
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
				SNP,RHO,XX,GI,PARS,converge,I_np,upPARS_2,
				max_iter,eps,gr_eps,conv_eps,show);
		}
	} else {
		if( arma::any(PARS.subvec(GI.at(0,0),GI.at(1,1)) == 0.0) ){
			PARS.zeros();
			Rcpp_NB_iPARS(iBETA,iPHI,GI,PARS);
			if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize BETA,PHI ...\n";
			
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)).zeros();
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
				SNP,RHO,XX,GI,PARS,converge,I_np,upPARS_2,
				max_iter,eps,gr_eps,conv_eps,show);
		}
	}

	if( arma::any(upETA == 1.0) ){
		if( arma::any(PARS.subvec(GI.at(3,0),GI.at(3,1)) == 0.0 && upETA == 1.0) ){
			if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize BETA,PHI,KAPPA,ETA w/o ASREC\n";
			
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
				SNP,RHO,XX,GI,PARS,converge,I_np,upPARS_2,
				max_iter,eps,gr_eps,conv_eps,show);
		}
		
		if( useASREC && any_PHASE ){
			if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize BETA,PHI,KAPPA,ETA w/ ASREC\n";
			
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
				SNP,RHO,XX,GI,PARS,converge,I_np,upPARS_2,
				max_iter,eps,gr_eps,conv_eps,show);
		}
	}
	
	if( useASREC && any_PHASE && arma::any(PARS.subvec(GI.at(5,0),GI.at(5,1)) == 0.0 && upALPHA == 1.0) ){
		if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize BETA,PHI,KAPPA,ETA,PSI,ALPHA\n";
		
		upPARS_2 = upPARS;
		
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
			SNP,RHO,XX,GI,PARS,converge,I_np,upPARS_2,
			max_iter,eps,gr_eps,conv_eps,show);
	}
	
	if( show ) Rcpp::Rcout << "\n\t###------- Joint Optimize FINAL\n";
	if( useASREC && any_PHASE ){
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
			SNP,RHO,XX,GI,PARS,converge,I_np,upPARS,
			max_iter,eps,gr_eps,conv_eps,false);
	} else {
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
			SNP,RHO,XX,GI,PARS,converge,I_np,upPARS,
			max_iter,eps,gr_eps,conv_eps,false);
	}
	
	if( show ){
		Rcpp::Rcout << "\tFINAL convergence status = ";
		if( converge == 1 ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
	}
	
}

void Rcpp_CSeQTL_gs_signif(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::vec& PARS_BPK,const arma::umat& GI,
	const arma::mat& I_np,const arma::vec& upPARS_full,const arma::uword& QQ,
	const arma::uword& np,arma::mat& LRT,arma::vec& ETA,const bool& TRECASE,
	const arma::uword& max_iter,const double& eps,const double& gr_eps,
	const double& conv_eps,const bool& show){
	
	// Note: TRECASE = true means we test each ct assuming its cis-eQTL, we'll combine TReC and ASE models
	// LRT = (LRT statistic,DF,pvalue)
	LRT.zeros();
	arma::uword converge = 0, ct;
	arma::vec full_PARS = PARS_BPK, tmp_PARS = PARS_BPK,
		tmp_ETA = arma::ones<arma::vec>(QQ), upPARS = upPARS_full;
	arma::uvec PHASE_2 = PHASE;
	if( !TRECASE ) PHASE_2.zeros();
	
	if( show ){
		if( TRECASE ){
			Rcpp::Rcout << "\n############\n#### TReCASE Section\n############\n";
		} else {
			Rcpp::Rcout << "\n############\n#### TReC Section\n############\n";
		}
		Rcpp_CSeQTL_upPARS(QQ,GI,upPARS_full);
	}

	// Get full TReC model MLEs
	if( show ){
		if( arma::any(PHASE_2 != 0.0) ){
			Rcpp::Rcout << "\n####### Optimize BETA,PHI,KAPPA,ETA,PSI!\n";
		} else {
			Rcpp::Rcout << "\n####### Optimize BETA,PHI,KAPPA,ETA!\n";
		}
	}
	
	upPARS = upPARS_full;
	upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	
	Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE_2,
		SNP,RHO,XX,GI,full_PARS,converge,I_np,upPARS,
		max_iter,eps,gr_eps,conv_eps,show);
	double full_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
		LBC,PHASE_2,SNP,RHO,XX,GI,full_PARS), larger_LL, reduced_LL;
	
	if( !TRECASE ){
		ETA = arma::exp(full_PARS.subvec(GI.at(3,0),GI.at(3,1)));
	} else {
		ETA.ones();
	}
	
	// Get null model MLEs
	for(ct = 0; ct < QQ; ct++){
		
		if( TRECASE ){ // Assuming ct's eQTL is cis-eQTL, the rest are trans-eQTL
			
			// -----
			// Larger model estimation process: estimate all parameters 
			// -----
			//		1) Start with beta,phi,kappa,eta,psi estimates (full_PARS)
			//		2) Update beta,phi,kappa,eta,psi,alpha
			if( show ) Rcpp::Rcout << "#---------- Larger model ct = " << ct + 1 << "\n";
			
			// upETA = upETA_full; upALPHA = upALPHA_full; upALPHA.at(ct) = 0.0;
			upPARS = upPARS_full;
			upPARS.at(GI.at(5,0) + ct) = 0.0;
			
			// start PARS with BETA,PHI,KAPPA,ETA,PSI estimate
			tmp_PARS = full_PARS;
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE_2,
				SNP,RHO,XX,GI,tmp_PARS,converge,I_np,upPARS,
				max_iter,eps,gr_eps,conv_eps,show);
			larger_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
				LBC,PHASE_2,SNP,RHO,XX,GI,tmp_PARS);
			tmp_ETA = arma::exp(tmp_PARS.subvec(GI.at(3,0),GI.at(3,1)));
			ETA.at(ct) = tmp_ETA.at(ct);
			
			// -----
			// Reduced model estimation process:
			// -----
			// upETA = upETA_full; upETA.at(ct) = 0;
			upPARS = upPARS_full;
			upPARS.at(GI.at(3,0) + ct) = 0;
			
			//		1) Start with beta,phi,kappa (init_PARS)
			tmp_PARS = PARS_BPK;
			
			//		2) Update beta,phi,kappa,eta,psi(constrained) and set all alpha = 0
			if( show ) Rcpp::Rcout << "#---------- Reduced model ct = " << ct + 1 << "\n";
			upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE_2,
				SNP,RHO,XX,GI,tmp_PARS,converge,I_np,upPARS,
				max_iter,eps,gr_eps,conv_eps,show);
			
			//		3) Update beta,phi,kappa,eta,psi(constrained) and any alphas
			// upALPHA = upALPHA_full % upETA;
			upPARS.subvec(GI.at(5,0),GI.at(5,1)) = 
				upPARS_full.subvec(GI.at(5,0),GI.at(5,1)) % upPARS.subvec(GI.at(3,0),GI.at(3,1));
			
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE_2,
				SNP,RHO,XX,GI,tmp_PARS,converge,I_np,upPARS,
				max_iter,eps,gr_eps,conv_eps,show);
			
			reduced_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
				LBC,PHASE_2,SNP,RHO,XX,GI,tmp_PARS);
			
		} else {
			
			// Larger model
			larger_LL = full_LL;
			
			if( show ) Rcpp::Rcout << "#---------- Reduced model ct = " << ct + 1 << "\n";
			
			// Reduced model
			// upETA = upETA_full; upETA.at(ct) = 0; 
			upPARS = upPARS_full;
			upPARS.at(GI.at(3,0) + ct) = 0.0;
			upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			
			tmp_PARS = PARS_BPK;
			if( arma::all(upPARS.subvec(GI.at(3,0),GI.at(3,1)) == upPARS_full.subvec(GI.at(3,0),GI.at(3,1))) ){
				reduced_LL = larger_LL;
			} else {
				Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE_2,
					SNP,RHO,XX,GI,tmp_PARS,converge,I_np,upPARS,
					max_iter,eps,gr_eps,conv_eps,show);
				reduced_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
					LBC,PHASE_2,SNP,RHO,XX,GI,tmp_PARS);
			}
			
		}
		
		if( show ){
			Rcpp::Rcout << "\tFull LL = " << full_LL
				<< "; Larger LL = " << larger_LL
				<< "; Reduced LL = " << reduced_LL << "\n";
		}
		
		LRT.at(ct,0) = 2.0 * (larger_LL - reduced_LL);
		if( LRT.at(ct,0) < 0 ) LRT.at(ct,0) = 0.0;
		LRT.at(ct,1) = 1;
		LRT.at(ct,2) = 1.0 - R::pchisq(LRT.at(ct,0),LRT.at(ct,1),1,0);
		
	}
	
	if( show ){
		if( TRECASE ){
			LRT.print("\n####### TReCASE eQTL Significance LRT results:");
		} else {
			LRT.print("\n####### TReC eQTL Significance LRT results:");
		}
	}
	
}

arma::vec Rcpp_CSeQTL_fullOptTest(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,const arma::mat& I_np,
	const arma::uword& QQ,const arma::uword& np,const arma::vec& iBETA,
	const double& iPHI,const arma::vec& upPARS_0,arma::vec& PARS,
	bool& useASREC,bool& useTRIM,const arma::uword& max_iter,const double& eps,
	const double& gr_eps,const double& conv_eps,const bool& show){
	
	// Run full model optimization, determine values of upKAPPA/upETA/upALPHA based on Hessian
	arma::uword converge = 0;
	double mu_thres = 2.0, hess_shift = 1e-5;
	arma::vec upKAPPA = arma::zeros<arma::vec>(QQ),
		upETA = upKAPPA,upALPHA = upETA,
		PARS_BPK = arma::zeros<arma::vec>(np),
		upPARS = upPARS_0,upPARS_2 = upPARS;
	
	// Set upKAPPA
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ){
		upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	}
	// upETA = upKAPPA % upETA_0
	upPARS.subvec(GI.at(3,0),GI.at(3,1)) %= upKAPPA;
	upETA = upPARS.subvec(GI.at(3,0),GI.at(3,1));
	// upALPHA = upALPHA_0 % upETA
	upPARS.subvec(GI.at(5,0),GI.at(5,1)) %= upPARS.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upPARS.subvec(GI.at(5,0),GI.at(5,1));
	
	// Full model optimization
	upPARS_2 = upPARS;
	upPARS_2.subvec(GI.at(3,0),GI.at(3,1)).zeros();
	upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	Rcpp_CSeQTL_optPARAMS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,SNP*0,RHO,
		XX,converge,QQ,GI,I_np,iBETA,iPHI,upPARS_2,PARS,
		false,max_iter,eps,gr_eps,conv_eps,show);
	PARS_BPK = PARS;
	Rcpp_CSeQTL_optPARAMS(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,
		XX,converge,QQ,GI,I_np,iBETA,iPHI,upPARS,PARS,
		useASREC,max_iter,eps,gr_eps,conv_eps,show);
	
	// Calculate and check mu_total = mu_Aq + mu_Bq > 0
	if( show ) Rcpp::Rcout << "fullOptTest(): check mu_total ...\n";
	arma::mat MU = Rcpp_CSeQTL_MU(GI,PARS);
	arma::vec mu_total = arma::sum(MU,0).t(), min_mu = arma::min(MU,0).t();
	bool check_mu = arma::any(mu_total > 0.0 && mu_total < mu_thres);
	bool small_mu = arma::any(min_mu < 1.0 && upETA == 1.0);
	if( check_mu ){
		if( show ) Rcpp::Rcout << "mu_total = " << mu_total.t();
		arma::uvec mu_idx = arma::find(mu_total > 0.0 && mu_total < mu_thres);
		double min_mu_total = arma::min(mu_total(mu_idx));
		
		if( arma::any(mu_total == min_mu_total && upALPHA == 1.0) ){
			if( show ) Rcpp::Rcout << "***-------(1) Small MU, constrain alpha\n";
			arma::vec upALPHA_2 = upALPHA;
			arma::uvec check_idx = arma::find(mu_total == min_mu_total && upALPHA_2 == 1.0);
			upALPHA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)) = upALPHA_2;
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		} else if( arma::any(mu_total == min_mu_total && upETA == 1.0) ){
			arma::vec upETA_2 = upETA;
			arma::uvec check_idx = arma::find(mu_total == min_mu_total && upETA_2 == 1.0);
			upETA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
			if( show ) Rcpp::Rcout << "***------- Small MU, constrain eta[" << check_idx.at(0) + 1 << "]\n";
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		} else if( arma::any(mu_total == min_mu_total && upKAPPA == 1.0) ){
			arma::uvec check_idx = arma::find(mu_total == min_mu_total && upKAPPA == 1.0);
			if( arma::any(check_idx == 0) ){
				// arma::mat upPARAM = arma::zeros<arma::mat>(QQ,3);
				// upPARAM.col(0).ones();
				// upPARAM.at(0,0) = 0;
				upPARS_2 = upPARS;
				upPARS_2.at(0) = 0;
				if( show ) Rcpp::Rcout << "***------- Small ref MU, constrain kappa[" << check_idx.at(0) + 1 << "]\n";
				// return upPARAM;
				return upPARS_2;
			} else {
				arma::vec upKAPPA_2 = upKAPPA;
				upKAPPA_2(check_idx).zeros();
				if( show ) Rcpp::Rcout << "***------- Small nonref MU, constrain kappa[" << check_idx.at(0) + 1 << "]\n";
				PARS.zeros(); // reset to 0 to re-optimize
				upPARS_2 = upPARS;
				upPARS_2.subvec(GI.at(2,0),GI.at(2,1)) = upKAPPA_2.subvec(1,QQ-1);
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			}
		}
	}
	
	// Check 0 < min(mu_Aq,mu_Ba) < 1
	if( show ) Rcpp::Rcout << "fullOptTest(): check min(mu_Aq,mu_Bq) ...\n";
	if( small_mu ){
		arma::vec upETA_2 = upETA;
		arma::uvec check_idx = arma::find( min_mu < 1.0 && upETA_2 == 1.0 );
		arma::uvec check_idx2 = arma::find( arma::min(min_mu(check_idx)) == min_mu );
		upETA_2(check_idx2).zeros();
		upPARS_2 = upPARS;
		upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
		PARS = PARS_BPK;
		if( show ) Rcpp::Rcout << "***------- Small allelic MU: Constrain eta[" << check_idx2.at(0) + 1 << "]\n";
		return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
			GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
			gr_eps,conv_eps,show);
	}
	
	// If converge = 1, we're done. Otherwise need to check for negative or very large variances
	if( converge == 1 ){
		return upPARS;
	}
	
	// Calculate/Approximate Hessian and possible the covariance
	if( show ) Rcpp::Rcout << "fullOptTest(): check LL/GRAD/hess/covar ...\n";
	double LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,GI,PARS);
	arma::vec GRAD = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,GI,upPARS,PARS);
	arma::mat hess = arma::zeros<arma::mat>(np,np), covar = hess;
	hess = Rcpp_CSeQTL_calc_HESS(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,GI,PARS,upPARS,I_np,hess_shift);
	bool rcond_nz = false;
	Rcpp_CSeQTL_hessBR(hess,GI,upPARS,rcond_nz,false);
	arma::uvec nz = arma::find( hess.diag() != 0.0 );
	arma::mat var_eqtl = arma::zeros<arma::mat>(3,QQ);
	if( rcond_nz ){ 
		covar.submat(nz,nz) = arma::inv(-1.0 * hess.submat(nz,nz));
		var_eqtl.at(0,0) = covar.at(0,0);
		if( QQ > 1 ) var_eqtl(0,arma::span(1,QQ-1)) = covar.submat(GI.at(2,0),GI.at(2,0),GI.at(2,1),GI.at(2,1)).diag().t();
		var_eqtl.row(1) = covar.submat(GI.at(3,0),GI.at(3,0),GI.at(3,1),GI.at(3,1)).diag().t();
		var_eqtl.row(2) = covar.submat(GI.at(5,0),GI.at(5,0),GI.at(5,1),GI.at(5,1)).diag().t();
	}
	
	// If we can't invert the hessian and upALPHA == 1 ...
	if( !rcond_nz ){
		if( useASREC && arma::any(upALPHA == 1.0) && arma::any(PHASE == 1.0) ){
			// Set upALPHA = 0
			upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
			return upPARS;
		} else {
			Rcpp::stop("Need to code for this");
		}
	}
	
	if( show ){
		Rcpp::Rcout << "\tLL = " << LL 
			<< "; nGRAD = " << Rcpp_norm(GRAD) << "\n";
		Rcpp::Rcout << "\trcond_nz = " << rcond_nz << "\n";
		
		if( rcond_nz ) var_eqtl.print("var_eqtl = ");
	}
	
	// Determine if Var(log_PHI) or Var(log_PSI) are negative or PHI/PSI are too small
	if( rcond_nz ){
		if( upPARS.at(GI.at(1,0)) == 1.0 ){
			double logPHI_thres = -9.0;
			if( covar.at(GI.at(1,0),GI.at(1,0)) < 0.0 || PARS.at(GI.at(1,0)) < logPHI_thres ){
				upPARS_2 = upPARS;
				upPARS_2.at(GI.at(1,0)) = 0.0;
				PARS.zeros();
				PARS.subvec(GI.at(0,0),GI.at(0,1)) = iBETA;
				PARS.at(GI.at(1,0)) = -arma::datum::inf;
				if( show ) Rcpp::Rcout << "***------- Small PHI: Fit TReC with Poisson!\n";
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			}
		}
		
		if( upPARS.at(GI.at(4,0)) == 1.0 ){
			double logPSI_thres = -9.0;
			if( covar.at(GI.at(4,0),GI.at(4,0)) < 0.0 || PARS.at(GI.at(4,0)) < logPSI_thres ){
				upPARS_2 = upPARS;
				upPARS_2.at(GI.at(4.0)) = 0.0;
				PARS.zeros();
				PARS.subvec(GI.at(0,0),GI.at(0,1)) = iBETA;
				PARS.at(GI.at(4,0)) = -arma::datum::inf;
				if( show ) Rcpp::Rcout << "***------- Small PSI: Fit ASReC with Binomial!\n";
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			}
		}
	}
	
	// If TReC was trimmed and only KAPPA being assessed,
	//	end this optimization using untrimmed TReC
	
	// Check magnitude of variances
	if( rcond_nz && show ){
		var_eqtl.print("var_eqtl = ");
		Rcpp::Rcout << "fullOptTest(): var[log(psi)] = "
			<< covar.at(GI.at(4,0),GI.at(4,0)) << "\n";
	}
	
	// Check ALPHA variances
	if( rcond_nz && arma::any(upALPHA == 1.0) ){
		arma::vec var_alpha = var_eqtl.row(2).t(),
			var_eta = var_eqtl.row(1).t();
		bool neg_alpha = arma::any(var_alpha < 0.0),
			neg_eta = arma::any(var_eta < 0.0);
		arma::vec upALPHA_2 = upALPHA;
		
		if( neg_alpha ){ // Check for negative alpha variances!
			if( neg_eta ){
				arma::uvec check_idx = arma::find(var_eta < 0.0);
				if( arma::any(check_idx == 0) ){
					arma::vec upETA_2 = upKAPPA;
					// arma::vec upETA_2 = arma::zeros<arma::vec>(QQ);
					// upETA_2.at(0) = upPARS.at(0);
					// if( QQ > 1 ) upETA_2.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
					upETA_2.at(0) = 0.0;
					upPARS_2 = upPARS;
					upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
					PARS = PARS_BPK;
					if( show ) Rcpp::Rcout << "***------- Negative Variance: constrain eta[1] and start over\n";
					return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
						GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
						gr_eps,conv_eps,show);
				}
			}
			arma::uvec check_idx = arma::find(arma::min(var_alpha) == var_alpha);
			upALPHA_2(check_idx).zeros();
			if( show ) Rcpp::Rcout << "***------- Negative Variance: constrain alpha[" << check_idx.at(0) + 1 << "]\n";
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)) = upALPHA_2;
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
		
		if( neg_eta ){ // Check for negative eta variances!
			arma::uvec check_idx = arma::find(var_eta < 0.0);
			arma::vec upETA_2 = upETA;
			if( check_idx.n_elem == 1 ){
				upETA_2(check_idx).zeros();
				upPARS_2 = upPARS;
				upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
				PARS = PARS_BPK;
				if( show ) Rcpp::Rcout << "***------- Negative Variance: constrain eta[" << check_idx.at(0) + 1 << "]\n";
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			} else if( arma::any(check_idx == 0) ){
				upETA_2.at(0) = 0.0;
				upPARS_2 = upPARS;
				upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
				PARS = PARS_BPK;
				if( show ) Rcpp::Rcout << "***------- Negative Variance: constrain eta[" << 1 << "]\n";
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			} else {
				arma::uvec check_idx2 = arma::find(arma::min(var_eta) == var_eta);
				upETA_2(check_idx2).zeros();
				upPARS_2 = upPARS;
				upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
				PARS = PARS_BPK;
				if( show ) Rcpp::Rcout << "***------- Negative Variance: constrain eta[" << check_idx2.at(0) + 1 << "]\n";
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			}
		}
		
		// Look for largest variance among ALPHA and ETA
		double max_var_eta = arma::max(var_eta),
			max_var_alpha = arma::max(var_alpha),
			max_var = max_var_eta;
		if( max_var_alpha > max_var_eta ) max_var = max_var_alpha;
		if( max_var == max_var_alpha ){
			arma::uvec check_idx = arma::find(var_alpha == max_var);
			upALPHA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)) = upALPHA_2;
			if( show ) Rcpp::Rcout << "***------- Largest Variance: constrain alpha[" << check_idx.at(0) + 1 << "]\n";
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		} else {
			arma::vec upETA_2 = upETA;
			arma::uvec check_idx = arma::find(var_eta == max_var);
			upETA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
			if( show ) Rcpp::Rcout << "***------- Largest Variance: constrain eta[" << check_idx.at(0) + 1 << "]\n";
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
		
	}
	
	// Check PSI variance
	if( false && covar.at(GI.at(4,0),GI.at(4,0)) < 0.0 ){ 
		if( show ) Rcpp::Rcout << "fullOptTest(): check log_PSI variance ...\n";
		useASREC = false;
		if( show ) Rcpp::Rcout << "***------- Don't use ASREC, potential disagreement with TReC\n";
		PARS = PARS_BPK;
		arma::vec upETA_2 = upKAPPA;
		// return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
			// GI,I_np,QQ,np,iBETA,iPHI,upKAPPA,upETA_2,upALPHA*0,PARS,useASREC,useTRIM,max_iter,eps,
			// gr_eps,conv_eps,show);
	}
	
	// Check ETA variances: (1) negative eta variance, (2) small allelic MU, (3) negative psi variance
	if( rcond_nz && arma::any(upETA == 1.0) ){
		if( show ) Rcpp::Rcout << "fullOptTest(): check log_ETA variance ...\n";
		arma::vec var_eta = var_eqtl.row(1).t(),
			var_kap = var_eqtl.row(0).t();
		bool neg_eta = arma::any(var_eta < 0.0);
		arma::vec upETA_2 = upETA;
		
		if( neg_eta ){ // Check for negative variances first
			arma::uvec check_idx = arma::find(var_eta < 0.0); 
			if( check_idx.n_elem == 1 ){
				upETA_2(check_idx).zeros();
				if( show ) Rcpp::Rcout << "***------- Negative Variance: Constrain eta[" << check_idx.at(0) + 1 << "]\n";
			} else if( arma::any(check_idx == 0) ){ // if first cell type has negative ETA variance, constrain
				upETA_2.at(0) = 0.0;
				if( show ) Rcpp::Rcout << "***------- Negative Variance: Constrain eta[" << 1 << "]\n";
			} else if( arma::any(var_kap < 0.0 && var_eta < 0.0) ){ // multiple negative variances
				arma::uvec check_idx2 = arma::find(var_kap < 0.0 && var_eta < 0.0);
				arma::uvec check_idx3 = arma::find(arma::min(min_mu(check_idx2)) == min_mu);
				upETA_2(check_idx3).zeros();
				if( show ) Rcpp::Rcout << "***------- Negative Variance: Constrain eta[" << check_idx3.at(0) + 1 << "]\n";
			} else { // multiple neg variance, remove smallest
				arma::uvec check_idx2 = arma::find(arma::min(var_eta) == var_eta);
				upETA_2(check_idx2).zeros();
				if( show ) Rcpp::Rcout << "***------- Negative Variance: Constrain eta[" << check_idx2.at(0) + 1 << "]\n";
			}
			upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
			PARS = PARS_BPK;
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
		
		if( false ){ // Disable this code, seems unnecessary
			// Look for largest variance among ETA and KAPPA
			double max_var_eta = arma::max(var_eta),
				max_var_kap = arma::max(var_kap),
				max_var = max_var_eta;
			if( max_var_kap > max_var_eta ) max_var = max_var_kap;
			if( max_var == max_var_kap ){
				arma::uvec check_idx = arma::find(var_kap == max_var);
				arma::vec upKAPPA_2 = upKAPPA;
				upKAPPA_2(check_idx).zeros();
				if( show ) Rcpp::Rcout << "***------- Largest Variance: constrain kappa[" << check_idx.at(0) + 1 << "]\n";
				upPARS_2 = upPARS;
				upPARS_2.at(0) = upKAPPA_2.at(0);
				if( QQ > 1 ) upPARS_2.subvec(GI.at(2,0),GI.at(2,1)) = upKAPPA_2.subvec(1,QQ-1);
				PARS.zeros();
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			} else {
				arma::vec upETA_2 = upETA;
				arma::uvec check_idx = arma::find(var_eta == max_var);
				upETA_2(check_idx).zeros();
				if( show ) Rcpp::Rcout << "***------- Largest Variance: constrain eta[" << check_idx.at(0) + 1 << "]\n";
				upPARS_2 = upPARS;
				upPARS_2.subvec(GI.at(3,0),GI.at(3,1)) = upETA_2;
				PARS = PARS_BPK;
				return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
					GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
					gr_eps,conv_eps,show);
			}
			
		}
		
	}
	
	if( false && useTRIM && rcond_nz && arma::all(upETA == 0.0) 
		&& covar.at(GI.at(1,0),GI.at(1,0)) < 0.0 ){
		// Want to use untrimmed TREC
		useTRIM = false;
		return arma::ones<arma::vec>(np);
	}
	
	// Check KAPPA final, biggest variance set cell type to 0
	if( QQ > 1 && rcond_nz && arma::any(mu_total > 0.0 && upKAPPA == 1.0) 
		&& arma::all(upETA == 0.0) ){
		if( show ) Rcpp::Rcout << "fullOptTest(): check log_KAP variance ...\n";
		
		arma::vec var_kap = var_eqtl.row(0).t();
		arma::vec upKAPPA_2 = upKAPPA;
		PARS.zeros();
		bool neg_var = arma::any(var_kap < 0.0);
		
		if( neg_var ){
			arma::uvec check_idx = arma::find(arma::min(var_kap) == var_kap);
			upKAPPA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.at(0) = upKAPPA_2.at(0);
			if( QQ > 1 ) upPARS_2.subvec(GI.at(2,0),GI.at(2,1)) = upKAPPA_2.subvec(1,QQ-1);
			if( show ) Rcpp::Rcout << "***------- Negative variance, constrain kappa[" << check_idx.at(0) + 1 << "]\n";
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
		
		// If mu_total[1] isn't the largest, we'll need to end this function and swap cell types before proceeding
		arma::uvec max_ct = arma::find(mu_total == arma::max(mu_total));
		if( max_ct.at(0) != 0 ){
			// Need to swap cell types!
			return upPARS;
		}
		
		if( false ){ // Disable this code, seems unnecessary
			arma::uvec check_idx = arma::find(arma::max(var_kap) == var_kap);
			upKAPPA_2(check_idx).zeros();
			upPARS_2 = upPARS;
			upPARS_2.at(0) = upKAPPA_2.at(0);
			if( QQ > 1 ) upPARS_2.subvec(GI.at(2,0),GI.at(2,1)) = upKAPPA_2.subvec(1,QQ-1);
			if( show ) Rcpp::Rcout << "***------- Big variance, constrain kappa[" << check_idx.at(0) + 1 << "]\n";
			return Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
				GI,I_np,QQ,np,iBETA,iPHI,upPARS_2,PARS,useASREC,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
		
	}
	
	return upPARS;
}

arma::vec Rcpp_CSeQTL_run_fullModel(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,const arma::mat& I_np,
	const arma::uword& QQ,const arma::uword& np,const arma::vec& iBETA,
	const double& iPHI,const arma::vec& upPARS_0,bool& useASREC,
	const arma::uword& max_iter,const double& eps,const double& gr_eps,
	const double& conv_eps,const bool& show){
	
	bool PHASE_available = arma::any(PHASE == 1);
	if( show ){
		Rcpp::Rcout << "\n###### Check Full Model\n######\n";
		Rcpp::Rcout << "PHASE_available = " << PHASE_available << "\n";
	}
	useASREC = false;
	arma::vec upPARS = upPARS_0, PARS = upPARS_0,
		upKAPPA = arma::zeros<arma::vec>(QQ),
		upETA = upKAPPA, upALPHA = upETA;
	
	if( show ){
		Rcpp::Rcout << "### Model Check input \n";
		Rcpp::Rcout << "\tupPHI = " << upPARS_0.at(GI.at(1,0)) << "\n";
		upKAPPA.at(0) = upPARS_0.at(0);
		if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS_0.subvec(GI.at(2,0),GI.at(2,1));
		Rcpp::Rcout << "\tupKAPPA = " << upKAPPA.t();
		Rcpp::Rcout << "\tupETA = " << upPARS_0.subvec(GI.at(3,0),GI.at(3,1)).t();
		Rcpp::Rcout << "\tupPSI = " << upPARS_0.at(GI.at(4,0)) << "\n";
		Rcpp::Rcout << "\tupALPHA = " << upPARS_0.subvec(GI.at(5,0),GI.at(5,1)).t();
	}
	
	// Check kappa
	if( show ) Rcpp::Rcout << "### Check KAPPA ...\n";
	PARS.zeros();
	bool useTRIM = true;
	upPARS = upPARS_0;
	upPARS.subvec(GI.at(3,0),GI.at(3,1)).zeros(); // upETA = 0
	upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros(); // upALPHA = 0
	upPARS = Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,
		PHASE*0,SNP,RHO,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,
		PARS,useASREC,useTRIM,max_iter,eps,gr_eps,conv_eps,show);
	
	// Update upKAPPA/upETA/upALPHA
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upETA = upKAPPA % upPARS_0.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upETA % upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
	
	// Check eta without ASREC
	if( arma::any( upETA == 1.0) ){
		if( show ) Rcpp::Rcout << "### Check ETA w/o ASREC ...\n";
		upPARS.subvec(GI.at(3,0),GI.at(3,1)) = upETA;
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		PARS.zeros();
		upPARS = Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,
			PHASE*0,SNP,RHO,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,PARS,
			useASREC,useTRIM,max_iter,eps,gr_eps,conv_eps,show);
	}
	
	// Update upKAPPA/upETA/upALPHA
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upETA = upKAPPA % upPARS.subvec(GI.at(3,0),GI.at(3,1))
		% upPARS_0.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upETA % upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
	
	if( arma::all(upETA == 0) ) PHASE_available = false;
	useASREC = PHASE_available;
	
	// Check eta with ASREC
	if( PHASE_available && arma::any( upETA == 1.0 ) ){
		if( show ) Rcpp::Rcout << "### Check ETA w/ ASREC ...\n";
		useASREC = PHASE_available;
		upPARS.subvec(GI.at(3,0),GI.at(3,1)) = upETA;
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		PARS.zeros();
		upPARS = Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,
			PHASE,SNP,RHO,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,PARS,useASREC,
			useTRIM,max_iter,eps,gr_eps,conv_eps,show);
		
		// Update upKAPPA/upETA/upALPHA
		upKAPPA.at(0) = upPARS.at(0);
		if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
		upETA = upKAPPA % upPARS.subvec(GI.at(3,0),GI.at(3,1)) 
			% upPARS_0.subvec(GI.at(3,0),GI.at(3,1));
		upALPHA = upETA % upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
		
		// Check alpha
		if( useASREC && arma::any( upALPHA == 1.0 ) ){
			if( show ) Rcpp::Rcout << "### Check ALPHA ...\n";
			upPARS.subvec(GI.at(3,0),GI.at(3,1)) = upETA;
			upPARS.subvec(GI.at(5,0),GI.at(5,1)) = upALPHA;
			PARS.zeros();
			upPARS = Rcpp_CSeQTL_fullOptTest(TREC,LGX1,hap2,ASREC,LBC,
				PHASE,SNP,RHO,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,
				PARS,useASREC,useTRIM,max_iter,eps,gr_eps,conv_eps,show);
		}
		
	}
	
	// Update upKAPPA/upETA/upALPHA
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upETA = upKAPPA % upPARS.subvec(GI.at(3,0),GI.at(3,1)) 
		% upPARS_0.subvec(GI.at(3,0),GI.at(3,1));
	upALPHA = upETA % upPARS.subvec(GI.at(5,0),GI.at(5,1))
		% upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
	
	if( show ){
		Rcpp::Rcout << "### Model Check output \n";
		Rcpp::Rcout << "\tupPHI = " << upPARS.at(GI.at(1,0)) << "\n";
		Rcpp::Rcout << "\tupKAPPA = " << upKAPPA.t();
		Rcpp::Rcout << "\tupETA = " << upETA.t();
		Rcpp::Rcout << "\tupPSI = " << upPARS.at(GI.at(4,0)) << "\n";
		Rcpp::Rcout << "\tupALPHA = " << upALPHA.t();
		Rcpp::Rcout << "useASREC = ";
		if( useASREC ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
	}
	
	return upPARS;
	
}

arma::vec Rcpp_CSeQTL_swapCT(const arma::vec& TREC,const arma::vec& LGX1,
	arma::mat& RHO_copy,arma::uword& swap_CT,const arma::mat& XX,
	const arma::umat& GI,const arma::mat& I_np,const arma::uword& QQ,
	const arma::uword& np,const arma::vec& iBETA,const double& iPHI,
	const arma::vec& upPARS_0,bool& useTRIM,const arma::uword& max_iter,
	const double& eps,const double& gr_eps,const double& conv_eps,
	const bool& show){
	
	if( show ) Rcpp::Rcout << "###------- Run swapCT() ...\n";
	swap_CT = 0;
	arma::uword NN = TREC.n_elem;
	arma::vec PARS = arma::zeros<arma::vec>(np),
		mu_total = arma::zeros<arma::vec>(QQ),
		vz = arma::zeros<arma::vec>(NN),
		upKAPPA = mu_total, upPARS = upPARS_0;
	arma::uvec uvz = arma::zeros<arma::uvec>(NN);
	bool useASREC = false;
	
	if( QQ == 1 ){;
		return upPARS;
	}
	
	bool useTRIM_orig = useTRIM;
	while( true ){
		// Run full model checking
		if( show ){
			Rcpp::Rcout << "\tswapCT(): Run fullOptTest() ...\n";
		}
		PARS.zeros();
		upPARS = upPARS_0;
		upPARS.subvec(GI.at(3,0),GI.at(3,1)).zeros(); // Don't estimate ETA
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros(); // Don't estimate ALPHA
		upPARS = Rcpp_CSeQTL_fullOptTest(TREC,LGX1,vz,vz,vz,
			uvz,uvz,RHO_copy,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,
			PARS,useASREC,useTRIM,max_iter,eps,gr_eps,conv_eps,show);
		
		upKAPPA.at(0) = upPARS.at(0);
		if( QQ > 1 ){
			upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
		}
		
		if( useTRIM_orig == true && useTRIM == false ){
			return upPARS_0;
		}
		
		// Calculate MU total
		mu_total = arma::sum(arma::trans(Rcpp_CSeQTL_MU(GI,PARS)),1);
		if( show ) mu_total.t().print("\tswapCT(): MU_total = ");
		
		// Swap the original reference cell type with cell type with largest baseline total expression
		arma::uvec express_ct = arma::find( mu_total == arma::max(mu_total) );
		if( show ) Rcpp::Rcout << "\tswapCT(): upKAPPA = " << upKAPPA.t();
		// if( show ) Rcpp::Rcout << "\tswapCT = " << express_ct.at(0) + 1 << "\n";
		if( express_ct.at(0) == 0 || upKAPPA.at(0) == 1.0 ){
			// Break when the cell type w/ largest TReC is set to reference
			if( show ) Rcpp::Rcout << "\tswapCT(): Done ...\n";
			break;
		} else {
			if( show ) Rcpp::Rcout << "\tswapCT(): swap cell type 1 with cell type " 
				<< express_ct.at(0) + 1 << "\n";
			swap_CT = express_ct.at(0);
			RHO_copy.swap_cols(0,swap_CT);
		}
		
	}
	
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upPARS.subvec(GI.at(3,0),GI.at(3,1)) = upKAPPA;
	upPARS.subvec(GI.at(5,0),GI.at(5,1)) = upKAPPA;
	
	return upPARS;
}

void Rcpp_CSeQTL_TEST(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::vec& PARS_BPK,const arma::umat& GI,
	const arma::mat& I_np,const arma::vec& upPARS_full,const arma::uword& QQ,
	arma::vec& PARS_full,arma::vec& ETA,double& LL_out,arma::mat& LRT_trec,
	arma::mat& LRT_trecase,arma::mat& LRT_cistrans,const double& cistrans,
	const arma::uword& max_iter,const double& eps,const double& gr_eps,
	const double& conv_eps,const bool& show){
	
	// Constants
	arma::uword converge = 0, ct;
	arma::vec PARS_tmp = PARS_BPK, PARS_trec = PARS_BPK,
		upPARS = upPARS_full;
	double LL_trec, LL_full, larger_LL, reduced_LL, tmp_ETA = 0.0;
	if( show ) Rcpp_CSeQTL_upPARS(QQ,GI,upPARS_full);
	bool yes_PHASE = arma::any(PHASE == 1.0);
	
	// Get full TReC model
	if( show ) Rcpp::Rcout << "Get PARS_trec ...\n";
	upPARS = upPARS_full;
	upPARS.subvec(GI.at(4,0),GI.at(5,1)).zeros();
	PARS_trec = PARS_BPK;
	Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
		SNP,RHO,XX,GI,PARS_trec,converge,I_np,upPARS,
		max_iter,eps,gr_eps,conv_eps,false);
	LL_trec = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
		LBC,PHASE*0,SNP,RHO,XX,GI,PARS_trec);
	ETA = arma::exp(PARS_trec.subvec(GI.at(3,0),GI.at(3,1)));
	LL_out = LL_trec;
	
	// Get full TReC+ASE model
	if( yes_PHASE ){
		if( show ) Rcpp::Rcout << "Get PARS_full ...\n";
		upPARS = upPARS_full;
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		PARS_full = PARS_trec;
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
			SNP,RHO,XX,GI,PARS_full,converge,I_np,upPARS,
			max_iter,eps,gr_eps,conv_eps,false);
		upPARS = upPARS_full;
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
			SNP,RHO,XX,GI,PARS_full,converge,I_np,upPARS,
			max_iter,eps,gr_eps,conv_eps,false);
		LL_full = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
			LBC,PHASE,SNP,RHO,XX,GI,PARS_full);
	} else {
		LL_full = LL_trec;
		PARS_full = PARS_trec;
	}
	
	// Get null model MLEs
	for(ct = 0; ct < QQ; ct++){
		if( show ) Rcpp::Rcout << "# --- Test ct = " << ct + 1 << "...\n";
		
		// TReC-only
		larger_LL = LL_trec;
		upPARS = upPARS_full;
		upPARS.subvec(GI.at(4,0),GI.at(5,1)).zeros();
		if( upPARS.at(GI.at(3,0) + ct) == 0.0 ){
			reduced_LL = larger_LL;
		} else {
			upPARS.at(GI.at(3,0) + ct) = 0.0;
			PARS_tmp = PARS_trec;
			Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
				SNP,RHO,XX,GI,PARS_tmp,converge,I_np,upPARS,
				max_iter,eps,gr_eps,conv_eps,false);
			reduced_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
				LBC,PHASE*0,SNP,RHO,XX,GI,PARS_tmp);
		}
		
		LRT_trec.at(ct,0) = 2.0 * (larger_LL - reduced_LL);
		if( LRT_trec.at(ct,0) < 0 ) LRT_trec.at(ct,0) = 0.0;
		LRT_trec.at(ct,1) = 1.0;
		LRT_trec.at(ct,2) = 1.0 - R::pchisq(LRT_trec.at(ct,0),LRT_trec.at(ct,1),1,0);
		if( show ) Rcpp::Rcout << "\t TReC-only: Larger LL = "
			<< larger_LL << "; Reduced LL = " << reduced_LL << "\n";
		
		if( !yes_PHASE ){
			LRT_trecase.row(ct) = LRT_trec.row(ct);
			LRT_cistrans.at(ct,0) = 40.0;
			LRT_cistrans.at(ct,1) = 1.0;
			LRT_cistrans.at(ct,2) = 1.0 - R::pchisq(LRT_cistrans.at(ct,0),
				LRT_cistrans.at(ct,1),1,0);
		} else {
		
			// TReCASE (assume ct-th cell type is cis)
			upPARS = upPARS_full;
			PARS_tmp = PARS_full;
			if( upPARS.at(GI.at(5,0) + ct) == 0.0 ){
				larger_LL = LL_full;
			} else {
				upPARS.at(GI.at(5,0) + ct) = 0.0;
				Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
					SNP,RHO,XX,GI,PARS_tmp,converge,I_np,upPARS,
					max_iter,eps,gr_eps,conv_eps,false);
				larger_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
					LBC,PHASE,SNP,RHO,XX,GI,PARS_tmp);
			}
			tmp_ETA = std::exp(PARS_tmp.at(GI.at(3,0) + ct));
			
			LRT_cistrans.at(ct,0) = 2.0 * (LL_full - larger_LL);
			if( LRT_cistrans.at(ct,0) < 0 ) LRT_cistrans.at(ct,0) = 0.0;
			LRT_cistrans.at(ct,1) = 1.0;
			LRT_cistrans.at(ct,2) = 1.0 - R::pchisq(LRT_cistrans.at(ct,0),
				LRT_cistrans.at(ct,1),1,0);
			if( show ) Rcpp::Rcout << "\t Cis-Trans: Larger LL = "
				<< LL_full << "; Reduced LL = " << larger_LL << "\n";
			
			if( upPARS.at(GI.at(3,0) + ct) == 0.0 ){
				reduced_LL = larger_LL;
			} else {
				upPARS.at(GI.at(3,0) + ct) = 0.0;
				Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
					SNP,RHO,XX,GI,PARS_tmp,converge,I_np,upPARS,
					max_iter,eps,gr_eps,conv_eps,false);
				reduced_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
					LBC,PHASE,SNP,RHO,XX,GI,PARS_tmp);
			}
			
			LRT_trecase.at(ct,0) = 2.0 * (larger_LL - reduced_LL);
			if( LRT_trecase.at(ct,0) < 0 ) LRT_trecase.at(ct,0) = 0.0;
			LRT_trecase.at(ct,1) = 1.0;
			LRT_trecase.at(ct,2) = 1.0 - R::pchisq(LRT_trecase.at(ct,0),
				LRT_trecase.at(ct,1),1,0);
			if( show ) Rcpp::Rcout << "\t TReCASE: Larger LL = "
				<< larger_LL << "; Reduced LL = " << reduced_LL << "\n";
		
		}
		
		if( LRT_cistrans.at(ct,2) >= cistrans ){
			ETA.at(ct) = tmp_ETA;
		}
		
	}
	
	if( show ){
		LRT_trec.print("LRT_trec = ");
		LRT_trecase.print("LRT_trecase = ");
		LRT_cistrans.print("LRT_cistrans = ");
		Rcpp::Rcout << "ETA = " << ETA.t();
	}
	
}


// --------------------
// Cook's Distance, Optimization, Hypothesis testing Functions

arma::vec Rcpp_CSeQTL_cookD_once(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::mat& RHO,const arma::mat& XX,const arma::umat& GI,
	const arma::uword& np,const arma::mat& I_np,const arma::vec& upPARS,
	const arma::vec& iBETA,const double& iPHI,const arma::uword& max_iter,
	const double& eps,const double& gr_eps,const double& conv_eps,
	const bool& show,const int& ncores){
	
	// Note: SNP_0 should be a uvector of zeros b/c we're ignoring SNPs!!
	
	arma::uword NN = TREC.n_elem, PP = XX.n_cols, QQ = RHO.n_cols, ii,
		converge = 0;
	arma::vec PARS = arma::zeros<arma::vec>(np),
		vz = arma::zeros<arma::vec>(NN),
		cookd = vz,full_MU = vz,
		vec_zeros1 = arma::zeros<arma::vec>(NN-1);
	arma::uvec uvz = arma::zeros<arma::uvec>(NN),
		 seq_idx = uvz, SNP = uvz,
		 uvec_zeros1 = arma::zeros<arma::uvec>(NN-1);
	double PHI;
	
	// TReC optimization
	Rcpp_NB_iPARS(iBETA,iPHI,GI,PARS);
	vz.zeros();
	Rcpp_CSeQTL_BFGS(TREC,LGX1,vz,vz,vz,uvz,uvz,
		RHO,XX,GI,PARS,converge,I_np,upPARS,
		max_iter,eps,gr_eps,conv_eps,show);
	arma::vec BETA = PARS.subvec(GI.at(0,0),GI.at(0,1)),
		KAPPA = arma::ones<arma::vec>(QQ);
	if(QQ > 1) KAPPA.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
	PHI = std::exp(PARS.at(1,0));
	
	// Calculate full_MU, aka log(MU) estimated from all subjects
	arma::vec ETA_KAP_TREC = KAPPA % arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1)));
	full_MU = XX * BETA + arma::log(RHO * KAPPA);
	full_MU = arma::exp(full_MU);
	arma::vec var_subj = full_MU + PHI * full_MU % full_MU;
	
	for(ii = 0; ii < NN; ii++){
		seq_idx.at(ii) = ii;
	}
	
	// assuming the number of parameters is covariates + baseline expression
	arma::uword pp = PP + QQ - 1; 
	
	// Get model predictions after removing i-th subject
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(ncores,show,NN,pp,np,I_np,GI,QQ,seq_idx,cookd,\
		max_iter,eps,gr_eps,conv_eps,TREC,LGX1,vec_zeros1,uvec_zeros1,XX,RHO,\
		iBETA,iPHI,upPARS,var_subj,full_MU)
	#endif
	for(arma::uword ii2 = 0; ii2 < NN; ii2++){
		if( ncores == 1 && show ){
			if( (ii2 + 1) % 5 == 0 ) Rcpp::Rcout << ".";
			if( (ii2 + 1) % 100 == 0 || (ii2 + 1) == NN ){
				Rcpp::Rcout << ii2+1 << " out of " << NN << "\n";
			}
		}
		
		arma::uword converge2 = 0;
		arma::vec tmp_PARS = arma::zeros<arma::vec>(np),
			tmp_vec = arma::zeros<arma::vec>(NN);
		arma::uvec sub_idx = arma::find(seq_idx != ii2); // subset index
		
		Rcpp_NB_iPARS(iBETA,iPHI,GI,tmp_PARS);
		Rcpp_CSeQTL_BFGS(TREC(sub_idx),LGX1(sub_idx),vec_zeros1,
			vec_zeros1,vec_zeros1,uvec_zeros1,uvec_zeros1,
			RHO.rows(sub_idx.t()),XX.rows(sub_idx.t()),GI,
			tmp_PARS,converge2,I_np,upPARS,max_iter,
			eps,gr_eps,conv_eps,false);
		
		arma::vec tmp_KAPPA = arma::ones<arma::vec>(QQ);
		if(QQ > 1) tmp_KAPPA.subvec(1,QQ-1) = arma::exp(tmp_PARS.subvec(GI.at(2,0),GI.at(2,1)));
		
		tmp_vec = XX * tmp_PARS.subvec(GI.at(0,0),GI.at(0,1)) + arma::log(RHO * tmp_KAPPA); // red_logMU
		tmp_vec = arma::exp(tmp_vec);
		tmp_vec = full_MU - tmp_vec;
		tmp_vec %= tmp_vec;
		tmp_vec /= var_subj;
		cookd.at(ii2) = arma::sum(tmp_vec) / pp;
	}

	return cookd;
}

void Rcpp_CSeQTL_gs_cistrans(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::vec& PARS_trec,const arma::umat& GI,
	const arma::mat& I_np,const arma::vec& upPARS_0,const arma::uword& QQ,
	const arma::uword& np,arma::mat& LRT,const arma::uword& max_iter,
	const double& eps,const double& gr_eps,const double& conv_eps,const bool& show){
	
	if( show ) Rcpp::Rcout << "\n############\n#### Cis/Trans Section\n############\n";
	
	// LRT = (LRT statistic,DF,pvalue)
	LRT.zeros();
	arma::uword converge = 0, ct;
	
	// If no phasing available ...
	if( arma::all(PHASE == 0) ){
		for(ct = 0; ct < QQ; ct++){
			LRT.at(ct,1) = 1;
		}
		return;
	}
	
	arma::vec init_PARS = PARS_trec, full_PARS = init_PARS,
		upPARS = upPARS_0, null_PARS = init_PARS, 
		upALPHA_red = arma::zeros<arma::vec>(QQ);
	
	// 1st: Start with PARS_trec (BETA,PHI,KAPPA,ETA) w/o ASREC
	// 2nd: Optimize over (BETA,PHI,KAPPA,ETA,PSI) w/ ASREC, set all ALPHA = 1
	if( show ) Rcpp::Rcout << "\n####### Optimize BETA,PHI,KAPPA,ETA,PSI (set ALPHA=1)!\n";
	upPARS = upPARS_0;
	upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
		SNP,RHO,XX,GI,init_PARS,converge,I_np,upPARS,
		max_iter,eps,gr_eps,conv_eps,show);
	
	// Set PARS
	null_PARS = init_PARS; full_PARS = init_PARS;
	
	// Get full model MLEs
	if( show ) Rcpp::Rcout << "\n####### Optimize BETA,PHI,KAPPA,ETA,PSI,ALPHA!\n";
	Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
		SNP,RHO,XX,GI,full_PARS,converge,I_np,upPARS_0,
		max_iter,eps,gr_eps,conv_eps,show);
	double full_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,
		ASREC,LBC,PHASE,SNP,RHO,XX,GI,full_PARS);

	// Get null model MLEs
	double null_LL;
	for(ct = 0; ct < QQ; ct++){
		if( show ) Rcpp::Rcout << "#---------- Test model ct = " << ct + 1 << "\n";
		
		null_PARS = init_PARS; // reset PARS
		upALPHA_red = upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
		upALPHA_red.at(ct) = 0.0;
		
		if( arma::all(upALPHA_red == upPARS_0.subvec(GI.at(5,0),GI.at(5,1))) ){
			// aka if we can't test haplotype counts (aka estimate ALPHA), we'll report the TReC result
			// LRT.at(ct,0) = 40.0;
			
			// aka if we can't test haplotype counts (aka estimate ALPHA), we'll report the TReCASE result
			LRT.at(ct,0) = 0.0;
			
			LRT.at(ct,1) = 1.0;
			LRT.at(ct,2) = 1.0 - R::pchisq(LRT.at(ct,0),LRT.at(ct,1),1,0);
			continue;
		}
		
		// Optimize null model, H_0: alpha_ct = 1
		upPARS = upPARS_0;
		upPARS.subvec(GI.at(5,0),GI.at(5,1)) = upALPHA_red;
		Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE,
			SNP,RHO,XX,GI,null_PARS,converge,I_np,upPARS,
			max_iter,eps,gr_eps,conv_eps,show);
		
		// Calculate output values
		null_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,
			LBC,PHASE,SNP,RHO,XX,GI,null_PARS);
		LRT.at(ct,0) = 2.0 * (full_LL - null_LL);
		if( LRT.at(ct,0) < 0 ) LRT.at(ct,0) = 0.0;
		LRT.at(ct,1) = arma::sum(upPARS_0.subvec(GI.at(5,0),GI.at(5,1)) - upALPHA_red);
		LRT.at(ct,2) = 1.0 - R::pchisq(LRT.at(ct,0),LRT.at(ct,1),1,0);
		
	}
	
	if( show ) LRT.print("Cis/Trans LRT results:");
}

void Rcpp_CSeQTL_genesnp(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,const arma::mat& I_np,
	const arma::uword& QQ,const arma::uword& np,const arma::vec& iBETA,
	const double& iPHI,const arma::vec& upPARS_0,double& LL_full,
	arma::vec& PARS_full,arma::vec& ETA,arma::mat& LRT_trec,
	arma::mat& LRT_trecase,arma::mat& LRT_cistrans,arma::mat& MU_full,
	const double& cistrans,const arma::uword& max_iter,const double& eps,
	const double& gr_eps,const double& conv_eps,const bool& show){
	
	arma::uword converge = 0;
	arma::vec PARS_BPK = arma::zeros<arma::vec>(np),
		upPARS = upPARS_0, upPARS_2 = upPARS;
	
	// Run full model checking, which parameters to estimate and 
	//	optimize over and should we use ASREC
	bool useASREC = arma::any(PHASE == 1);
	upPARS = Rcpp_CSeQTL_run_fullModel(TREC,LGX1,hap2,ASREC,LBC,
		PHASE,SNP,RHO,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,useASREC,
		max_iter,eps,gr_eps,conv_eps,show);
		// aka if model checking shows ASREC + TREC don't agree 
		// with each other, use TREC only, no cistrans or trecase analyses
	
	// Get PARS_BPK after model checking
	upPARS_2 = upPARS;
	upPARS_2.subvec(GI.at(3,0),GI.at(5,1)).zeros();
	Rcpp_NB_iPARS(iBETA,iPHI,GI,PARS_BPK);
	Rcpp_CSeQTL_BFGS(TREC,LGX1,hap2,ASREC,LBC,PHASE*0,
		SNP,RHO,XX,GI,PARS_BPK,converge,I_np,upPARS_2,
		max_iter,eps,gr_eps,conv_eps,show);
	
	// Hypothesis Testing: TReC, TReCASE, Cis/Trans
	Rcpp_CSeQTL_TEST(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,
		PARS_BPK,GI,I_np,upPARS,QQ,PARS_full,ETA,LL_full,LRT_trec,
		LRT_trecase,LRT_cistrans,cistrans,max_iter,eps,
		gr_eps,conv_eps,show);
	
	MU_full = Rcpp_CSeQTL_MU(GI,PARS_full);
}

// [[Rcpp::export]]
Rcpp::List Rcpp_CSeQTL_FDR(const arma::vec& pvalues,
	const double& lambda = 0.5,const int& ncores = 1){
	
	double prop_one, prop_null, np, max_FDR, alpha,
		FDR, num_sig, max_pvalue;
	arma::vec uniq_pvalues = arma::sort(arma::unique(pvalues));
	arma::uword num_test = pvalues.n_elem,
		num_uniq = uniq_pvalues.n_elem, bb;
	arma::vec qvalues = arma::zeros<arma::vec>(num_test);
	
	// Calculate some constants
	prop_one = arma::sum( pvalues == 1.0 ) * 1.0 / num_test;
	max_pvalue = arma::max(pvalues);
	prop_null = 1.0 / (1.0 - lambda) * 
		arma::sum(pvalues >= lambda && pvalues < 1.0) * 
		1.0 / num_test + prop_one;
		if( prop_null > 1.0 ) prop_null = 1.0;
	np = num_test * prop_null;
	max_FDR = max_pvalue * prop_null;
	
	arma::mat res = arma::zeros<arma::mat>(num_uniq,4);
		// (alpha,FDR,sig_tests,exp_false_pos)
		res.col(0) = uniq_pvalues;
	
	for(arma::uword aa = 0; aa < num_uniq; aa++){
		bb = num_uniq - 1 - aa;
		alpha = res.at(bb,0);
		num_sig = arma::sum(pvalues <= alpha);
		FDR = np * alpha / num_sig;
		if( bb == num_uniq - 1 ){
			FDR = max_FDR;
		} else {
			if( FDR > res.at(bb+1,1) ) FDR = res.at(bb+1,1);
		}
		res.at(bb,1) = FDR;
	}
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(static) \
		num_threads(ncores) \
		shared(res,pvalues,num_test,qvalues)
	#endif
	for(arma::uword aa = 0; aa < num_test; aa++){
		arma::mat sub_res = res.rows(arma::find(res.col(0) >= pvalues.at(aa)));
		qvalues.at(aa) = arma::min(sub_res.col(1));
	}
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(static) \
		num_threads(ncores) \
		shared(res,num_uniq,qvalues)
	#endif
	for(arma::uword aa = 0; aa < num_uniq; aa++){
		res.at(aa,2) = arma::sum(qvalues <= res.at(aa,1)); // num sig_tests
	}
	res.col(3) = num_test * res.col(1);
	
	return Rcpp::List::create(
		Rcpp::Named("prop_one",prop_one),
		Rcpp::Named("prop_null",prop_null),
		Rcpp::Named("qvalues",
			Rcpp::NumericVector(qvalues.begin(),qvalues.end())),
		Rcpp::Named("res",res));
}


// --------------------
// CSeQTL profile likelihood related functions

arma::vec Rcpp_CSeQTL_calc_pGRAD(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::umat& GI,
	const arma::vec& upPARS,const arma::vec& PARS,const arma::uword& index){
	
	arma::vec GRAD = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,
		PHASE,SNP,RHO,XX,GI,upPARS,PARS);
	GRAD.at(index) = 0.0;
	return GRAD;
}

void Rcpp_CSeQTL_pBFGS(const arma::vec& TREC,const arma::vec& LGX1,
	const arma::vec& hap2,const arma::vec& ASREC,const arma::vec& LBC,
	const arma::uvec& PHASE,const arma::uvec& SNP,const arma::mat& RHO,
	const arma::mat& XX,const arma::umat& GI,arma::vec& PARS,
	arma::uword& converge,const arma::mat& I_np,const arma::vec& upPARS,
	const arma::uword& index,const arma::uword& max_iter,
	const double& eps,const bool& show){
	
	converge = 0;
	arma::uword iter = 0, jj, uu,
		reset_Bk = 0, np = PARS.n_elem, QQ = RHO.n_cols;
	
	// Initialize parameters
	Rcpp_CSeQTL_control_PAR(PARS,QQ,GI,upPARS);
	// if( show ) Rcpp::Rcout << "\t##\n\tiPARS = " << PARS.t();
	if( show ) Rcpp::Rcout << "iPARS = " << PARS.t();
	
	arma::mat inv_Bk = I_np, ISYT = inv_Bk;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk, gr_k = curr_xk,
		p_k = curr_xk, s_k = curr_xk, y_k = curr_xk,
		ETA_ASREC = arma::zeros<arma::vec>(QQ);
	double old_LL, new_LL, inv_norm_p_k, tmp_alpha, ys,
		fnscale = -1.0, curr_LL = 0.0;
	
	while(iter < max_iter){
		// Calculate Direction p_k
		gr_k = fnscale * Rcpp_CSeQTL_calc_pGRAD(TREC,hap2,ASREC,PHASE,
			SNP,RHO,XX,GI,upPARS,xk,index);
		
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));

		// Line search for new xk
		uu = 0;
		old_LL = fnscale * Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,
			ASREC,LBC,PHASE,SNP,RHO,XX,GI,xk);
		
		if(old_LL < 0) break; // Less than 0 b/c fnscale = -1
		
		for(jj = 0; jj <= 30; jj++){
			tmp_alpha = inv_norm_p_k / std::pow(4,jj);
			new_xk = xk + tmp_alpha * p_k;
			new_LL = fnscale * Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,
				ASREC,LBC,PHASE,SNP,RHO,XX,GI,new_xk);
			if(new_LL < old_LL){ // minimizing
				s_k = tmp_alpha * p_k;
				y_k = fnscale * Rcpp_CSeQTL_calc_pGRAD(TREC,hap2,ASREC,PHASE,
					SNP,RHO,XX,GI,upPARS,new_xk,index) - gr_k;
				ys = arma::dot(y_k,s_k);
				if( ys > 0.0 ){
					ISYT = I_np - (s_k * y_k.t()) / ys;
					inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
				}
				xk = new_xk;
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if(show) printR_obj("Reset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if(show) printR_obj("Failed line search");
				break;
			}
		}
		
		if( reset_Bk > 5 ) break;
		
		// Check Convergence
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				Rcpp_norm(curr_xk - xk) < eps ){
				gr_k = Rcpp_CSeQTL_calc_pGRAD(TREC,hap2,ASREC,PHASE,
					SNP,RHO,XX,GI,upPARS,xk,index);
				if( Rcpp_norm(gr_k) < eps ){
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	}
	
	// Update parameters
	PARS = xk;
	
	if( show ){
		// Calculate LL, GRAD
		old_LL = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,RHO,XX,GI,PARS);
		gr_k = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,GI,upPARS,PARS);
		
		Rcpp::Rcout << "\tIter = " << iter+1 
			<< "; LL = " << old_LL << "\n";
		Rcpp::Rcout << "\tPARS = " << PARS.t();
		Rcpp::Rcout << "\tConvergence Indicators: \n"
			<< "\t   NormGrad = " << Rcpp_norm(gr_k)
			<< "; NormIBkGrad = " << Rcpp_norm(inv_Bk * gr_k) << "\n";
		Rcpp::Rcout << "\tGRAD = " << gr_k.t();
		arma::vec KAPPA = arma::ones<arma::vec>(QQ), upKAPPA = KAPPA;
		upKAPPA.at(0) = upPARS.at(0);
		if(QQ > 1){
			KAPPA.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
			upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
		}
		Rcpp::Rcout << "\tupKAPPA = " << upKAPPA.t();
		Rcpp::Rcout << "\tupETA = " << upPARS.subvec(GI.at(3,0),GI.at(3,1)).t();
		Rcpp::Rcout << "\tBETA = " << PARS.subvec(GI.at(0,0),GI.at(0,1)).t();
		Rcpp::Rcout << "\tPHI = " << std::exp(PARS.at(GI.at(1,0))) << "\n";
		Rcpp::Rcout << "\tKAPPA = " << KAPPA.t();
		Rcpp::Rcout << "\tETA = " << arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))).t();
		Rcpp::Rcout << "\tPSI = " << std::exp(PARS.at(GI.at(4,0))) << "\n";
		Rcpp::Rcout << "\tALPHA = " << arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1)).t());
		Rcpp::Rcout << "Convergence Status = ";
		if( converge == 1 ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
		arma::mat MU = Rcpp_CSeQTL_MU(GI,PARS);
		MU.print("MU = ");
	}
	
}

// [[Rcpp::export]]
arma::vec Rcpp_CSeQTL_profile_PAR(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& upPARS,
	const arma::uword& index,const arma::vec& bounds,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-10,const int& ncores = 1,const bool& show = true){
	
	// Constants
	arma::uword jj, NN = TREC.n_elem, PP = XX.n_cols, QQ = RHO.n_cols,
		nvalues = bounds.n_elem;
	arma::vec LGX1 = arma::lgamma(1.0 + TREC),
		LBC = arma::zeros<arma::vec>(NN), OO = LBC;
	arma::umat GI = Rcpp_calc_GI(PP,QQ);
	arma::uword np = GI.at(5,1)+1;
	arma::mat I_np = arma::eye<arma::mat>(np,np),
		mat_PARS = arma::zeros<arma::mat>(nvalues,np);
	for(jj = 0; jj < NN; jj++){
		LBC.at(jj) = R::lchoose( ASREC.at(jj), hap2.at(jj) );
	}
	arma::vec nb_PARS = Rcpp_NB_reg_one(TREC,XX,OO,max_iter,eps,false),
		iBETA = nb_PARS.subvec(0,PP - 1);
	double iPHI = std::exp(nb_PARS.at(PP));
	arma::vec iPARS = arma::zeros<arma::vec>(np),
		vec_LL = arma::zeros<arma::vec>(nvalues);
	Rcpp_NB_iPARS(iBETA,iPHI,GI,iPARS);
	
	// arma::uword thres = 0;
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(ncores,show,I_np,upPARS,\
		GI,TREC,LGX1,hap2,ASREC,LBC,PHASE,SNP,XX,RHO,\
		index,max_iter,eps,bounds,nvalues,vec_LL,iPARS,mat_PARS)
	#endif
	for(arma::uword ii = 0; ii < nvalues; ii++){
		if( ncores == 1 && show ){ // && (ii + 0.0) / bounds.n_elem * 100.0 >= thres){
			Rcpp::Rcout << ".";
			if( (ii + 1) == bounds.n_elem ) Rcpp::Rcout << "\n";
		}
		
		arma::vec PARS = iPARS,
			upPARS_2 = upPARS;
		arma::uword converge = 0;
		
		// Step by step optimization
		PARS.at(index) = bounds.at(ii);
		upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(3,0),GI.at(3,1)).zeros();
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		Rcpp_CSeQTL_pBFGS(TREC,LGX1,hap2,ASREC,LBC,
			PHASE * 0,SNP,RHO,XX,GI,PARS,converge,I_np,
			upPARS_2,index,max_iter,eps,false);
		
		PARS.at(index) = bounds.at(ii);
		upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		Rcpp_CSeQTL_pBFGS(TREC,LGX1,hap2,ASREC,LBC,
			PHASE * 0,SNP,RHO,XX,GI,PARS,converge,I_np,
			upPARS_2,index,max_iter,eps,false);
		
		PARS.at(index) = bounds.at(ii);
		upPARS_2 = upPARS;
			upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		Rcpp_CSeQTL_pBFGS(TREC,LGX1,hap2,ASREC,LBC,
			PHASE,SNP,RHO,XX,GI,PARS,converge,I_np,
			upPARS_2,index,max_iter,eps,false);
		
		PARS.at(index) = bounds.at(ii);
		Rcpp_CSeQTL_pBFGS(TREC,LGX1,hap2,ASREC,LBC,
			PHASE,SNP,RHO,XX,GI,PARS,converge,I_np,
			upPARS,index,max_iter,eps,false);
		
		vec_LL.at(ii) = Rcpp_CSeQTL_calc_LL(TREC,LGX1,hap2,ASREC,LBC,
			PHASE,SNP,RHO,XX,GI,PARS);
		mat_PARS.row(ii) = PARS.t();
	}
	
	// Output
	return vec_LL;
}

// [[Rcpp::export]]
arma::mat Rcpp_CSeQTL_profile_PARS(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& upPARS,
	const arma::vec& bounds,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-10,const int& ncores = 1,const bool& show = true){
	
	// Constants
	arma::uword index, PP = XX.n_cols, QQ = RHO.n_cols,
		nvalues = bounds.n_elem;
	arma::umat GI = Rcpp_calc_GI(PP,QQ);
	arma::uword np = GI.at(5,1) + 1;
	arma::mat mat_LL = arma::zeros<arma::mat>(nvalues,np);
	
	for(index = 0; index < np; index++){
		if( upPARS.at(index) == 1 ){
			if( show ) Rcpp::Rcout << "PAR[" << index+1 << "]";
			mat_LL.col(index) = Rcpp_CSeQTL_profile_PAR(TREC,hap2,ASREC,
				PHASE,SNP,RHO,XX,upPARS,index,bounds,max_iter,eps,ncores,show);
		}
	}
	
	return mat_LL;
}


// --------------------
// Main CSeQTL workflow functions

// [[Rcpp::export]]
arma::mat Rcpp_CSeQTL_cooksD(const arma::vec& TREC,const arma::mat& RHO,
	const arma::mat& XX,const double& trim_thres,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-10,const double& mad_const = 1.4826,const int& ncores = 1,
	const double& gr_eps = 1e-2,const double& conv_eps = 5e-5,const bool& show = true){
	
	// Calculate constants
	arma::uword NN = XX.n_rows, PP = XX.n_cols, QQ = RHO.n_cols, converge = 0;
	arma::vec LGX1 = arma::lgamma(TREC + 1.0),
		KAPPA = arma::ones<arma::vec>(QQ),
		vec_zeros = arma::zeros<arma::vec>(NN);
	arma::umat GI = Rcpp_calc_GI(PP,QQ);
	arma::uword np = GI.at(5,1) + 1;
	arma::uvec uvec_zeros = arma::zeros<arma::uvec>(NN);
	arma::mat RHO_copy = RHO, I_np = arma::eye<arma::mat>(np,np);
	
	// Get negative binomial initialized parameters
	if( show ) Rcpp::Rcout << "   cooksD: Initiate NB pars ...\n";
	arma::vec iPARS = Rcpp_NB_reg_one(TREC,XX,vec_zeros,max_iter,eps,show),
		iBETA = iPARS.subvec(0,PP-1),
		PARS = arma::zeros<arma::vec>(np),
		upPARS = PARS;
	double iPHI = std::exp(iPARS.at(PP));
	
	// Check if swapping reference cell type is needed
	arma::uword swap_CT = 0;
	bool useTRIM = false;
	upPARS.ones();
	upPARS = Rcpp_CSeQTL_swapCT(TREC,LGX1,RHO_copy,swap_CT,
		XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,useTRIM,max_iter,
		eps,gr_eps,conv_eps,show);
	
	// Calculate Cook's Distance
	bool show2 = show && ncores == 1;
	arma::vec cooksd = Rcpp_CSeQTL_cookD_once(TREC,LGX1,RHO_copy,XX,GI,
		np,I_np,upPARS,iBETA,iPHI,max_iter,eps,gr_eps,conv_eps,show2,ncores);
	
	// Calculate parameters
	Rcpp_NB_iPARS(iBETA,iPHI,GI,PARS);
	Rcpp_CSeQTL_optPARAMS(TREC,LGX1,vec_zeros,vec_zeros,vec_zeros,
		uvec_zeros,uvec_zeros,RHO_copy,XX,converge,QQ,GI,I_np,iBETA,iPHI,
		upPARS,PARS,false,max_iter,eps,gr_eps,conv_eps,show2);
	arma::vec BETA = PARS.subvec(GI.at(0,0),GI.at(0,1));
	if(QQ > 1) KAPPA.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
	
	// Calculate log_subj_mu
	arma::vec log_subj_mu = XX * BETA + arma::log(RHO_copy * KAPPA);
	
	// Output matrix with cooksd, pred_mu, pred_trec, trim_trec
	arma::mat out = arma::zeros<arma::mat>(NN,4);
	out.col(0) = cooksd;
	out.col(1) = arma::exp(log_subj_mu);
	out.col(2) = arma::round(out.col(1)) + 1.0;
	
	// Standardize cooksd aka (xx - median(xx)) / mad(xx)
	cooksd -= arma::median(cooksd); // xx - median(xx)
	cooksd /= mad_const * arma::median(arma::abs(cooksd)); // mad(xx)
	
	// Calculate trimmed TREC
	arma::vec trim_TREC = TREC;
	arma::uvec idx = arma::find(cooksd >= trim_thres); // isolate influential subjects
	arma::vec sub_TREC = out.col(2);
	trim_TREC(idx) = sub_TREC(idx); // impute influential subjects with predicted TReC
	out.col(3) = trim_TREC;
	
	if( show ) Rcpp::Rcout << "   Finished cooksD ...\n";
	
	return out;
}

// [[Rcpp::export]]
Rcpp::List Rcpp_CSeQTL_BFGS_smart(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE_0,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& upPARS_0,
	const bool& iFullModel = true,const bool& trim = true,const double& trim_thres = 10,
	const bool& hypotest = false,const bool& swap = true,const arma::uword& numAS = 5,
	const arma::uword& numASn = 10,const arma::uword& numAS_het = 5,
	const double& mad_const = 1.4826,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-10,const int& ncores = 1,const double& gr_eps = 1e-2,
	const double& conv_eps = 5e-5,const double& hess_shift = 1e-5,const bool& show = true){
	
	arma::uword ii, NN = XX.n_rows, PP = XX.n_cols, QQ = RHO.n_cols;
	arma::umat GI = Rcpp_calc_GI(PP,QQ);
	arma::uword np = GI.at(5,1)+1;
	arma::mat I_np = arma::eye<arma::mat>(np,np);
	arma::vec iPARS = arma::zeros<arma::vec>(PP + 1),
		iBETA = arma::zeros<arma::vec>(PP),PARS = arma::zeros<arma::vec>(np),
		OO = arma::zeros<arma::vec>(NN),LGX1 = OO, LBC = OO,
		upPARS = upPARS_0, upPARS_2 = upPARS,
		KAPPA = arma::ones<arma::vec>(QQ), upKAPPA = KAPPA;
	
	// Define upKAPPA
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	// upETA = upKAPPA % upETA_0
	upPARS.subvec(GI.at(3,0),GI.at(3,1)) %= upKAPPA;
	// upALPHA = upETA % upALPHA_0
	upPARS.subvec(GI.at(5,0),GI.at(5,1)) %= upPARS.subvec(GI.at(3,0),GI.at(3,1));
	
	arma::uvec PHASE = PHASE_0;
	PHASE(arma::find(ASREC < numAS)).zeros();
	// Count number of subjects with at least numAS ASREC
	arma::uword num_subj = arma::sum(ASREC >= numAS),hap_cnt;
	if( num_subj < numASn ){
		PHASE.zeros();
		if( show ) Rcpp::Rcout << "###------- ASREC check, insufficient counts ...\n";
	} else {
		hap_cnt = arma::sum(PHASE == 1 && ASREC >= numAS && (SNP == 1 || SNP == 2));
		if( hap_cnt < numAS_het ){
			PHASE.zeros();
			if( show ) Rcpp::Rcout << "###------- ASREC check, insufficient counts for SNP ...\n";
		}
	}
	
	if( arma::all(SNP == 0) ){
		upPARS.subvec(GI.at(3,0),GI.at(3,1)).zeros();
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	}
	if( arma::all(PHASE == 0) ){
		upPARS.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	}
	
	// Get LGX1 and LBC
	if( show ) Rcpp::Rcout << "###------- Get LGX1 and LBC ...\n";
	LGX1 = arma::lgamma(TREC + 1.0);
	LBC = arma::lgamma(ASREC + 1.0) - arma::lgamma(hap2 + 1.0)
		- arma::lgamma(ASREC - hap2 + 1.0);
	
	// Get iBETA and iPHI
	if( show ) Rcpp::Rcout << "###------- Get iBETA and iPHI ...\n";
	iPARS = Rcpp_NB_reg_one(TREC,XX,OO,max_iter,eps,show);
	iBETA = iPARS.subvec(0,PP-1);
	double iPHI = std::exp(iPARS.at(PP));
	PARS.subvec(0,PP) = iPARS;
	
	// trim TREC
	arma::vec final_TREC = TREC,final_LGX1 = LGX1,
		COOKSD = arma::zeros<arma::vec>(NN);
	bool useTRIM = trim;
	if( trim ){
		if( show ) Rcpp::Rcout << "###------- Trim TReC ...\n";
		arma::mat cooksd_res = Rcpp_CSeQTL_cooksD(TREC,RHO,XX,
			trim_thres,max_iter,eps,mad_const,ncores,gr_eps,
			conv_eps,show);
		final_TREC = cooksd_res.col(3);
		final_LGX1 = arma::lgamma(final_TREC + 1.0);
		iPARS = Rcpp_NB_reg_one(final_TREC,XX,OO,max_iter,eps,false);
		iBETA = iPARS.subvec(0,PP-1);
		iPHI = std::exp(iPARS.at(PP));
		PARS.subvec(0,PP) = iPARS;
		COOKSD = cooksd_res.col(0);
	}
	
	// Check reference cell type
	arma::mat RHO_copy = RHO;
	arma::uword swap_CT = 0;
	if( swap ){
		upPARS = Rcpp_CSeQTL_swapCT(final_TREC,final_LGX1,RHO_copy,
			swap_CT,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,useTRIM,max_iter,eps,
			gr_eps,conv_eps,show);
		
		// If we tried trimmed TREC but failed to converge, reset to untrimmed TREC
		if( trim == true && useTRIM == false ){
			if( show ) Rcpp::Rcout << "*** DON'T USE TRIMMED TREC, RESET to ORIG TREC!\n";
			
			// Set TREC and LGX1
			final_TREC = TREC;
			final_LGX1 = LGX1;
			
			// Get iBETA,iPHI
			iPARS = Rcpp_NB_reg_one(TREC,XX,OO,max_iter,eps,false);
			iBETA = iPARS.subvec(0,PP-1);
			iPHI = std::exp(iPARS.at(PP));
			PARS.subvec(0,PP) = iPARS;
			
			// Re-run swap
			// upKAPPA.ones();
			upPARS.at(0) = 1.0;
			if( QQ > 1 ) upPARS.subvec(GI.at(2,0),GI.at(2,1)).ones();
			RHO_copy = RHO;
			upPARS = Rcpp_CSeQTL_swapCT(final_TREC,final_LGX1,RHO_copy,
				swap_CT,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,useTRIM,max_iter,eps,
				gr_eps,conv_eps,show);
		}
	}
	// upETA = upKAPPA % upETA;
	upKAPPA.at(0) = upPARS.at(0);
	if( QQ > 1 ) upKAPPA.subvec(1,QQ-1) = upPARS.subvec(GI.at(2,0),GI.at(2,1));
	upPARS.subvec(GI.at(3,0),GI.at(3,1)) = upKAPPA 
		% upPARS_0.subvec(GI.at(3,0),GI.at(3,1));
	// upALPHA = upETA % upALPHA;
	upPARS.subvec(GI.at(5,0),GI.at(5,1)) = upPARS.subvec(GI.at(3,0),GI.at(3,1))
		% upPARS_0.subvec(GI.at(5,0),GI.at(5,1));
	
	// Run full model checking: Check ref cell type total expression and ref cell type and allele expression
	bool useASREC = arma::any(PHASE == 1);
	if( iFullModel ){
		if( show ) Rcpp::Rcout << "###------- Run full model checking ...\n";
		upPARS = Rcpp_CSeQTL_run_fullModel(final_TREC,final_LGX1,hap2,ASREC,LBC,
			PHASE,SNP,RHO_copy,XX,GI,I_np,QQ,np,iBETA,iPHI,upPARS,useASREC,
			max_iter,eps,gr_eps,conv_eps,show);
	} else {
		if( show ) Rcpp::Rcout << "###------- No full model checking ...\n";
	}
	arma::uvec final_PHASE = arma::zeros<arma::uvec>(NN);
	if( useASREC ) final_PHASE = PHASE;
	
	// Optimize
	if( show ) Rcpp::Rcout << "###------- Run formal optimization ...\n";
	arma::vec PARS_BPK = PARS, PARS_trec = PARS;
	arma::uword converge = 0;
	
	upPARS_2 = upPARS;
		upPARS_2.subvec(GI.at(3,0),GI.at(3,1)).zeros();
		upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	Rcpp_CSeQTL_optPARAMS(final_TREC,final_LGX1,hap2,ASREC,LBC,
		final_PHASE*0,SNP,RHO_copy,XX,converge,QQ,GI,I_np,iBETA,iPHI,
		upPARS_2,PARS,useASREC,max_iter,eps,gr_eps,conv_eps,show);
	PARS_BPK = PARS;
	
	upPARS_2 = upPARS;
		upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
	Rcpp_CSeQTL_optPARAMS(final_TREC,final_LGX1,hap2,ASREC,LBC,
		final_PHASE*0,SNP,RHO_copy,XX,converge,QQ,GI,I_np,iBETA,iPHI,
		upPARS_2,PARS,useASREC,max_iter,eps,gr_eps,conv_eps,show);
	PARS_trec = PARS;
	
	Rcpp_CSeQTL_optPARAMS(final_TREC,final_LGX1,hap2,ASREC,LBC,
		final_PHASE,SNP,RHO_copy,XX,converge,QQ,GI,I_np,iBETA,iPHI,
		upPARS,PARS,useASREC,max_iter,eps,gr_eps,conv_eps,show);
	// PARS is the full model
	
	// Calculate gradient and log likelihood
	arma::vec gk = Rcpp_CSeQTL_calc_GRAD(final_TREC,hap2,ASREC,
		final_PHASE,SNP,RHO_copy,XX,GI,upPARS,PARS),
		BETA = PARS.subvec(GI.at(0,0),GI.at(0,1)),
		ETA = arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))),
		ALPHA = arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1)));
	if(QQ > 1) KAPPA.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
	arma::vec ETA_KAP_TREC = ETA % KAPPA;
	double LL = Rcpp_CSeQTL_calc_LL(final_TREC,final_LGX1,hap2,ASREC,
		LBC,final_PHASE,SNP,RHO_copy,XX,GI,PARS),xi;
	
	arma::vec log_OF = arma::log(RHO_copy * KAPPA);
	arma::vec log_subj_mu = XX * BETA + log_OF;
	for(ii = 0; ii < NN; ii++){
		xi = arma::dot(RHO_copy.row(ii).t(),ETA_KAP_TREC) 
			/ arma::dot(RHO_copy.row(ii).t(),KAPPA);
		if( SNP.at(ii) == 0 ){
			continue;
		} else if( SNP.at(ii) == 1 || SNP.at(ii) == 2 ){
			log_subj_mu.at(ii) += std::log( (1.0 + xi) / 2.0 );
		} else if( SNP.at(ii) == 3 ){
			log_subj_mu.at(ii) += std::log( xi );
		} else {
			continue;
		}
	}
	arma::mat MU = Rcpp_CSeQTL_MU(GI,PARS);
	
	// Calculate approximated hessian matrix
	arma::mat hess = Rcpp_CSeQTL_calc_HESS(final_TREC,hap2,ASREC,
		final_PHASE,SNP,RHO_copy,XX,GI,PARS,upPARS,I_np,hess_shift);
	
	// Calculate variances
	arma::mat covar = arma::zeros<arma::mat>(np,np), 
		eqtl_vars = arma::zeros<arma::mat>(3,QQ);
	arma::uvec nz = arma::find( hess.diag() != 0.0 );
	double rcond_num = arma::rcond(hess.submat(nz,nz));
	if( rcond_num > 0.0 ){
		covar.submat(nz,nz) = arma::inv(-1.0 * hess.submat(nz,nz));
		// Kappa, Eta, Alpha variances
		eqtl_vars.at(0,0) = covar.at(0,0);
		if(QQ > 1) eqtl_vars(0,arma::span(1,QQ-1)) = arma::diagvec(covar.submat(GI.at(2,0),
			GI.at(2,0),GI.at(2,1),GI.at(2,1))).t();
		eqtl_vars.row(1) = arma::diagvec(covar.submat(GI.at(3,0),GI.at(3,0),GI.at(3,1),GI.at(3,1))).t();
		eqtl_vars.row(2) = arma::diagvec(covar.submat(GI.at(5,0),GI.at(5,0),GI.at(5,1),GI.at(5,1))).t();
	} else {
		Rcpp::stop("Hessian is approximately singular!");
	}
	double norm_iHG = Rcpp_norm(covar * gk);
	
	// Hypothesis Testing
	arma::mat LRT_trec = arma::zeros<arma::mat>(QQ,3),
		LRT_trecase = LRT_trec,LRT_cistrans = LRT_trec,
		mat_ETA = arma::zeros<arma::mat>(2,QQ),
		PVAL = LRT_trec;
	arma::vec tmp_ETA = arma::zeros<arma::vec>(QQ);
	if( hypotest && rcond_num > 0.0 && arma::any(upPARS.subvec(GI.at(3,0),GI.at(3,1)) == 1.0) ){
		// TReC eQTL significance
		upPARS_2 = upPARS;
		upPARS_2.subvec(GI.at(5,0),GI.at(5,1)).zeros();
		Rcpp_CSeQTL_gs_signif(final_TREC,final_LGX1,hap2,ASREC,LBC,
			final_PHASE*0,SNP,RHO_copy,XX,PARS_BPK,GI,I_np,upPARS_2,
			QQ,np,LRT_trec,tmp_ETA,false,max_iter,eps,gr_eps,conv_eps,show);
		mat_ETA.row(0) = tmp_ETA.t();
		PVAL.col(0) = LRT_trec.col(2);
		
		// TReCASE eQTL significance test: Assuming given cell type 
		//		is cis eQTL while all others are trans-eQTL
		if( useASREC == false ){
			LRT_trecase = LRT_trec;
			mat_ETA.row(1) = mat_ETA.row(0);
		} else {
			Rcpp_CSeQTL_gs_signif(final_TREC,final_LGX1,hap2,ASREC,
				LBC,final_PHASE,SNP,RHO_copy,XX,PARS_BPK,GI,I_np,
				upPARS,QQ,np,LRT_trecase,tmp_ETA,true,
				max_iter,eps,gr_eps,conv_eps,show);
			mat_ETA.row(1) = tmp_ETA.t();
		}
		PVAL.col(1) = LRT_trecase.col(2);
		
		// Cis/Trans test
		if( useASREC == false ){
			// report gene/SNP pair as trans eQTL if significant
			LRT_cistrans.col(0).fill(40.0);
			LRT_cistrans.col(1).fill(1);
			LRT_cistrans.col(2).fill(0);
			// aka report TReC estimates
		} else {
			Rcpp_CSeQTL_gs_cistrans(final_TREC,final_LGX1,hap2,ASREC,
				LBC,final_PHASE,SNP,RHO_copy,XX,PARS_trec,GI,I_np,upPARS,
				QQ,np,LRT_cistrans,max_iter,eps,gr_eps,conv_eps,show);
		}
		PVAL.col(2) = LRT_cistrans.col(2);
		
	} else {
		if( show ) Rcpp::Rcout << "No hypothesis testing ...\n";
	}
	
	if( swap_CT != 0 ){
		if( show ) Rcpp::Rcout << "Swap back cell types " << swap_CT + 1 << " and 1\n";
		PVAL.swap_rows(swap_CT,0); mat_ETA.swap_cols(swap_CT,0);
		MU.swap_cols(swap_CT,0);
		
		// upPARAM.swap_rows(swap_CT,0);
		// upPARS_2 = upPARS;
		// upPARS_2.at(0) = upPARS.at(GI.at(2,0) + swap_CT);
		// upPARS_2.at(GI.at(2,0) + swap_CT) = upPARS.at(0);
		// upPARS_2.at(GI.at(3,0)) = upPARS.at(GI.at(3,0) + swap_CT);
		// upPARS_2.at(GI.at(3,0) + swap_CT) = upPARS.at(GI.at(3,0));
		// upPARS_2.at(GI.at(5,0)) = upPARS.at(GI.at(5,0) + swap_CT);
		// upPARS_2.at(GI.at(5,0) + swap_CT) = upPARS.at(GI.at(5,0));
		// upPARS = upPARS_2;
		
		// eqtl_vars.swap_cols(swap_CT,0);
		LRT_trec.swap_rows(swap_CT,0);
		LRT_trecase.swap_rows(swap_CT,0);
		LRT_cistrans.swap_rows(swap_CT,0);
	}
	
	Rcpp::List OPT = Rcpp::List::create(
		Rcpp::Named("iBETA",Rcpp::NumericVector(iBETA.begin(),iBETA.end())),
		Rcpp::Named("iPHI",iPHI),Rcpp::Named("LL",LL),
		Rcpp::Named("PARS",Rcpp::NumericVector(PARS.begin(),PARS.end())),
		Rcpp::Named("GRAD",Rcpp::NumericVector(gk.begin(),gk.end())),
		Rcpp::Named("HESS",hess),Rcpp::Named("norm_GRAD",Rcpp_norm(gk)),
		Rcpp::Named("eqtl_vars",eqtl_vars),
		Rcpp::Named("norm_iHG",norm_iHG),Rcpp::Named("converge",converge),
		Rcpp::Named("upPARS",Rcpp::NumericVector(upPARS.begin(),upPARS.end())));
	
	Rcpp::List HYPO = Rcpp::List::create(
		Rcpp::Named("LRT_trec",LRT_trec),Rcpp::Named("LRT_trecase",LRT_trecase),
		Rcpp::Named("LRT_cistrans",LRT_cistrans),Rcpp::Named("PVAL",PVAL),
		Rcpp::Named("mat_ETA",mat_ETA),Rcpp::Named("swap_CT",swap_CT + 1));
	
	Rcpp::List SUBJ = Rcpp::List::create(
		Rcpp::Named("COOKSD",Rcpp::NumericVector(COOKSD.begin(),COOKSD.end())),
		Rcpp::Named("finTREC",Rcpp::NumericVector(final_TREC.begin(),final_TREC.end())),
		Rcpp::Named("log_OF",Rcpp::NumericVector(log_OF.begin(),log_OF.end())),
		Rcpp::Named("log_subj_mu",Rcpp::NumericVector(log_subj_mu.begin(),log_subj_mu.end())));
	
	if( show ) Rcpp::Rcout << "\n#####\n";
	
	return Rcpp::List::create(
		Rcpp::Named("OPT",OPT),Rcpp::Named("HYPO",HYPO),
		Rcpp::Named("SUBJ",SUBJ),Rcpp::Named("MU",MU),
		Rcpp::Named("ETA",mat_ETA));
	
}

// [[Rcpp::export]]
Rcpp::List Rcpp_CSeQTL_GS(const arma::mat& XX,const arma::mat& TREC_0,
	const arma::umat& SNP,const arma::mat& hap2,const arma::mat& ASREC,
	const arma::umat& PHASE_0,const arma::mat& RHO,const arma::umat& GS_index,
	const bool& trim = true,const bool& swapCT = true,const double& trim_thres = 10,
	const arma::uword& numAS = 5,const arma::uword& numASn = 5,
	const arma::uword& numAS_het = 5,const double& cistrans = 0.01,
	const double& mad_const = 1.4826,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-10,const double& gr_eps = 1e-2,const double& conv_eps = 5e-5,
	const bool& show = true,const bool& prompt = true,const int& ncores = 1){
	
	arma::uword NN = XX.n_rows, GG = TREC_0.n_rows, QQ = RHO.n_cols,
		PP = XX.n_cols, SS = SNP.n_rows, gg, num_snps, num_subj;
	arma::umat GI = Rcpp_calc_GI(PP,QQ);
	
	// Run TReC trimming
	arma::mat TREC = TREC_0;
	arma::vec CD = arma::zeros<arma::vec>(NN);
	if( trim ){
		if( prompt ){
			Rcpp::Rcout << "Running TReC Trimming ...\n";
			Rcpp::Rcout << "\tTrim threshold = " << trim_thres << "...\n";
		}
		arma::mat out_trim = arma::zeros<arma::mat>(NN,4);
		for(gg = 0; gg < GG; gg++){
			out_trim = Rcpp_CSeQTL_cooksD(TREC_0.row(gg).t(),RHO,XX,trim_thres,
				max_iter,eps,mad_const,ncores,gr_eps,conv_eps,false);
			CD = out_trim.col(0);
			TREC.row(gg) = out_trim.col(3).t();
			if( prompt && show ){
				Rcpp::Rcout << "\tProgress: " << gg + 1 << " out of " << GG << " done\n";
			}
		}
		if( prompt ){
			Rcpp::Rcout << "\tFinished Trimming ...\n";
		}
	}
	
	// Make LGX1 and LBC matrices
	if( prompt ) Rcpp::Rcout << "Making LGX1 and LBC matrices ...\n";
	arma::mat LGX1 = arma::zeros<arma::mat>(GG,NN),
		LBC = arma::zeros<arma::mat>(GG,NN);
	for(gg = 0; gg < GG; gg++){
		LGX1.row(gg) = arma::lgamma(TREC.row(gg) + 1.0);
		LBC.row(gg) = arma::lgamma(ASREC.row(gg) + 1.0) 
			- arma::lgamma(hap2.row(gg) + 1.0) 
			- arma::lgamma(ASREC.row(gg) - hap2.row(gg) + 1.0);
	}

	// Make PHASE matrix, based on PHASE_0, enforces threshold criteria from asSeq
	if( prompt ) Rcpp::Rcout << "Last QC check on ASREC ...\n";
	arma::umat PHASE = PHASE_0;
	arma::uvec tmp_PHASE = arma::zeros<arma::uvec>(NN);
	for(gg = 0; gg < GG; gg++){
		// Count number of subjects with at least numAS ASREC
		num_subj = arma::sum(ASREC.row(gg).t() >= numAS && PHASE_0.row(gg).t() == 1);
		if( num_subj < numASn ){
			PHASE.row(gg).zeros();
		} else { // Set PHASE = 0 for ASREC < numAS
			tmp_PHASE = PHASE.row(gg).t();
			if( arma::any( ASREC.row(gg).t() < numAS ) ){
				tmp_PHASE(arma::find(ASREC.row(gg).t() < numAS)).zeros();
				PHASE.row(gg) = tmp_PHASE.t();
			}
		}
	}
	arma::uword numASREC = arma::sum(arma::sum(PHASE,1) > 0);
	if( prompt ) Rcpp::Rcout << "Algorithm Specifications:"
		<< "\n   Number of threads = " << ncores
		<< "\n   Number of covariates (including intercept) = " << PP
		<< "\n   Number of genes = " << GG
		<< "\n   Number of snps across genes = " << SS
		<< "\n   Number of subjects = " << NN
		<< "\n   Number of cell types = " << QQ
		<< "\n   Number of genes w/ proper ASREC = " << numASREC
		<< "\n";
	
	arma::uword np = GI.at(5,1)+1, aa = 0, bb;
	arma::mat PARS_full = arma::zeros<arma::mat>(SS,np),
		ETA = arma::zeros<arma::mat>(SS,QQ),
		LRT_trec = ETA,LRT_trecase = LRT_trec,LRT_cistrans = LRT_trec,
		I_np = arma::eye<arma::mat>(np,np),
		MU_A = ETA,MU_B = ETA,RHO_copy = RHO;
	
	arma::vec TREC_gg = arma::zeros<arma::vec>(NN),
		LGX1_gg = TREC_gg, hap2_gg = TREC_gg, ASREC_gg = TREC_gg,
		LBC_gg = TREC_gg, iBETA_gg = arma::zeros<arma::vec>(PP),
		OO = arma::zeros<arma::vec>(NN),LL = arma::zeros<arma::vec>(SS);
	double iPHI_gg;
	arma::uvec PHASE_gg = arma::zeros<arma::uvec>(NN),
		no_PHASE = arma::zeros<arma::uvec>(NN);
	
	// Store NegBinom Regression results
	arma::mat iBETA = arma::zeros<arma::mat>(GG,PP);
	arma::vec iPHI = arma::zeros<arma::vec>(GG),
		tmp_eta = arma::ones<arma::vec>(QQ),tmp_alpha = tmp_eta;
	// Negative binomial regression
	Rcpp_CSeQTL_BETA_PHI(XX,TREC,iBETA,iPHI,show,ncores);
	
	// Loop over GENE/SNP pairs
	arma::uword count_snps = 0;
	for(gg = 0; gg < GG; gg++){
		if( prompt && show ){
			if( (gg + 1) % 5 == 0 ) Rcpp::Rcout << ".";
			if( (gg + 1) % 100 == 0 ) Rcpp::Rcout << gg + 1;
			if( ( (gg + 1) % 200 == 0 ) || ( (gg + 1) == GG ) ) Rcpp::Rcout << "(SNPs = " << count_snps << ")\n";
		}
		
		TREC_gg = TREC.row(gg).t(); LGX1_gg = LGX1.row(gg).t();
		hap2_gg = hap2.row(gg).t(); ASREC_gg = ASREC.row(gg).t();
		LBC_gg = LBC.row(gg).t(); PHASE_gg = PHASE.row(gg).t();
		iBETA_gg = iBETA.row(gg).t(); iPHI_gg = iPHI.at(gg);
		num_snps = GS_index.at(gg,1) - GS_index.at(gg,0) + 1;
		
		arma::umat SNP_gg = SNP.rows(GS_index.at(gg,0),GS_index.at(gg,1));
		bb = aa + num_snps - 1;
		arma::vec LL_gg = LL.subvec(aa,bb), upPARS = arma::ones<arma::vec>(np);
		arma::mat PARS_full_gg = PARS_full.rows(aa,bb),
			ETA_gg = ETA.rows(aa,bb),LRT_trec_gg = LRT_trec.rows(aa,bb),
			LRT_trecase_gg = LRT_trecase.rows(aa,bb),
			LRT_cistrans_gg = LRT_cistrans.rows(aa,bb),
			MU_A_gg = MU_A.rows(aa,bb),MU_B_gg = MU_B.rows(aa,bb);
		
		// Assess kappa estimation b/c ignoring SNP but accounting for RHO
		arma::uword swap_CT = 0;
		if( swapCT && QQ > 1 ){
			bool useTRIM = trim;
			upPARS = Rcpp_CSeQTL_swapCT(TREC_gg,LGX1_gg,RHO_copy,
				swap_CT,XX,GI,I_np,QQ,np,iBETA_gg,iPHI_gg,upPARS,useTRIM,
				max_iter,eps,gr_eps,conv_eps,false);
			
			// If we tried trimmed TREC but failed to converge, reset to untrimmed TREC
			if( trim == true && useTRIM == false ){
				if( show ) Rcpp::Rcout << "*** DON'T USE TRIMMED TREC, RESET to ORIG TREC!\n";
				
				// Set TREC and LGX1
				TREC_gg = TREC_0.row(gg).t();
				LGX1_gg = arma::lgamma(TREC_gg + 1.0);
				
				// Get iBETA,iPHI
				arma::vec iPARS = Rcpp_NB_reg_one(TREC_gg,XX,OO,max_iter,eps,false);
				iBETA_gg = iPARS.subvec(0,PP-1);
				iPHI_gg = std::exp(iPARS.at(PP));
				
				// Re-run swap
				upPARS.ones();
				RHO_copy = RHO;
				upPARS = Rcpp_CSeQTL_swapCT(TREC_gg,LGX1_gg,RHO_copy,
					swap_CT,XX,GI,I_np,QQ,np,iBETA_gg,iPHI_gg,upPARS,useTRIM,
					max_iter,eps,gr_eps,conv_eps,false);
			}
			
		}
		
		#ifdef _OPENMP
		# pragma omp parallel for schedule(dynamic) \
			num_threads(ncores) \
			shared(cistrans,max_iter,eps,gr_eps,conv_eps,prompt,show,ncores,\
			np,I_np,NN,GI,QQ,TREC_gg,LGX1_gg,hap2_gg,ASREC_gg,LBC_gg,\
			PHASE_gg,SNP_gg,XX,RHO_copy,iBETA_gg,iPHI_gg,upPARS,\
			PARS_full_gg,ETA_gg,LRT_trec_gg,LRT_trecase_gg,LRT_cistrans_gg,\
			LL_gg,MU_A_gg,MU_B_gg,num_snps,gg,numAS,numAS_het)
		#endif
		for(arma::uword ss = 0; ss < num_snps; ss++){
			if( prompt && show && ncores == 1 ){
				if( (ss + 1) % 5 == 0 ) Rcpp::Rcout << ".";
				if( (ss + 1) % 100 == 0 || (ss + 1) == num_snps ) Rcpp::Rcout << ss + 1;
				if( (ss + 1) % 300 == 0 || (ss + 1) == num_snps ) Rcpp::Rcout << "\n";
			}
			
			arma::vec PARS = arma::zeros<arma::vec>(np),
				tmp_ETA = arma::zeros<arma::vec>(QQ);
			arma::uvec PHASE_snp = arma::zeros<arma::uvec>(NN);
			arma::mat tmp_LRT_trec = arma::zeros<arma::mat>(QQ,3),
				tmp_LRT_trecase = tmp_LRT_trec, tmp_LRT_cistrans = tmp_LRT_trec,
				MU_full = arma::zeros<arma::mat>(2,QQ);
			double tmp_LL;
			
			// Determine if haplotype counts should be used
			if( arma::any(PHASE_gg == 1) ){
				// Count num subjects with PHASE == 1 & ASREC >= numAS & (SNP == 1 | SNP == 2)
				arma::uvec tmp_SNP = SNP_gg.row(ss).t();
				arma::uword hap_cnt = arma::sum( PHASE_gg == 1 
					&& ASREC_gg >= numAS 
					&& (tmp_SNP == 1 || tmp_SNP == 2) );
				if( hap_cnt >= numAS_het ) PHASE_snp = PHASE_gg;
			}
			
			// Run analysis
			bool show2 = show && ncores == 1;
			Rcpp_CSeQTL_genesnp(TREC_gg,LGX1_gg,hap2_gg,ASREC_gg,LBC_gg,
				PHASE_snp,SNP_gg.row(ss).t(),RHO,XX,GI,I_np,QQ,np,iBETA_gg,
				iPHI_gg,upPARS,tmp_LL,PARS,tmp_ETA,tmp_LRT_trec,tmp_LRT_trecase,
				tmp_LRT_cistrans,MU_full,cistrans,max_iter,eps,gr_eps,
				conv_eps,show2);
			
			// Calculate log likelihood for TReC model
			LL_gg.at(ss) = tmp_LL;
			
			// Store results
			PARS_full_gg.row(ss) 		= PARS.t();
			ETA_gg.row(ss) 					= tmp_ETA.t();
			LRT_trec_gg.row(ss) 		= tmp_LRT_trec.col(0).t();
			LRT_trecase_gg.row(ss) 	= tmp_LRT_trecase.col(0).t();
			LRT_cistrans_gg.row(ss) = tmp_LRT_cistrans.col(0).t();
			
			MU_A_gg.row(ss) = MU_full.row(0);
			MU_B_gg.row(ss) = MU_full.row(1);
		}
		
		// Reorder ETA, LRT matrices
		if( swap_CT != 0 ){
			ETA_gg.swap_cols(swap_CT,0);
			LRT_trec_gg.swap_cols(swap_CT,0);
			LRT_trecase_gg.swap_cols(swap_CT,0);
			LRT_cistrans_gg.swap_cols(swap_CT,0);
			RHO_copy.swap_cols(swap_CT,0); // swap back the original matrix
			MU_A_gg.swap_cols(swap_CT,0);
			MU_B_gg.swap_cols(swap_CT,0);
		}
		
		// Save SNP results to final output
		PARS_full.rows(aa,bb) 		= PARS_full_gg;
		ETA.rows(aa,bb) 					= ETA_gg;
		LRT_trec.rows(aa,bb) 			= LRT_trec_gg;
		LRT_trecase.rows(aa,bb) 	= LRT_trecase_gg;
		LRT_cistrans.rows(aa,bb) 	= LRT_cistrans_gg;
		LL.subvec(aa,bb) = LL_gg;
		MU_A.rows(aa,bb) = MU_A_gg;
		MU_B.rows(aa,bb) = MU_B_gg;
		aa = bb + 1;
		count_snps += num_snps;
	}
	
	// Proportion subjects trimmed
	double propTRIM = ( 1.0 * arma::sum(TREC_0.row(0).t() != TREC.row(0).t()) ) / NN;
	
	// Output
	return Rcpp::List::create(Rcpp::Named("iBETA",iBETA),
		Rcpp::Named("iPHI",Rcpp::NumericVector(iPHI.begin(),iPHI.end())),
		Rcpp::Named("CD",Rcpp::NumericVector(CD.begin(),CD.end())),
		Rcpp::Named("PARS",PARS_full),Rcpp::Named("MU_A",MU_A),
		Rcpp::Named("MU_B",MU_B),Rcpp::Named("ETA",ETA),
		Rcpp::Named("LRT_trec",LRT_trec),Rcpp::Named("LRT_trecase",LRT_trecase),
		Rcpp::Named("LRT_cistrans",LRT_cistrans),Rcpp::Named("propTRIM",propTRIM),
		Rcpp::Named("LL",Rcpp::NumericVector(LL.begin(),LL.end())));
	
}






