#include "../inst/include/pack_GRAD.h"
#include "../inst/include/pack_HESS.h"

// --------------------
// Hessian related functions

arma::mat Rcpp_CSeQTL_HESS(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::vec& BETA,
	const double& PHI,const double& PSI,const arma::vec& KAPPA2Q,
	const arma::vec& ETA_TREC,const arma::vec& ALPHA,const arma::umat& GI,
	const arma::vec& upPARS,const arma::mat& I_np,const double& shift,
	const bool& show = true){
	
	arma::uword np = GI.at(5,1) + 1, ii, QQ = RHO.n_cols;
	arma::mat HESS = arma::zeros<arma::mat>(np,np);
	arma::vec PARS = arma::zeros<arma::vec>(np),PARS_2 = PARS,
		tmp_vec = PARS, GRAD_xx = PARS;

	if( show ) Rcpp::Rcout << "upPARS = " << upPARS.t();
	GRAD_xx = Rcpp_CSeQTL_GRAD(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,
		BETA,PHI,PSI,KAPPA2Q,ETA_TREC,ALPHA,GI,upPARS);
	// bool check_phase = arma::any(PHASE == 1);
	
	PARS.subvec(GI.at(0,0),GI.at(0,1)) = BETA;
	PARS.at(GI.at(1,0)) = std::log(PHI);
	if(QQ > 1) PARS.subvec(GI.at(2,0),GI.at(2,1)) = arma::log(KAPPA2Q);
	PARS.subvec(GI.at(3,0),GI.at(3,1)) = arma::log(ETA_TREC);
	PARS.at(GI.at(4,0)) = std::log(PSI);
	PARS.subvec(GI.at(5,0),GI.at(5,1)) = arma::log(ALPHA);
	
	for(ii = 0; ii < np; ii++){
		if( upPARS.at(ii) != 0.0 ){
			PARS_2 = PARS + shift * I_np.col(ii) % upPARS;
			tmp_vec = Rcpp_CSeQTL_calc_GRAD(TREC,hap2,ASREC,PHASE,SNP,
				RHO,XX,GI,upPARS,PARS_2) - GRAD_xx;
			tmp_vec /= shift;
			if( show ) tmp_vec.t().print("hess_vec = ");
			// HESS.col(ii) = tmp_vec;
			HESS(arma::span(ii,np - 1),ii) = tmp_vec.subvec(ii,np - 1);
			if(ii < np - 1){
				HESS(ii,arma::span(ii+1,np - 1)) = tmp_vec.subvec(ii+1,np - 1).t();
			}
		}
	}
	
	return HESS;
}

arma::mat Rcpp_CSeQTL_calc_HESS(const arma::vec& TREC,const arma::vec& hap2,
	const arma::vec& ASREC,const arma::uvec& PHASE,const arma::uvec& SNP,
	const arma::mat& RHO,const arma::mat& XX,const arma::umat& GI,const arma::vec& PARS,
	const arma::vec& upPARS,const arma::mat& I_np,const double& shift){
	// shift = 1e-5
	
	return Rcpp_CSeQTL_HESS(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,
		PARS.subvec(GI.at(0,0),GI.at(0,1)),
		std::exp(PARS.at(GI.at(1,0))),
		std::exp(PARS.at(GI.at(4,0))),
		arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1))),
		arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1))),
		arma::exp(PARS.subvec(GI.at(5,0),GI.at(5,1))),
		GI,upPARS,I_np,shift,false);
}

void Rcpp_CSeQTL_hessBR(const arma::mat& hess,const arma::umat& GI,
	const arma::vec& upPARS,bool& rcond_nz,const bool& show){
	
	// Get the variances of the bottom right of the hessian
	arma::uword PP = GI.at(0,1) + 1, np = GI.at(5,1) + 1,
		nKEPA = np - PP - 1, QQ = GI.at(5,1) - GI.at(5,0) + 1;
	
	// BETA + PHI parts
	arma::mat A = hess.submat(GI.at(0,0),GI.at(0,0),GI.at(1,1),GI.at(1,1)),
		inv_A = A;
	inv_A.zeros();
	arma::uvec nz_A = arma::find( upPARS.subvec(GI.at(0,0),GI.at(1,1)) == 1.0 );
	if( arma::rcond(A.submat(nz_A,nz_A)) <= 0.0 ){
		rcond_nz = false;
		return;
	}
	
	// BETA + PHI vs KAPPA + ETA + PSI + ALPHA
	arma::mat B = arma::zeros<arma::mat>(PP + 1,nKEPA);
	
	// KAPPA + ETA + PSI + ALPHA
	arma::mat D = arma::zeros<arma::mat>(nKEPA,nKEPA),iKEPA = D;
	if( QQ == 1 ){
		B = hess.submat(GI.at(0,0),GI.at(3,0),GI.at(1,1),GI.at(5,1));
		D = hess.submat(GI.at(3,0),GI.at(3,0),GI.at(5,1),GI.at(5,1));
	} else {
		B = hess.submat(GI.at(0,0),GI.at(2,0),GI.at(1,1),GI.at(5,1));
		D = hess.submat(GI.at(2,0),GI.at(2,0),GI.at(5,1),GI.at(5,1));
	}
	
	if( show ){
		Rcpp::Rcout << "===== Rcpp_CSeQTL_hessBR() =====\n";
		A.print("A = ");
		B.print("B = ");
		D.print("D = ");
	}
	
	inv_A.submat(nz_A,nz_A) = arma::inv(A.submat(nz_A,nz_A));
	
	iKEPA = D - B.t() * inv_A * B;
	if( show ) iKEPA.print("D - t(B) * inv(A) * B = ");
	arma::uvec nz = arma::find( iKEPA.diag() != 0.0 );
	rcond_nz = arma::rcond(iKEPA.submat(nz,nz)) > 0.0;
	
	if( rcond_nz ){ // If above rcond checks pass, check the full hess matrix
		arma::uvec nz_full = arma::find( hess.diag() != 0.0 );
		double rcond_full = arma::rcond(hess.submat(nz_full,nz_full));
		rcond_nz = rcond_full > 0.0;
		if( show ) Rcpp::Rcout << "Rcond hess_nz = " << rcond_full << "\n";
	}
	
}

