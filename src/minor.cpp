#include "../inst/include/minor.h"

// Minor Functions

double Rcpp_norm(const arma::vec& a){
	return arma::norm(a);
}

// [[Rcpp::export]]
arma::umat Rcpp_calc_GI(const arma::uword& PP,const arma::uword& QQ){
	// Setting GI
	arma::uword jj, index = 0;
	arma::umat GI = arma::zeros<arma::umat>(6,2);
	for(jj = 0; jj < 6; jj++){
		if(jj == 0){ 						// BETA
			GI.at(jj,0) = index;
			index += PP - 1;
			GI.at(jj,1) = index;
		} else if(jj == 1){ 		// PHI
			GI.at(jj,0) = index;
			GI.at(jj,1) = index;
		} else if(jj == 2){ 		// KAPPA2Q
			if(QQ > 1){
				GI.at(jj,0) = index;
				index += QQ - 2;
				GI.at(jj,1) = index;
			} else {
				GI.row(jj) = GI.row(jj-1);
			}
		} else if(jj == 3){ 		// ETA
			GI.at(jj,0) = index;
			index += QQ - 1;
			GI.at(jj,1) = index;
		} else if(jj == 4){ 		// PSI
			GI.at(jj,0) = index;
			GI.at(jj,1) = index;
		} else { 								// ALPHA
			GI.at(jj,0) = index;
			index += QQ - 1;
			GI.at(jj,1) = index;
		}
		
		if(jj == 2 && QQ == 1) continue;
		index++;
	}
	
	return GI;
}

// [[Rcpp::export]]
arma::mat Rcpp_CSeQTL_MU(const arma::umat& GI,
	const arma::vec& PARS){
	
	arma::vec ETA = arma::exp(PARS.subvec(GI.at(3,0),GI.at(3,1)));
	arma::uword QQ = ETA.n_elem, ct;
	arma::vec KAPPA = arma::ones<arma::vec>(QQ);
	if(QQ > 1) KAPPA.subvec(1,QQ-1) = arma::exp(PARS.subvec(GI.at(2,0),GI.at(2,1)));
	arma::mat MU = arma::zeros<arma::mat>(2,QQ);
	
	for(ct = 0; ct < QQ; ct++){
		MU.at(0,ct) = std::exp(PARS.at(0)) / 2.0 * KAPPA.at(ct);
		MU.at(1,ct) = MU.at(0,ct) * ETA.at(ct);
	}
	
	return MU;
}

double Rcpp_calc_MAF(const arma::vec& SNP,const bool& phasing = true,
	const bool& show = true){
	
	// If phasing == true, genotypes are coded 0/1/2/3 for AA/AB/BA/BB
	// Otherwise, genotypes are coded 0/1/2 for AA/AB/BB
	
	arma::uword num_counts, gg;
	double MAF;
	
	if( phasing ){
		num_counts = 4;
	} else {
		num_counts = 3;
	}
	
	arma::vec counts = arma::zeros<arma::vec>(num_counts);
	for(gg = 0; gg < num_counts; gg++){
		arma::uvec idx = arma::find(SNP == gg);
		counts.at(gg) = idx.n_elem;
	}
	
	if( show ) counts.t().print("counts = ");
	
	if( phasing ){
		MAF = counts.at(1) + counts.at(2) + 2*counts.at(3);
	} else {
		MAF = counts.at(1) + 2*counts.at(2);
	}
	
	MAF /= 2 * arma::sum(counts);
	return MAF;
}

// [[Rcpp::export]]
arma::vec Rcpp_calc_MAF_all(const arma::mat& SNP,
	const bool& phasing = true,const bool& show = true,
	const int& ncores = 1){
	
	arma::uword num_snps = SNP.n_rows;
	arma::vec MAF = arma::zeros<arma::vec>(num_snps);
	
	#ifdef _OPENMP
	# pragma omp parallel for schedule(dynamic) \
		num_threads(ncores) \
		shared(SNP,phasing,num_snps,ncores,show)
	#endif
	for(arma::uword ss = 0; ss < num_snps; ss++){
		
		if( show && ncores == 1 ){
			if( (ss+1) % 10 == 0 )
				Rcpp::Rcout << ".";
			if( (ss+1) % 200 == 0 || (ss+1) == num_snps )
				Rcpp::Rcout << ss+1 << " out of " << num_snps << "\n";
		}
		
		MAF.at(ss) = Rcpp_calc_MAF(SNP.row(ss).t(),phasing,false);
	}
	return MAF;
}

