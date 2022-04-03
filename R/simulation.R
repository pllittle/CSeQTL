
# ----------
# CSeQTL Simulation Functions
# ----------

calc_POWER_2 = function(PVAL,alpha = 0.05){
	PVAL = PVAL[!is.na(PVAL)]
	mean(PVAL <= alpha)
}

#' @title gen_true_RHO
#' @description Simulates cell type proportions under three scenarios.
#' @param wRHO An integer taking values 1, 2, or 3 for one of the three
#'	scenarios. For \code{wRHO} equal to 2 or 3, \code{QQ} needs to be 
#'	set to 3 cell types.
#' @param NN Positive integer for sample size.
#' @param QQ Positive integer for number of cell types.
#' @param RHO Default to \code{NULL} leads to simulating cell 
#'	type proportions. If a matrix of proportions is supplied, 
#'	this function will append row/column names.
#' @export
gen_true_RHO = function(wRHO = 1,NN,QQ,RHO = NULL){
	if(FALSE){
		wRHO = 3; NN = 3e2; QQ = 3; RHO = NULL
	}
	
	if( is.null(RHO) ){
		if( wRHO == 1 ){
			true_RHO = matrix(exp(4*runif(NN*QQ,-1,1)),NN,QQ)
		} else if( wRHO == 2 ){
			if( QQ != 3 ) stop("For wRHO = 2, set QQ = 3")
			true_RHO = matrix(NA,NN,QQ)
			true_RHO[,1] = rbeta(n = NN,shape1 = 10,shape2 = 24)
			true_RHO[,2] = 0.85 + -0.76 * true_RHO[,1] + 
				-0.03 * true_RHO[,1]^2 + rnorm(NN,0,0.02)
			true_RHO[,2] = abs(true_RHO[,2])
			true_RHO[,3] = 1 - rowSums(true_RHO[,1:2])
			neg_idx = which(true_RHO[,3] <= 0)
			if( length(neg_idx) > 0 ) true_RHO[neg_idx,3] = 0
			true_RHO = true_RHO / rowSums(true_RHO)
		} else if( wRHO == 3 ){
			RHO = gen_true_RHO(wRHO = 2,NN = NN,QQ = QQ)
			alpha = 0.01
			for(CT in seq(QQ)){
				one_rho = RHO[,CT]
				idx_high = which(one_rho > quantile(one_rho,1-alpha))
				idx_low = which(one_rho < quantile(one_rho,alpha))
				one_rho[idx_high] = runif(length(idx_high),0.7,0.9)
				one_rho[idx_low] = runif(length(idx_low),0.0,0.1)
				RHO[,CT] = one_rho
			}
			RHO = RHO / rowSums(RHO)
			true_RHO = RHO
		}
		
	} else {
		NN = nrow(RHO)
		QQ = ncol(RHO)
		true_RHO = RHO
	}
	
	true_RHO = true_RHO / rowSums(true_RHO)
	if( any(true_RHO == 0) ){
		true_RHO = true_RHO + 1e-5
	}
	true_RHO = true_RHO / rowSums(true_RHO)
	
	if( is.null(colnames(true_RHO)) ){
		if( QQ > 1 ){
			colnames(true_RHO) = paste0("CT",seq(QQ))
		} else {
			colnames(true_RHO) = "MARG"
		}
	}
	rownames(true_RHO) = paste0("Subj",seq(NN))
	
	if(FALSE){
		boxplot(true_RHO)
		
	}
	
	true_RHO
}
gen_gene_TREC_ASREC = function(XX,RHO,gene_BETA,gene_PHI,
	gene_SNP_eQTL,gene_KAPPA,gene_ETA,gene_ALPHA,gene_PSI,prob_phased){
	
	if(FALSE){
		XX = XX; RHO = true_RHO; gene_BETA = true_BETA; gene_PHI = true_PHI
		gene_SNP_eQTL = true_SNP; gene_KAPPA = c(1,0.1,1); gene_ETA = c(1,1,1)
		gene_ALPHA = true_ALPHA; gene_PSI = true_PSI; prob_phased = 0.05
	}
	
	NN = nrow(XX); dat = smart_df(XX); dat$Z = gene_SNP_eQTL
	dat$mu = NA; dat$total = NA; dat$total_phased = NA
	dat$hapA = 0; dat$hap2 = 0; dat$LBC = 1; dat$phased = NA;
	# hapA = 1st haplotype, hap2 = 2nd haplotype
	dat$pp = NA; dat$log_mu = NA;
	KAPPA_ETA = gene_KAPPA * gene_ETA
	RHO_KAPPA = c(RHO %*% gene_KAPPA) # \sum_(q=1)^Q rho_iq * kappa_gq
	vec_log_mu_AA = c(as.matrix(XX) %*% gene_BETA) + log(RHO_KAPPA)
	vec_XI_TREC = c(RHO %*% KAPPA_ETA) / RHO_KAPPA
	vec_XI_ASREC = c(RHO %*% (gene_ALPHA * KAPPA_ETA)) / RHO_KAPPA
	
	# double check calculations
	# ii = 1; RHO[ii,]; 
	# RHO_KAPPA[ii]; sum(RHO[ii,] * gene_KAPPA)
	# vec_XI_TREC[ii]; sum(RHO[ii,] * gene_KAPPA * gene_ETA) / sum(RHO[ii,] * gene_KAPPA)
	
	genotypes = seq(0,3)
	for(zz in genotypes){
		idx = which(dat$Z == zz)
		if( zz == 0 ){
			dat$log_mu[idx] = vec_log_mu_AA[idx]
			dat$pp[idx] = 0.5
		} else if( zz %in% c(1,2) ){
			dat$log_mu[idx] = vec_log_mu_AA[idx] + log( (1 + vec_XI_TREC[idx]) / 2 )
			dat$pp[idx] = vec_XI_ASREC[idx] / (1 + vec_XI_ASREC[idx])
		} else {
			dat$log_mu[idx] = vec_log_mu_AA[idx] + log( vec_XI_TREC[idx] )
			dat$pp[idx] = 0.5
		}
	}
	
	dat$mu = exp(dat$log_mu)
	if( gene_PHI > 0.0 ){
		dat$total = rnbinom(n = NN,mu = dat$mu,size = 1/gene_PHI)
	} else {
		dat$total = rpois(n = NN,lambda = dat$mu)
	}
	dat$total_phased = rbinom(n = NN,size = dat$total,prob = prob_phased)
	dat$phased = ifelse(dat$total_phased > 0,1,0)
	for(Z in c(0,1,2,3)){
		# Z = 0
		tmp_index = which(dat$phased == 1 & dat$Z == Z)
		len_index = length(tmp_index)
		if(len_index > 0){
			tmp_ASREC = dat$total_phased[tmp_index]
			tmp_prob = dat$pp[tmp_index]
			if(Z == 1){ 				# AB
				if( gene_PSI > 0.0 ){
					dat$hap2[tmp_index] = rbetabinom(n = len_index,prob = tmp_prob,size = tmp_ASREC,theta = 1/gene_PSI)
				} else {
					dat$hap2[tmp_index] = rbinom(n = len_index,prob = tmp_prob,size = tmp_ASREC)
				}
				dat$hapA[tmp_index] = tmp_ASREC - dat$hap2[tmp_index]
			} else if(Z == 2){ 	# BA
				if( gene_PSI > 0.0 ){
					dat$hapA[tmp_index] = rbetabinom(n = len_index,prob = tmp_prob,size = tmp_ASREC,theta = 1/gene_PSI)
				} else {
					dat$hapA[tmp_index] = rbinom(n = len_index,prob = tmp_prob,size = tmp_ASREC)
				}
				dat$hap2[tmp_index] = tmp_ASREC - dat$hapA[tmp_index]
			} else { 						# AA, BB
				if( gene_PSI > 0.0 ){
					dat$hapA[tmp_index] = rbetabinom(n = len_index,prob = tmp_prob,size = tmp_ASREC,theta = 1/gene_PSI)
				} else {
					dat$hapA[tmp_index] = rbinom(n = len_index,prob = tmp_prob,size = tmp_ASREC)
				}
				dat$hap2[tmp_index] = tmp_ASREC - dat$hapA[tmp_index]
			}
		}
	}

	dat$LBC = lchoose(dat$total_phased,dat$hap2)
	dat$LGX1 = lgamma(dat$total + 1)
	dat = dat[,which(!(names(dat) %in% c("Z")))]
	dat[,c("total","LGX1","total_phased","hap2",
		"LBC","phased","log_mu","mu","pp")]
}

#' @title CSeQTL_dataGen
#' @description Simulates a gene/SNP pair with baseline covariates \code{XX},
#'	cell type compositions \code{true_RHO}, phased SNP genotypes \code{true_SNP},
#'	and total (TReC) and allele-specific read counts (ASReC) contained in \code{dat}.
#' 
#' @param NN Positive integer for sample size.
#' @param MAF Positive numeric value between 0 and 1 for the minor 
#'	allele frequency to simulate phased SNP genotypes assuming Hardy-Weinberg.
#' @param true_BETA0 A positive numeric value denoting the reference cell type
#'	and reference base's expression multiplied by two and log transformed.
#'	For example, if the TReC for reference base and cell type is 500, then
#'	\code{true_BETA0 = log{2 * 500}}.
#' @param true_KAPPA A numeric vector denoting the baseline fold change in TReC
#'	between a cell type and reference. By definition, the first element is 1.
#' @param true_ETA A numeric vector where each element denotes the fold change 
#'	in TReC between the non-reference and reference base in a cell type.
#' @param true_PHI A non-negative numeric value denoting the over-dispersion term
#'	associated with TReC. If \code{true_PHI > 0}, TReC is simulated with the 
#'	negative binomial. If \code{true_PHI = 0}, TReC is simulated with the poisson.
#' @param true_PSI A non-negative numeric value denoting the over-dispersion term
#'	associated with ASReC. If \code{true_PSI > 0}, ASReC is simulated with the
#'	beta-binomial, otherwise it is simulated with the binomial distribution.
#' @param prob_phased A positive numeric value denoting the simulated proportion of 
#'	simulated TReC that are ASReC.
#' @param true_ALPHA By default, it is set to \code{NULL} setting each cell 
#'	type with an eQTL to be cis-eQTL. Otherwise, a positive numeric vector
#'	of fold changes between TReC eQTL effect sizes and ASReC eQTL effect sizes.
#' @param batch A numeric value set to 1 by default to allow underlying batch effects. 
#'	Set to zero to eliminate batch effects.
#' @param RHO A numeric matrix of cell type proportions where each row sums to one.
#'	If set to \code{NULL}, a matrix of cell type proportions will be simulated.
#' @param cnfSNP A boolean value where \code{TRUE} re-arranges simulated SNPs to
#'	correlate with baseline bulk expression. When fitting the marginal model 
#'	(not accounting for cell type proportions) and in the presence of cell 
#'	type-specific differentiated expression, a marginal eQTL may be incorrectly inferred.
#' @param show A boolean value to display verbose output and plot intermediate 
#'	simulated results.
#'
#' @export
CSeQTL_dataGen = function(NN,MAF,true_BETA0 = log(1e3),true_KAPPA,true_ETA,
	true_PHI = 0.1,true_PSI = 0.05,prob_phased = 0.05,true_ALPHA = NULL,
	batch = 1,RHO = NULL,cnfSNP = FALSE,show = TRUE){
	
	if(FALSE){
		NN = 3e2; MAF = 0.4; true_BETA0 = log(1e3); true_KAPPA = c(1,2,2)
		true_ETA = rep(1,3); true_PHI = 0.1; RHO = NULL; cnfSNP = TRUE; show = !TRUE
	 
	}
	
	if( !(batch %in% c(0,1)) ) stop("batch should be 0 or 1")
	
	# Simulate SNP eQTL
	if( show ) cat(sprintf("%s: Simulate phased SNPs ...\n",date()))
	geno_probs 	= c((1-MAF)^2,MAF*(1-MAF),MAF*(1-MAF),MAF^2)
	genotypes 	= seq(0,3) # 0/1/2/3 <=> AA/AB/BA/BB
	true_SNP 		= sample(genotypes,NN,replace = TRUE,prob = geno_probs)
	
	if( show ){
		# cat("Simulate SNP eQTL ...\n")
		print(smart_table(true_SNP))
	}
	
	# Simulate RHO
	if( show ) cat(sprintf("%s: Simulate cell type proportions ...\n",date()))
	QQ = length(true_ETA)
	true_RHO = gen_true_RHO(NN = NN,QQ = QQ,RHO = RHO) 
	# generate cell type proportion matrix, each cell type has sufficient variability
	if(show && QQ > 1) boxplot(true_RHO,xlab = "Cell Type",ylab = "Proportion")
	
	if( cnfSNP && all(true_ETA == 1) ){
		if( show ) cat(sprintf("%s: Re-arrange simulated SNPs to induce confounding ...\n",date()))
		# true_RHO = sim$true_RHO; true_KAPPA = c(1,2,3); NN = 3e2; true_SNP = sim$true_SNP
		OO = log(c(true_RHO %*% true_KAPPA)); OO[1:10]
		true_SNP2 = rep(NA,NN); agg_rank = 0
		for(mc in c(0,1,2)){
			# mc = 0
			if( mc == 0 ){
				num_geno = length(which(true_SNP == 0))
			} else if( mc == 1 ){
				num_geno = length(which(true_SNP %in% c(1,2)))
			} else {
				num_geno = length(which(true_SNP == 3))
			}
			agg_rank = agg_rank + num_geno
			idx = which(is.na(true_SNP2) & rank(OO) <= agg_rank)
			length(idx)
			if( mc == 0 ){
				true_SNP2[idx] = mc
			} else if( mc == 1 ){
				tmp_genos = true_SNP[true_SNP %in% c(1,2)]
				true_SNP2[idx] = sample(tmp_genos,length(idx))
			} else {
				true_SNP2[idx] = 3
			}
		}
		if( rbinom(1,1,0.5) == 1 ) true_SNP2 = 3 - true_SNP2
		true_SNP = true_SNP2
	}
	
	# Simulate covariates
	if( show ) cat(sprintf("%s: Simulate covariates ...\n",date()))
	XX = matrix(NA,NN,5)
	XX[,1] = 1; # intercept
	log_lib_size = rgamma(n = NN,shape = 600,rate = 100) # log(library size)
		XX[,2] = scale(log_lib_size)[,1]
		# XX[,2] = log_lib_size
	XX[,3] = rbinom(NN,1,0.5) 			# categorical variable
	XX[,4] = runif(NN,-1,1)   			# uniform dist. variable
		XX[,4] = scale(XX[,4])[,1]
	# XX[,5] = scale(rnorm(NN))[,1]		# norm dist. variable
	XX[,5] = rnorm(NN)

	# Set parameters
	vGAMMA = c(-0.5,0.2,-0.2) * batch
	true_BETA = c(true_BETA0,0.5,vGAMMA);  # regression covariates
	# KAPPA:  log fold change on A allele between q-th and 1st cell type
	# ETA:    log fold change between B and A allele for q-th cell type, eQTL effect size
	if( is.null(true_ALPHA) ){
		true_ALPHA = rep(1,QQ)
	}
	eta_one_idx = which(true_ETA == 1)
	true_ALPHA[eta_one_idx] = 1
	
	if( show ){
		cat(sprintf("%s: True parameter values ...\n",date()))
		cat(sprintf("   BETA = (%s)\n",paste(round(true_BETA,3),collapse=",")))
		cat(sprintf("   PHI = %s\n",paste(true_PHI,collapse=",")))
		cat(sprintf("   KAPPA = (%s)\n",paste(true_KAPPA,collapse=",")))
		cat(sprintf("   ETA = (%s)\n",paste(true_ETA,collapse=",")))
		cat(sprintf("   PSI = (%s)\n",paste(true_PSI,collapse=",")))
		cat(sprintf("   ALPHA = (%s)\n",paste(true_ALPHA,collapse=",")))
	}
	
	# Simulate outcomes
	if( show ) cat(sprintf("%s: Simulate TReC and ASReC ...\n",date()))
	XX2 = XX; cont_vars = c(2,4,5)
	XX2[,cont_vars] = apply(XX2[,cont_vars],2,function(xx) scale(xx)[,1])
	XX = XX2
	dat = gen_gene_TREC_ASREC(XX = XX,
		RHO = true_RHO,gene_BETA = true_BETA,gene_PHI = true_PHI,
		gene_SNP_eQTL = true_SNP,gene_KAPPA = true_KAPPA,gene_ETA = true_ETA,
		gene_ALPHA = true_ALPHA,gene_PSI = true_PSI,prob_phased = prob_phased)
	dat$phase0 = rep(0,NN)
	dat$log_lib_size = log_lib_size
	PP = ncol(XX); GI = Rcpp_calc_GI(PP = PP,QQ = QQ); np = GI[6,2]+1
	I_np = diag(np)
	
	# Rcpp code to initiate negative binom params
	vz = rep(0,NN)
	iPARS = Rcpp_NB_reg_one(YY = dat$total,XX = XX,
		OO = vz,max_iter = 4e3,eps = 1e-5,show = FALSE)
	iBETA = iPARS[seq(PP)]; iPHI = exp(iPARS[PP+1])
	
	dat$geno_col = rep("",NN)
	dat$geno_col[true_SNP == 0] = rgb(1,0,0,0.75)
	dat$geno_col[true_SNP %in% c(1,2)] = rgb(0,1,0,0.75)
	dat$geno_col[true_SNP == 3] = rgb(0,0,1,0.75)
	
	true_PARS = c(true_BETA,log(true_PHI))
	vname = c(paste0("beta",seq(PP)),"logPhi")
	if( QQ > 1 ){
		true_PARS = c(true_PARS,log(true_KAPPA[-1]))
		vname = c(vname,paste0("logK",seq(2,QQ)))
	}
	true_PARS = c(true_PARS,log(true_ETA),log(true_PSI),log(true_ALPHA))
	vname = c(vname,paste0("logE",seq(QQ)),"logPsi",paste0("logA",seq(QQ)))
	names(true_PARS) = vname
	
	# Calculate TReC offset term
	OF = log(true_RHO %*% true_KAPPA)
	for(ii in seq(NN)){
		tmp_xi = sum(true_RHO[ii,] * true_ETA * true_KAPPA) /
			sum(true_RHO[ii,] * true_KAPPA)
		
		if( true_SNP[ii] %in% c(1,2) ){
			OF[ii] = OF[ii] + log( ( 1 + tmp_xi ) / 2 )
		} else if( true_SNP[ii] == 3 ){
			OF[ii] = OF[ii] + log( tmp_xi )
		}
	}
	
	colnames(XX) = paste0("X",seq(ncol(XX)))
	list(true_PARS = true_PARS,true_SNP = true_SNP,dat = dat,XX = XX,
		true_RHO = true_RHO,QQ = QQ,GI = GI,np = np,I_np = I_np,
		iBETA = iBETA,iPHI = iPHI,MU = Rcpp_CSeQTL_MU(GI = GI,PARS = true_PARS),
		true_OF = OF,vname = vname)
}

#' @title CSeQTL_oneExtremeSim
#' @description Performs a simulation with multiple replicates
#'	for a pre-specified set of arguments.
#' 
#' @param wRHO Takes integer values 1, 2, or 3 to simulate three 
#'	scenarios of cell type proportions.
#' @param noiseRHO A positive numeric value to purposely distort simulated 
#'	cell type compositions.
#' @param RR A positive integer for number of replicates to generate and analyze.
#' @param vec_MARG A boolean vector for marginal and/or cell type-specific
#'	analyses to be run. By default, both sets of analyses are run.
#' @param vec_TRIM A boolean vector for whether or not analyses with 
#'	trimmed outcomes are included. By default, both sets of analyses are run.
#' @param vec_PERM A boolean vector for whether or not permuted SNP analyses
#'	are included. By default, both sets of analyses are run.
#' @param thres_TRIM A positive numeric value to perform subject outcome trimming.
#'	Subjects with standardized Cooks' Distances greater than the threshold are trimmed.
#' @param ncores A positive integer specifying the number of threads available
#'	to decrease computational runtime.
#' @param sim_fn Character value specifying the full path and filename
#'	to store intermediate simulation replicates should errors arise.
#' @inheritParams CSeQTL_dataGen
#'
#' @export
CSeQTL_oneExtremeSim = function(NN,MAF,true_BETA0,true_KAPPA,true_ETA,
	true_PHI = 0.1,wRHO,noiseRHO = 0,RR = 5e1,vec_MARG = c(TRUE,FALSE),
	vec_TRIM = c(TRUE,FALSE),vec_PERM = c(TRUE,FALSE),thres_TRIM = 10,
	ncores = 1,sim_fn){
	
	if(FALSE){
		NN = 3e2; MAF = 0.2; true_BETA0 = 5.3; 
		true_PHI = 0.2; 
		RR = 1e1
		
		kk = 6; ee = 7; ww = 3; nw = 2
		sim_fn = file.path(simext_dir,sprintf("sim_kk%s_ee%s_ww%s_nw%s",kk,ee,ww,nw))
		ncores = smart_ncores()
		
		true_KAPPA = mat_true_KAPPA[kk,]; true_KAPPA
		true_ETA = mat_true_ETA[ee,]; true_ETA
		wRHO = wRHO_set[ww]; wRHO
		noiseRHO = noiseRHO_set[nw]; noiseRHO
		vec_MARG = c(TRUE,FALSE)
		vec_TRIM = c(TRUE,FALSE)
		vec_PERM = c(TRUE,FALSE)
		thres_TRIM = 10
		
		# true_KAPPA = c(1,1e-2,1)
	}
	
	# Run replicates
	rr = 1; QQ = length(true_KAPPA); res = c()
	for(rr in seq(RR)){
		# rr = 4
		smart_progress(ii = rr,nn = RR,iter = 5,iter2 = 1e2)
		
		# Simulate dataset
		RHO = gen_true_RHO(wRHO = wRHO,NN = NN,QQ = QQ,RHO = NULL)
		sim = CSeQTL_dataGen(NN = NN,MAF = MAF,true_BETA0 = true_BETA0,
			true_KAPPA = true_KAPPA,true_ETA = true_ETA,true_PHI = true_PHI,
			true_ALPHA = NULL, # all cell types are cis
			RHO = RHO,show = FALSE)
		# true_SNP = sim$true_SNP; dat = sim$dat; XX2 = sim$XX2
		
		if(FALSE){ # Check what tRHO and nRHO look like
			# Generate true_RHO
			tRHO = gen_true_RHO(wRHO = 1,NN = NN,QQ = QQ,RHO = NULL)
			# Distort true_RHO, call it new_RHO
			## nRHO = tRHO + 3
			## nRHO = tRHO * 0.9
			# dd = 0.1; nRHO = tRHO + matrix(runif(NN*QQ,0,dd),NN,QQ)
			dd = 0.3; nRHO = tRHO * matrix(exp(runif(NN*QQ,-dd,dd)),NN,QQ)
			# Normalize new_RHO
			nRHO = nRHO / rowSums(nRHO)
			colnames(nRHO) = paste0("n",colnames(tRHO))
			# Compare true vs new RHO
			smart_scatter(MAT1 = tRHO,MAT2 = nRHO,diagOnly = TRUE)
			cor(tRHO,nRHO)
		
		}
		
		# Add noise to rho
		tRHO = sim$true_RHO
		tRHO = tRHO * exp(matrix(runif(NN*QQ,-1,1)*noiseRHO,NN,QQ))
		tRHO = tRHO / rowSums(tRHO)
		if(FALSE){
			smart_compMATs(MAT1 = sim$true_RHO,
				MAT2 = tRHO,xlab = "truth",ylab = "alt_RHO",
				show_plot = TRUE,show_corr = FALSE)
		}
		
		if(FALSE){# Add influential counts
			n_influ 		= sample(seq(3,5),1)
			idx_influ 	= sample(seq(NN),n_influ,replace = FALSE)
			low_trec 		= seq(1,min(sim$dat$total))
			high_trec		= seq(max(sim$dat$total),max(sim$dat$total) + 50)
			all_trec 		= sim$dat$total
			noise_trec 	= sapply(seq(n_influ),function(zz){
				xx = rbinom(1,1,0.5)
				if(xx == 0){
					yy = sample(low_trec,1)
				} else {
					yy = sample(high_trec,1)
				}
				yy
				})
			all_trec[idx_influ]		= noise_trec
			if(FALSE){
				par(mfrow=c(1,2))
				boxplot(sim$dat$total ~ sim$true_SNP)
				boxplot(all_trec ~ sim$true_SNP)
				par(mfrow=c(1,1))
			}
			sim$dat$total				= all_trec
		}
		
		# Run CSeQTL main function
		sim2_fn = sprintf("%s_rr%s.rds",sim_fn,rr)
		sim$RHO = tRHO
		saveRDS(sim,sim2_fn)
		
		# Store replicate's output, if successful then delete
		sink_fn = sprintf("%s.txt",sim_fn)
		sink(sink_fn)
		
		cat(sprintf("rr = %s\n",rr))
		fout = CSeQTL_full_analysis(TREC = sim$dat$total,
			hap2 = sim$dat$hap2,ASREC = sim$dat$total_phased,
			PHASE = sim$dat$phased,SNP = sim$true_SNP,
			RHO = sim$RHO,XX = sim$XX,log_lib_size = sim$dat$log_lib_size,
			vec_MARG = vec_MARG,vec_TRIM = vec_TRIM,vec_PERM = vec_PERM,
			thres_TRIM = thres_TRIM,ncores = ncores,show = TRUE)
		
		sink()
		unlink(sink_fn)
		unlink(sim2_fn)
		
		res = rbind(res,smart_df(REP = rr,NN = NN,
			MAF = MAF,thres_TRIM = thres_TRIM,fout$res))
		rm(sim,fout)
		rr = rr + 1
	}
	
	# Summarize
	ures = unique(res[,c("NN","MAF","thres_TRIM","MODEL","MARG",
		"TRIM","PERM","CELLTYPE","uASREC")])
	ures$kap_idx = paste(round(true_KAPPA,2),collapse="_")
	ures$eta_idx = paste(round(true_ETA,2),collapse="_")
	ures$wRHO = wRHO; ures$nRHO = noiseRHO; ures$POWER = NA
	ures$consETA = NA; ures$estZERO = NA
	for(ii in seq(nrow(ures))){
		# ii = 1
		idx = which(res$MODEL == ures$MODEL[ii]
			& res$MARG == ures$MARG[ii]
			& res$TRIM == ures$TRIM[ii]
			& res$PERM == ures$PERM[ii]
			& res$CELLTYPE == ures$CELLTYPE[ii]
			& res$uASREC == ures$uASREC[ii])
		ures$POWER[ii] = calc_POWER_2(PVAL = res$PVAL[idx])
		ures$consETA[ii] = mean(res$ETA[idx] == 1)
		ures$estZERO[ii] = mean(res$MU_A[idx] == 0)
	}
	# ures
	# ures[which(ures$PERM==TRUE),]
	
	list(res = res,ures = ures)
	
}
sim_untrim_v_trim = function(NN,MAF,true_BETA0,true_KAPPA,true_ETA,
	true_PHI = 0.1,RHO = NULL,PERMS = 5e1,thres_TRIM = 10,ncores = 1){
	
	if(FALSE){
		NN = 2e2; MAF = 0.3; true_BETA0 = log(200)
		true_KAPPA = c(1,1,1); true_ETA = c(3,1,1)
		true_PHI = 0.1
		RHO = gen_true_RHO(wRHO = 3,NN = NN,QQ = 3)
		PERMS = 5e2
		thres_TRIM = 50; ncores = smart_ncores()
		
	}
	
	# Demonstrate Type 1 Error from simulation under H_A with permutation
	
	# Simulate data
	cat(sprintf("%s: Simulate data ...\n",date()))
	sim = CSeQTL_dataGen(NN = NN,MAF = MAF,true_BETA0 = true_BETA0,
		true_KAPPA = true_KAPPA,true_ETA = true_ETA,true_PHI = true_PHI,
		prob_phased = 0.05,true_ALPHA = NULL,RHO = RHO,show = FALSE)
	
	# Simulate permuted SNPs
	cat(sprintf("%s: Generate permuted SNPs ...\n",date()))
	pSNP = matrix(NA,PERMS,NN)
	for(pp in seq(PERMS)){
		pSNP[pp,] = sample(sim$true_SNP)
	}
	
	# Check trim
	cooksd_res = Rcpp_CSeQTL_cooksD(TREC = sim$dat$total,RHO = sim$true_RHO,
		XX = sim$XX,trim_thres = thres_TRIM,ncores = ncores,show = FALSE)
	
	# Loop thru marg vs cell type-specific
	GS_index = t(c(0,PERMS-1))
	vec_TRIM = c(FALSE,TRUE)
	out = list()
	for(trim in vec_TRIM){
		# trim = vec_TRIM[2]; trim
		cat(sprintf("%s: trim = %s ...\n",date(),trim))
		
		gs_out = Rcpp_CSeQTL_GS(XX = sim$XX,TREC_0 = t(sim$dat$total),
			SNP = t(sim$true_SNP),hap2 = t(sim$dat$hap2),ASREC = t(sim$dat$total_phased),
			PHASE_0 = t(sim$dat$phased*0),RHO = sim$true_RHO,GS_index = t(c(0,0)),
			trim = trim,trim_thres = thres_TRIM,ncores = ncores)
		
		perm_gs_out = Rcpp_CSeQTL_GS(XX = sim$XX,TREC_0 = t(sim$dat$total),
			SNP = pSNP,hap2 = t(sim$dat$hap2),ASREC = t(sim$dat$total_phased),
			PHASE_0 = t(sim$dat$phased*0),RHO = sim$true_RHO,GS_index = GS_index,
			trim = trim,trim_thres = thres_TRIM,ncores = ncores)
		
		tmp_name = sprintf("trim_%s_%s",trim,thres_TRIM)
		out[[tmp_name]] = list(GS = gs_out,pGS = perm_gs_out)
		rm(gs_out,perm_gs_out)
	}
	
	out$trimRES = cooksd_res
	out
}



