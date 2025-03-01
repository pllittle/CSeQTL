# ----------
# CS_eQTL Functions
# ----------

chk_XX = function(XX){
	# XX = cbind(1,c(1,2))
	
	if( !is(XX,"matrix") ) stop("XX should be a matrix")
	if( !all(XX[,1] == 1) ) stop("XX's first column should be all ones")
	
	for(ii in seq(ncol(XX))){
		
		if( ii == 1 ) next
		u_XX = unique(XX[,ii])
		if( all(u_XX %in% c(0,1)) ) next
		
		mean_XX = mean(XX[,ii])
		if( abs(mean_XX) < 1e-10 ) next
		
		stop(sprintf("XX's column %s is not mean zero",ii))
		
	}
	
}

# ----------
# CSeQTL Analysis Functions
# ----------

#' @title CSeQTL_full_analysis
#' @description Performs marginal and cell type-specific eQTL analysis
#'	with both CSeQTL and OLS models for simulation and comparative purposes.
#'
#' @param TREC An integer vector of total read counts.
#' @param hap2 An integer vector of second haplotype counts
#' @param ASREC An integer vector of total haplotype counts
#' @param PHASE An integer vector of 0s and 1s denoting if a subject 
#'	has available haplotype counts.
#' @param SNP An integer vector of phased genotypes coded 0 (AA),
#'	1 (AB), 2 (BA), 3 (BB), and 5 (NA).
#' @param RHO A numeric matrix of cell type proportions. Rows 
#'	correspond to subjects and columns correspond to cell types.
#' @param XX A numeric design matrix of baseline covariates 
#'	including the intercept in the first column and centered 
#'	continuous covariates.
#' @param log_lib_size A positive numeric vector of log transformed
#'	library sizes per subject.
#' @inheritParams CSeQTL_oneExtremeSim
#' @inheritParams CSeQTL_dataGen
#' @return A R list containing multiple objects. \code{res} is a R dataframe
#'	containing the model fitted, marginal model indicator, TRIM indicator, 
#'	permutation indicator, cell types, utilized allele-specific reads indicator,
#'	A and B allele-specific expression, eta/eqtl fold change estimate, p-value.
#'	\code{out} contains lists of detailed estimates per model fit.
#' @export
CSeQTL_full_analysis = function(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,
	log_lib_size,vec_MARG = c(TRUE,FALSE),vec_TRIM = c(TRUE,FALSE),
	vec_PERM = c(TRUE,FALSE),thres_TRIM = 10,ncores = 1,
	show = TRUE){
	
	chk_XX(XX = XX)
	
	NN = nrow(XX); QQ = ncol(RHO)
	marg_RHO = matrix(1,NN,1)
	
	# Permuted SNPs
	perm_idx = sample(seq(NN))
	perm_SNP = SNP[perm_idx]
	
	run_TRIM = any(vec_TRIM == TRUE)
	run_MARG = any(vec_MARG == TRUE)
	
	# Calculate Cooks' Distance and standardize
	if( run_TRIM ){
		cooksd_res = Rcpp_CSeQTL_cooksD(TREC = TREC,RHO = RHO,XX = XX,
			trim_thres = thres_TRIM,ncores = ncores,show = FALSE)
		trim_TREC_cts = cooksd_res[,4]
		
		if( run_MARG ){
			cooksd_res = Rcpp_CSeQTL_cooksD(TREC = TREC,RHO = marg_RHO,XX = XX,
				trim_thres = thres_TRIM,ncores = ncores,show = FALSE)
			trim_TREC_marg_propNo = cooksd_res[,4]
		} else {
			trim_TREC_marg_propNo = NULL
		}
		
		trim_TREC = list(cts = trim_TREC_cts,propNo = trim_TREC_marg_propNo)
	} else {
		trim_TREC = NULL
	}
	
	# Optimize
	out = list(); res = c()
	for(model in c("cseqtl","ols")){
	for(marg in vec_MARG){
		# marg = FALSE
	for(trim in vec_TRIM){
	for(perm in vec_PERM){
		if( perm ){
			all_use_asrec = FALSE
		} else {
			all_use_asrec = c(TRUE,FALSE)
			if( all(PHASE == 0) ) all_use_asrec = FALSE
		}
	for(use_asrec in all_use_asrec){
		# model = c("cseqtl","ols","ols2")[3]; marg = c(TRUE,FALSE)[2]
		# trim = c(TRUE,FALSE)[2]; perm = c(TRUE,FALSE)[2]
		
		if( marg == TRUE ){
			final_RHO = marg_RHO
			colnames(final_RHO) = "MARG"
		} else {
			final_RHO = RHO
		}
		if( trim == TRUE ){
			if( marg == TRUE ){
				final_TREC = trim_TREC_marg_propNo
			} else {
				final_TREC = trim_TREC_cts
			}
		} else {
			final_TREC = TREC
		}
		if( perm ){
			final_SNP 	= perm_SNP
			final_PHASE = PHASE*0
		} else {
			final_SNP 	= SNP
			if( use_asrec == TRUE ){
				final_PHASE = PHASE
			} else {
				final_PHASE = PHASE*0
			}
		}
		
		out_nm = sprintf("trim%s_perm%s_marg%s_%s_asrec%s",
			trim,perm,marg,model,use_asrec)
		out_nm
		
		if( model == "cseqtl" ){
			GS_index 	= matrix(c(0,0),1,2); GS_index
			if( show ) message(sprintf("trim = %s; perm = %s; useASREC = %s\n",trim,perm,use_asrec),appendLF = FALSE)
			out_list = Rcpp_CSeQTL_GS(XX = XX,TREC_0 = t(final_TREC),SNP = t(final_SNP),
				hap2 = t(hap2),ASREC = t(ASREC),PHASE_0 = t(final_PHASE),RHO = final_RHO,
				GS_index = GS_index,trim = FALSE,show = show,prompt = show,ncores = 1)
			
			log_ETA = log(out_list$ETA[1,])
			MU = rbind(out_list$MU_A,out_list$MU_B); MU
			LRT_trec = cbind(out_list$LRT_trec[1,],1,1-pchisq(out_list$LRT_trec[1,],df=1))
			colnames(LRT_trec) = paste0(c("LRT","DF","P"),"_trec")
			LRT_trecase = cbind(out_list$LRT_trecase[1,],1,1-pchisq(out_list$LRT_trecase[1,],df=1))
			colnames(LRT_trecase) = paste0(c("LRT","DF","P"),"_trecase")
			LRT_cistrans = cbind(out_list$LRT_cistrans[1,],1,1-pchisq(out_list$LRT_cistrans[1,],df=1))
			colnames(LRT_cistrans) = paste0(c("LRT","DF","P"),"_cistrans")
			tmp_res = cbind(LRT_trec,LRT_trecase,LRT_cistrans)
			rownames(tmp_res) = colnames(final_RHO)
			tmp_res = smart_df(tmp_res)
			tmp_res$P_final = apply(tmp_res[,c("P_trec","P_trecase","P_cistrans")],1,
				function(xx) ifelse(xx[3] < 0.05,xx[1],xx[2]))
			# if( show ) print(tmp_res)
			out[[out_nm]] = list(res = tmp_res,log_ETA = log_ETA,MU = MU)
			
			res = rbind(res,smart_df(MODEL = model,
				MARG = marg,TRIM = trim,
				PERM = perm,CELLTYPE = colnames(final_RHO),
				uASREC = use_asrec,MU_A = out_list$MU_A[1,],
				MU_B = out_list$MU_B[1,],ETA = out_list$ETA[1,],
				# PVAL = tmp_res$P_final
				PVAL = tmp_res$P_trecase
				))
			rm(tmp_res,log_ETA,MU,out_list)
			
		} else if( model %in% c("ols") && use_asrec == FALSE ){
			dat = smart_df(total = TREC,log_lib_size = log_lib_size)
			
			if( model == "ols" ){
				out_ols_list = CSeQTL_linearTest(input = dat,
					XX = XX,RHO = final_RHO,SNP = final_SNP,
					trim = trim,thres_TRIM = thres_TRIM,
					show_plot = FALSE)
				out[[out_nm]] = out_ols_list
				
				res = rbind(res,smart_df(MODEL = model,
					MARG = marg,TRIM = trim,
					PERM = perm,CELLTYPE = colnames(final_RHO),
					uASREC = use_asrec,
					MU_A = NA,MU_B = NA,ETA = out_ols_list$out_df$EST,PVAL = out_ols_list$out_df$PVAL))
				rm(out_ols_list)
			
			}
			
		}
		
	}}}}}
	if( show ) print(res)
	
	return(list(res = res,out = out,trim_TREC = trim_TREC))
}
check_cont = function(xx){
	if( all(xx == 1) ){
		FALSE # intercept
	} else if( all(xx %in% c(0,1)) ){
		FALSE # categorical variables
	} else {
		TRUE
	}
}
proc_LRT = function(LRT,DF=1){
	tmp_df = smart_df(LRT = c(LRT),DF = DF)
	tmp_df$PVAL = 1 - pchisq(tmp_df$LRT,df = tmp_df$DF)
	tmp_df
}
polish_res = function(gs_out,celltypes,rdnum = 5){
	LRT_nms = c("LRT","DF","PVAL")
	ETA_trec = apply(cbind(c(gs_out$MU_A),c(gs_out$MU_B)),1,function(xx)
		ifelse(xx[1] == xx[2],1,xx[2]/xx[1]))
	list(test = round(cbind(smart_names(proc_LRT(gs_out$LRT_trec),
			ROW = celltypes,COL = paste0(LRT_nms,"_trec"))[,-2],
		smart_names(proc_LRT(gs_out$LRT_trecase),
			ROW = celltypes,COL = paste0(LRT_nms,"_trecase"))[,-2],
		smart_names(proc_LRT(gs_out$LRT_cistrans),
			ROW = celltypes,COL = paste0(LRT_nms,"_cistrans"))[,-2],
		MU_A = c(gs_out$MU_A),MU_B = c(gs_out$MU_B),
		ETA_trec = ETA_trec,
		ETA_final = c(gs_out$ETA)),rdnum),
		LL = gs_out$LL)
}

#' @title CSeQTL_linearTest
#' @description Runs marginal and cell type-specific 
#'	analysis using ordinary least squares.
#'
#' @param input A data.frame containing columns \code{total}
#'	(total read counts) and \code{log_lib_size} (log transformed
#'	library size).
#' @param YY Default is \code{NULL}. By default, the outcome
#'	for OLS is the inverse rank quantile normalized TReC 
#'	after library size correction. Otherwise, the user can input
#'	their own transformed outcome as a numeric vector.
#' @param MARG Boolean value. Set to \code{TRUE} to fit the
#'	marginal OLS model. Default is set to \code{FALSE} to fit
#'	the cell type-specific interaction OLS model.
#' @param impute_geno Boolean value. Default is set to \code{FALSE}
#'	to only analyze subjects with non-missing genotype.
#'	If \code{TRUE}, missing genotypes are imputed with the 
#'	mean of non-missing genotypes mimicing \code{MatrixEQTL}.
#' @param trim Boolean value set to \code{FALSE} by default to prevent
#'	outcome trimming. If \code{TRUE}, the OLS model will be fitted
#'	without covariates containing SNP genotype to calculate each 
#'	subject's Cooks' distance.
#' @param show_plot Boolean value set to \code{TRUE} to visualize
#'	boxplot or interactions.
#' @param main_plot Character string for the visual's main title.
#' @param CTs Set to \code{NULL} by default. If \code{NULL} and 
#'	\code{show_plot = TRUE}, all cell type interaction plots are shown.
#'	Otherwise a subset of cell types can be displayed using an integer vector
#'	of columns or character vector of column names of \code{RHO} can 
#'	be provided.
#' @inheritParams CSeQTL_full_analysis
#' @return A R list containing \code{lm()} output for \code{lm_out}, 
#'	a R dataframe for \code{out_df} containing regression estimates, 
#'	standard errors, p-values. \code{res_trim} provides a summary of trimmed 
#'	results over a grid of cut-off values and number of samples trimmed.
#'	\code{cooksd} is a numeric vector of median shifted and MAD scaled 
#'	Cook's distances per sample. \code{prop_trim}, a numeric value for 
#'	number of samples with outcome values trimmed for the user-specified
#'	\code{thres_TRIM} value.
#'
#' @export
CSeQTL_linearTest = function(input,XX,RHO,SNP,YY = NULL,MARG = FALSE,
	impute_geno = FALSE,trim = FALSE,thres_TRIM = 10,show_plot = TRUE,
	main_plot = "",CTs = NULL){
	
	if(FALSE){
		input = dat; XX = XX2; RHO = true_RHO; SNP = true_SNP
		RHO = matrix(1,nrow(tRHO),1)
		
		input = sim$dat; XX = sim$XX; RHO = sim$true_RHO; SNP = sim$true_SNP
		
		input = tmp_dat; XX = inputs$XX; RHO = inputs$RHO; SNP = SNP
		
		YY = NULL; MARG = FALSE; impute_geno = FALSE
		trim = FALSE; thres_TRIM = 10; show_plot = FALSE
		
	}
	
	# Code to analyze data with linear model with 
	#		genotype and cell type proportions with interaction
	
	if( MARG == TRUE ){
		RHO = matrix(rowSums(RHO),ncol = 1)
		colnames(RHO) = "Bulk"
	}
	
	chk_XX(XX = XX)
	
	# Constants
	NN = nrow(RHO)
	QQ = ncol(RHO)
	nX = ncol(XX)
	num_covars = nX + 1 + (QQ - 1) + (QQ - 1)
	colnames(XX) = paste0("X",seq(ncol(XX)))
	
	# Contrast matrix
	KK = matrix(0,QQ,num_covars)
	if( QQ > 1 ){
		colnames(KK) = c(colnames(XX),"SNP",colnames(RHO)[-1],
			paste0("SNP_CT.",colnames(RHO)[-1]))
	} else {
		colnames(KK) = c(colnames(XX),"SNP")
	}
	KK[,"SNP"] = 1 # eQTL 1st cell type
	if( QQ > 1 ){
		if( QQ == 2 ){
			KK[-1,grepl("SNP_CT",colnames(KK))] = 1
		} else {
			diag(KK[-1,grepl("SNP_CT",colnames(KK))]) = 1
		}
	}
	
	# SNP code as 0/1/2
	SNP2 = SNP
	SNP2[SNP2 == 2] = 1
	SNP2[SNP2 == 3] = 2
	SNP2[SNP2 == 5] = NA
	na_idx = is.na(SNP2)
	if( impute_geno ){
		SNP2[na_idx] = mean(SNP2,na.rm = TRUE)
	}
	nm_idx = !is.na(SNP2)
	if( is.null(YY) ){
		input = input[nm_idx,,drop = FALSE]
	} else {
		YY = YY[nm_idx]
	}
	SNP2 	= SNP2[nm_idx]
	RHO 	= RHO[nm_idx,,drop = FALSE]
	XX 		= XX[nm_idx,,drop = FALSE]
	# table(SNP,SNP2)
	
	# Calculate interactions and final design matrix
	NN = length(SNP2)
	if(QQ == 1){
		iXX = smart_df(SNP = SNP2)
	} else if(QQ >= 2){
		iXX = apply(RHO[,-1,drop = FALSE],2,function(xx) xx * SNP2)
		iXX = smart_df(iXX)
		names(iXX) = paste0("SNP_CT",seq(2,QQ))
		iXX = smart_df(SNP = SNP2,RHO[,-1,drop = FALSE],iXX)
	}
	rownames(iXX) = NULL
	fin_XX = smart_df(XX[,-1,drop = FALSE],iXX)
	rm(iXX)
	PP = ncol(XX)
	
	# Calculate transformed outcome
	if( is.null(YY) ){
		# First correct for library size
		YY = input$total / exp(input$log_lib_size)
		
		# Inverse rank quantile normalize
		YY = qnorm(p = rank(YY) / (NN + 1),mean = 0,sd = 1)
	} else {
		YY = YY
	}
	
	# Calculate standardized cooks' distance
	if( QQ == 1 ){
		tmp_df = smart_df(XX[,-1,drop = FALSE])
	} else {
		tmp_df = smart_df(XX[,-1,drop = FALSE],RHO[,-1,drop = FALSE])
	}
	lm_out 	= lm(YY ~ .,data = tmp_df)
	resi 		= residuals(lm_out)
	pred 		= predict(lm_out)
	cooksd 	= cooks.distance(lm_out)
	cooksd 	= (cooksd - median(cooksd)) / mad(cooksd)
	prop_trim = mean(cooksd >= thres_TRIM)
	
	# Determine final outcome (trimmed or untrimmed)
	YY_final = YY
	idx_high_cooksd = which(cooksd >= thres_TRIM)
	if( trim ){
		YY_trim = YY
		if( length(idx_high_cooksd) > 0 ){
			YY_trim[idx_high_cooksd] = pred[idx_high_cooksd]
		}
		YY_final = YY_trim
		
		# Assess prop subjs trimmed
		res_trim = c()
		for(thres in seq(5,30)){
			idx = which(cooksd >= thres)
			tmp_df = smart_df(trim_thres = thres,
				prop_trim = length(idx) / NN)
			res_trim = rbind(res_trim,tmp_df)
		}
	
	} else {
		res_trim = NULL
		prop_trim = 0
	}
	
	# Full model
	lm_full = lm(YY_final ~ .,data = fin_XX)
	# summary(lm_full)
	pred = predict(lm_full)
	resi = residuals(lm_full)
	
	# Hypothesis testing
	tout = glht(model = lm_full,linfct = KK,alternative = "two.sided")
	out_df = smart_df(CT = colnames(RHO))
	tout2 = summary(tout,test = adjusted("none"))$test
	out_df$EST 	= as.numeric(tout2$coefficients)
	out_df$SE 	= as.numeric(tout2$sigma)
	out_df$PVAL = as.numeric(tout2$pvalues)
	rm(tout,tout2)
	
	if( show_plot ){ # Visualize interaction plots
		# Regress out baseline covariates
		lm_out = lm(YY_final ~ .,data = smart_df(XX[,-1]))
		# summary(lm_out)
		SNP_cols = c("deepskyblue2","darkorange","forestgreen")
		idx_012 = which(SNP2 %in% c(0,1,2))
		
		oma = rep(0,4)
		if( main_plot != "" ) oma = c(0,0,2,0)
		
		oldpar = par(no.readonly = TRUE)
		on.exit(par(oldpar))
		
		par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)
		if( QQ == 1 ){
			SNP2_char = "black"
			SNP2_char[SNP2 == 0] = "AA"
			SNP2_char[SNP2 == 1] = "AB"
			SNP2_char[SNP2 == 2] = "BB"
			boxplot(lm_out$residuals[idx_012] ~ SNP2_char[idx_012],
				col = SNP_cols,pch = 16,xlab = "Genotype",
				ylab = "Expression (confounder-adjusted)",
				cex.lab = 1.3,main = sprintf("P = %s",round(out_df$PVAL,3)),
				cex.axis = 1.3)
			
		} else {
			alpha = 0.5
			cex_high_cook = 2
			if( is.null(CTs) ){
				fRHO = RHO
			} else {
				fRHO = RHO[,CTs,drop = FALSE]
			}
			nplots = ncol(fRHO) + 2
			num_rows = ceiling(sqrt(nplots))
			num_cols = ceiling(nplots/num_rows)
			tmp_vec = rep(0,(num_rows+1)*num_cols)
			tmp_vec[seq(nplots)] = seq(nplots)
			mat = matrix(tmp_vec,nrow = num_rows + 1,
				ncol = num_cols,byrow = TRUE)
			mat[nrow(mat),] = nplots + 1
			# mat
			hts = rep(2.5,num_rows+1); hts[num_rows+1] = 1
			hts = hts / sum(hts); hts
			layout(mat,heights = hts)
			# layout.show(max(mat))
			par(mar = c(4.4,4,4,0.5),oma = oma)
			for(ct in colnames(fRHO)){
				# ct = colnames(fRHO)[1]
				sub_main = sprintf("P = %s",
					formatC(out_df$PVAL[out_df$CT==ct],3))
				plot(fRHO[,ct],lm_out$residuals,pch = 19,
					col = SNP_cols[SNP2+1],bty = "n",
					xlab = paste0(ct," Proportion"),
					ylab = "Expression",cex.axis = 1.3,
					cex.lab = 1.3,main = sub_main)
				for(gg in c(0,1,2)){
					gg_idx = which(SNP2 == gg)
					if( length(gg_idx) > 1 ){
						tmp_lm = lm(lm_out$residuals[gg_idx] ~ fRHO[gg_idx,ct])
						abline(a = tmp_lm$coefficients[1],
							b = tmp_lm$coefficients[2],lwd = 2,
							lty = 2,col = SNP_cols[gg+1])
					}
				}
				points(fRHO[idx_high_cooksd,ct],
					lm_out$residuals[idx_high_cooksd],
					pch = 1,lwd = 2,cex = cex_high_cook)
			}
			
			plot(pred,YY_final,bty = "n",pch = 19,
				xlab = "Predicted Y",ylab = "Observed Y",
				cex.axis = 1.3,cex.lab = 1.3,col = SNP_cols[SNP2+1])
			points(pred[idx_high_cooksd],YY_final[idx_high_cooksd],
				pch = 1,lwd = 2,cex = cex_high_cook)
			abline(a = 0,b = 1,lty = 2)
			
			plot(pred,resi,bty = "n",pch = 19,
				xlab = "Predicted Y",ylab = "Residuals",
				cex.axis = 1.3,cex.lab = 1.3,col = SNP_cols[SNP2+1])
			points(pred[idx_high_cooksd],resi[idx_high_cooksd],
				pch = 1,lwd = 2,cex = cex_high_cook)
			abline(h = 0,lty = 2)
			
			par(mar = rep(0,4))
			plot(0,0,xlim = c(0,1),ylim = c(0,1),type = "n",
				xlab = "",ylab = "",axes = FALSE)
			legend("center",legend = c("AA","AB","BB"),col = SNP_cols,
				pch = 16,cex = 2,bty = "n",ncol = 3,title = "Genotype")
			if( main_plot != "" ) mtext(main_plot,outer = TRUE,cex = 1.2)
		}
		par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1,oma = rep(0,4))
		
	}
	
	# Output
	return(list(lm_out = lm_full,out_df = out_df,
		res_trim = res_trim,cooksd = cooksd,
		prop_trim = prop_trim))
}

#' @title CSeQTL_run_MatrixEQTL
#' @param TREC A matrix of integer total read counts, rows are genes
#'	with row labels with syntax "gene_name:chrom:start:end", 
#'	columns are samples with column labels.
#' @param RD A positive numeric vector of library sizes per sample.
#' @param SNP An integer matrix, rows are genomic loci with row labels
#'	such as "chrom:pos", columns correspond to samples with column labels.
#' @param out_cis_fn A full path filename string to store
#'	MatrixEQTL output
#' @param cisDist A positive integer for number of SNPs to
#'	include relative to the gene body
#' @inheritParams CSeQTL_full_analysis
#' @return The MatrixEQTL output file/R dataframe generated containing genes, 
#'	SNPs, associated p-values, effect sizes, etc.
#' @export
CSeQTL_run_MatrixEQTL = function(TREC,RD,XX,SNP,out_cis_fn,cisDist = 1e6){
	
	normscore <- function(vec) {
		len  = length(na.omit(vec))+1
		rank = rank(na.omit(vec))
		ties = (rank - floor(rank)) > 0
		new.vec = vec[!is.na(vec)] 
		new.vec[!ties]=qnorm(rank[!ties]/len)
		new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
		vec[!is.na(vec)] = new.vec
		vec
	}

	snps = SNP
	snps[snps == 2] = 1
	snps[snps == 3] = 2
	snps[snps == 5] = NA
	colnames(snps) = colnames(TREC)
	snps = SlicedData$new(snps)
	
	gene = t(apply(TREC,1,function(xx) xx / RD))
	gene = t(apply(gene,1,normscore))
	gene = SlicedData$new(gene)
	
	chk_XX(XX = XX)
	
	cvrt = t(XX[,-1])
	colnames(cvrt) = colnames(TREC)
	cvrt = SlicedData$new(cvrt)
	snpspos = smart_df(snpid = rownames(SNP))
	snpspos$chr = sapply(snpspos$snpid,function(xx)
		strsplit(xx,":")[[1]][1],USE.NAMES = FALSE)
	snpspos$chr = gsub("chr","",snpspos$chr)
	snpspos$pos = as.integer(sapply(snpspos$snpid,function(xx)
		strsplit(xx,":")[[1]][2],USE.NAMES = FALSE))
	if( ncol(snpspos) != 3 ){
		print(str(snpspos))
		print(head(snpspos))
	}
	
	genepos = smart_df(geneid = rownames(TREC))
	genepos$chr = sapply(rownames(TREC),function(xx)
		strsplit(xx,":")[[1]][2],USE.NAMES = FALSE)
	genepos$chr = gsub("chr","",genepos$chr)
	genepos$left = sapply(rownames(TREC),function(xx)
		strsplit(xx,":")[[1]][3],USE.NAMES = FALSE)
	genepos$left = as.integer(genepos$left)
	genepos$right = sapply(rownames(TREC),function(xx)
		strsplit(xx,":")[[1]][4],USE.NAMES = FALSE)
	genepos$right = as.integer(genepos$right)
	
	me = Matrix_eQTL_main(snps = snps,gene = gene,
		cvrt = cvrt,output_file_name.cis = out_cis_fn,
		output_file_name = NULL,pvOutputThreshold = 1,
		useModel = modelLINEAR,
		snpspos = snpspos,genepos = genepos,
		pvOutputThreshold.cis = 1,errorCovariance = numeric(),
		cisDist = cisDist,verbose = TRUE,pvalue.hist = TRUE,
		min.pv.by.genesnp = FALSE,noFDRsaveMemory = FALSE)
	me_res = fread(out_cis_fn,data.table = FALSE)
	unlink(out_cis_fn)
	
	return(me_res)
}

#' @title OLS_sim
#' @inheritParams CSeQTL_oneExtremeSim
#' @inheritParams CSeQTL_linearTest
#' @param showRHO Boolean for seeing a preview of how 
#'	proportions are distributed.
#' @return Simulation results for OLS-based approach with per replicate
#'	and per cell type metrics including p-values, OLS effect sizes,
#'	and proportion of samples with trimmed outcomes.
#' @export
OLS_sim = function(NN,MAF,true_BETA0,true_KAPPA,true_ETA,wRHO,RR,
	trim = FALSE,thres_TRIM = 10,showRHO = TRUE){
	
	vec_cols = c("deepskyblue2","darkorange","forestgreen")
	QQ = length(true_KAPPA)
	
	# Generate specific RHO scenario
	if( showRHO ){
		RHO = gen_true_RHO(wRHO = wRHO,NN = NN,QQ = QQ,RHO = NULL)
		boxplot(RHO,col = vec_cols,pch = 16,
			xlab = "Cell Type",ylab = "Proportion",
			main = sprintf("Scenario %s",wRHO))
	}

	PVAL 	= matrix(NA,RR,QQ)
	EST		= matrix(NA,RR,QQ)
	SE		= matrix(NA,RR,QQ)
	EQTL 	= matrix(NA,RR,QQ)
	colnames(PVAL) 	= paste0("CT",seq(QQ))
	colnames(EST) 	= paste0("CT",seq(QQ))
	colnames(SE) 		= paste0("CT",seq(QQ))
	colnames(EQTL) 	= paste0("CT",seq(QQ))
	pTRIM = rep(NA,RR)

	for(rr in seq(RR)){
		smart_progress(ii = rr,nn = RR)
		
		# Sim proportions
		RHO = gen_true_RHO(wRHO = wRHO,NN = NN,QQ = QQ,RHO = NULL)
		
		# Sim covariates, SNP, TReC
		sim = CSeQTL_dataGen(NN = NN,MAF = MAF,true_BETA0 = true_BETA0,
			true_KAPPA = true_KAPPA,true_ETA = true_ETA,true_PHI = 0.1,
			RHO = RHO,show = FALSE)
		
		# Hypothesis testing
		ols_out = CSeQTL_linearTest(input = sim$dat,XX = sim$XX,
			RHO = sim$true_RHO,SNP = sim$true_SNP,MARG = FALSE,
			trim = trim,thres_TRIM = thres_TRIM,show_plot = FALSE)
		
		# Hypothesis test results
		PVAL[rr,] = ols_out$out_df$PVAL
		EST[rr,] 	= ols_out$out_df$EST
		SE[rr,]		= ols_out$out_df$SE
		EQTL[rr,] = ols_out$out_df$EST / ols_out$out_df$SE
		pTRIM[rr] = ols_out$prop_trim
		rm(RHO,sim,ols_out)
	}
	
	hist(pTRIM,breaks = 30,col = "deepskyblue",
		xlab = "Prop Subjs Trimmed",main = "",cex.lab = 1.3)
	
	# Check out p-value distribution and Type 1 Error
	oldpar = par(no.readonly = TRUE)
	on.exit(par(oldpar))
	
	par(mfrow = c(3,4),oma = c(0,0,2,0))
	for(ct in seq(QQ)){
		# ct = 1
		hist(PVAL[,ct],main = sprintf("%s: Error = %s",colnames(PVAL)[ct],
			round(mean(PVAL[,ct] <= 0.05),3)),freq = FALSE,cex.lab = 1.3,
			xlab = "P-value",col = "deepskyblue",breaks = 20)
		abline(h = 1,lty = 2,lwd = 2,col = "red")
		
		hist(EST[,ct],breaks = 20,freq = FALSE,xlab = "EST",
			main = colnames(EST)[ct],col = "red",xlim = range(c(EST)),
			cex.lab = 1.3)
		
		hist(SE[,ct],breaks = 20,freq = FALSE,xlab = "SE",
			main = colnames(SE)[ct],col = "red",xlim = range(c(SE)),
			cex.lab = 1.3)
		
		hist(EQTL[,ct],breaks = 20,freq = FALSE,xlab = "EST / SE",
			main = colnames(EQTL)[ct],col = "red",xlim = range(c(EQTL)),
			cex.lab = 1.3)
	}
	
	if( !trim ){
		plot_text = sprintf("Scenario %s;Trim = %s;Num Reps = %s",wRHO,trim,RR)
	} else {
		plot_text = sprintf("Scenario %s;Trim = %s(thres = %s);Num Reps = %s",
			wRHO,trim,thres_TRIM,RR)
	}
	mtext(plot_text,outer = TRUE,cex = 1.4)
	par(mfrow = c(1,1),oma = rep(0,4))
	
	return(list(PVAL = PVAL,EQTL = EQTL,pTRIM = pTRIM))
}

#' @title plot_RHO
#' @param ... Arguments for \code{plot(x,y,bty = "n",...)}.
#' @inheritParams CSeQTL_linearTest
#' @export
plot_RHO = function(RHO,main_plot = "",...){
	QQ = ncol(RHO)
	num_pairs = choose(QQ,2)
	num_rows = ceiling(sqrt(num_pairs))
	num_cols = ceiling(num_pairs / num_rows)
	
	oma = rep(0,4)
	if( main_plot != "" ) oma = c(0,0,2,0)
	
	oldpar = par(no.readonly = TRUE)
	on.exit(par(oldpar))
	
	par(mfrow = c(num_rows,num_cols),oma = oma,mar = c(4,4,0.1,0.1))
	for(CT_1 in seq(QQ)){
	for(CT_2 in seq(QQ)){
		if( CT_1 < CT_2 ){
			plot(smart_df(RHO[,c(CT_1,CT_2)]),bty = "n",
				...)
		}
	}}
	if( main_plot != "" ) mtext(main_plot,outer = TRUE,cex = 1.2)
	par(mfrow = c(1,1),oma = rep(0,4),mar = c(5,4,4,2)+0.1)
	
}

#' @title CSeQTL_smart
#' @description A function to understand the novelties within CSeQTL. 
#'	The user can experiment with trimming TReC, assess parameter estimation, 
#'	perform hypothesis testing, actively constrain cell type-specific parameters,
#'	when analyzing a single gene/SNP pair.
#'
#' @param upPHI A value of 0 or 1 indicating if a Poisson or
#'	Negative binomial distribution is fitted, respectively.
#' @param upKAPPA An integer vector of zeroes and ones where 
#'	\code{length(upKAPPA) == ncol(RHO)}. A requirement is \code{upKAPPA[1] = 1}.
#'	Cell types with indices equal to one have the baseline fold change
#'	between the q-th and reference cell type are estimated. Otherwise
#'	that cell type's parameter is constrained to zero.
#' @param upETA An integer vector of zeroes and ones where
#'	\code{length(upETA) == ncol(RHO)} indicating which
#'	cell types eQTL parameters are estimated or constrained to their null.
#' @param upPSI A value of 0 or 1 indicating if a Binomial or
#'	Beta-binomial distribution is fitted, respectively.
#' @param upALPHA An integer vector of zeroes and ones where
#'	\code{length(upALPHA) == ncol(RHO)} indicating which
#'	cell types' cis-trans parameters are estimated or 
#'	constrained to their null.
#' @param iFullModel A boolean that if set to \code{TRUE} will 
#'	determine the submodel of the full model, based by upKAPPA, 
#'	upETA, and upALPHA parameters, that can be estimated with stability.
#' @param trim Boolean value. If \code{TRUE}, the TReC model will be fitted
#'	without SNP genotype to calculate each subject's Cooks' distance.
#' @param hypotest A boolean to perform eQTL significance 
#'	testing and cis-trans eQTL testing.
#' @param swap A boolean to determine if the reference cell type should
#'	be swapped with the cell type with highest TReC across alleles.
#' @param numAS A positive integer to determine if a subject has
#'	enough total haplotype counts.
#' @param numASn A positive integer to determine how many subjects
#'	have at least \code{numAS} to use the haplotype counts.
#' @param numAS_het A positive integer to determine how many subjects
#'	with at least \code{numAS} are heterozygous (AB or BA). If
#'	\code{sum(PHASE == 1 & ASREC >= numAS & (SNP == 1 | SNP == 2)) >= numAS_het},
#'	those subjects haplotype counts will be used for TReCASE and cis/trans 
#'	testing and estimation.
#' @param cistrans_thres A numeric value to determine the 
#'	cis/trans test p-value cutoff.
#' @param gr_eps A numeric value to determine if convergence
#'	is achieved based on the L2 norm of the gradient.
#' @param conv_eps A numeric value to determine if convergence
#'	is achieved based on the L2 norm of the product of the 
#'	inverse hessian and gradient.
#' @inheritParams CSeQTL_full_analysis
#'
#' @export
CSeQTL_smart = function(TREC,hap2,ASREC,PHASE,SNP,RHO,XX,
	upPHI,upKAPPA,upETA,upPSI,upALPHA,iFullModel = FALSE,trim = FALSE,
	thres_TRIM = 10,hypotest = TRUE,swap = TRUE,numAS = 5,
	numASn = 5,numAS_het = 5,cistrans_thres = 0.01,gr_eps = 1e-2,
	conv_eps = 1e-3,ncores = 1,show = FALSE){
	
	if(FALSE){
		TREC = sim$dat$total; hap2 = sim$dat$hap2; ASREC = sim$dat$total_phased
		PHASE = sim$dat$phased; SNP = sim$true_SNP; RHO = sim$true_RHO; XX = sim$XX
		upPHI = 0; upPSI = 0
		QQ = ncol(RHO); upKAPPA = rep(1,QQ); upETA = upKAPPA; upALPHA = upKAPPA
		
		iFullModel = FALSE; trim = FALSE; thres_TRIM = 10; hypotest = TRUE;
		swap = TRUE; numAS = 5; numASn = 5; numAS_het = 5; ncores = 1
		show = !TRUE
		
	}
	
	# Input checks
	if( is.null(colnames(RHO)) ) stop("Set colnames(RHO)!")
	if( is.null(colnames(XX)) ) stop("Set colnames(XX)!")
	
	chk_XX(XX = XX)
	
	# Constants
	PP = ncol(XX)
	QQ = ncol(RHO)
	GI = Rcpp_calc_GI(PP = PP,QQ = QQ)
	log_phi_idx = GI[2,1] + 1
	kap_idx = seq(GI[3,1],GI[3,2]) + 1
	eta_idx = seq(GI[4,1],GI[4,2]) + 1
	log_psi_idx = GI[5,1] + 1
	alp_idx = seq(GI[6,1],GI[6,2]) + 1
	np = GI[6,2] + 1
	vname = c(colnames(XX),"logPHI")
	if( QQ > 1 ) vname = c(vname,paste0("logKAP.",colnames(RHO)[-1]))
	vname = c(vname,paste0("logETA.",colnames(RHO)))
	vname = c(vname,"logPSI")
	vname = c(vname,paste0("logALP.",colnames(RHO)))
	
	# Enforce nested constraints
	if( length(upKAPPA) != QQ ) stop("length(upKAPPA) != ncol(RHO)")
	if( length(upETA) != QQ ) 	stop("length(upETA) != ncol(RHO)")
	if( length(upALPHA) != QQ ) stop("length(upALPHA) != ncol(RHO)")
	if( !(upPHI %in% c(0,1)) ) stop("upPHI should be 0 or 1")
	if( !(upPSI %in% c(0,1)) ) stop("upPSI should be 0 or 1")
	upKAPPA[1] = 1
	upETA 	= upKAPPA * upETA
	upALPHA = upETA * upALPHA
	
	# Constraining parameters
	upPARS_0 = rep(1,np)
	names(upPARS_0) = vname
	upPARS_0[log_phi_idx] = upPHI
	upPARS_0[1] = upKAPPA[1]
	if( QQ > 1 ) upPARS_0[kap_idx] = upKAPPA[-1]
	upPARS_0[eta_idx] = upETA
	upPARS_0[log_psi_idx] = upPSI
	upPARS_0[alp_idx] = upALPHA
	
	# Fit model and hypothesis testing
	bs_out = Rcpp_CSeQTL_BFGS_smart(TREC = TREC,hap2 = hap2,
		ASREC = ASREC,PHASE_0 = PHASE,SNP = SNP,RHO = RHO,
		XX = XX,upPARS_0 = upPARS_0,iFullModel = iFullModel,
		trim = trim,trim_thres = thres_TRIM,hypotest = hypotest,
		swap = swap,numAS = numAS,numASn = numASn,numAS_het = numAS_het,
		gr_eps = gr_eps,conv_eps = conv_eps,hess_shift = 1e-5,
		ncores = ncores,show = show)
	
	# ETA
	colnames(bs_out$ETA) = colnames(RHO)
	rownames(bs_out$ETA) = c("TReC","TReCASE")
	
	# expected TReC per cell type and allele
	colnames(bs_out$MU) = colnames(RHO)
	rownames(bs_out$MU) = paste0("allele.",c("A","B"))
	# bs_out$MU

	# Hypothesis Testing (TReC-only)
	colnames(bs_out$HYPO$LRT_trec) = c("LRT","DF","PVAL")
	rownames(bs_out$HYPO$LRT_trec) = colnames(RHO)
	# bs_out$LRT_trec

	# Hypothesis Testing (TReCASE)
	colnames(bs_out$HYPO$LRT_trecase) = c("LRT","DF","PVAL")
	rownames(bs_out$HYPO$LRT_trecase) = colnames(RHO)
	# bs_out$LRT_trecase

	# Hypothesis Testing (cis-trans)
	colnames(bs_out$HYPO$LRT_cistrans) = c("LRT","DF","PVAL")
	rownames(bs_out$HYPO$LRT_cistrans) = colnames(RHO)
	# bs_out$LRT_cistrans

	# Final Results
	res = smart_df(CT = colnames(RHO),
		LRT_trec = bs_out$HYPO$LRT_trec[,1],
		PVAL_trec = bs_out$HYPO$LRT_trec[,3],
		LRT_trecase = bs_out$HYPO$LRT_trecase[,1],
		PVAL_trecase = bs_out$HYPO$LRT_trecase[,3],
		LRT_cistrans = bs_out$HYPO$LRT_cistrans[,1],
		PVAL_cistrans = bs_out$HYPO$LRT_cistrans[,3])
	res$cistrans = ifelse(res$PVAL_cistrans < cistrans_thres,"trans","cis")
	res$PVAL_eQTL = as.numeric(apply(res[,c("PVAL_trec","PVAL_trecase","cistrans")],1,
		function(xx) ifelse(xx[3] == "cis",xx[2],xx[1])))
	rownames(res) = NULL
	bs_out$res = res
	rm(res)
	
	# Naming objects
	names(bs_out$OPT$iBETA) 	= colnames(XX)
	names(bs_out$OPT$PARS) = vname
	cell_types = colnames(RHO)
	if( bs_out$HYPO$swap_CT != 1 ){
		swap_CT = bs_out$HYPO$swap_CT
		cell_types[c(1,swap_CT)] = cell_types[c(swap_CT,1)]
		if( show ) message(sprintf("Cell type %s was set as the reference\n",colnames(RHO)[swap_CT]),appendLF = FALSE)
		names(bs_out$OPT$PARS)[kap_idx] = paste0("logKAP.",cell_types[-1])
		names(bs_out$OPT$PARS)[eta_idx] = paste0("logETA.",cell_types)
		names(bs_out$OPT$PARS)[alp_idx] = paste0("logALP.",cell_types)
	}
	# cell_types
	vname2 = names(bs_out$OPT$PARS)
	names(bs_out$OPT$GRAD) = vname2
	names(bs_out$OPT$upPARS) = vname2
	bs_out$OPT$HESS = smart_names(MAT = bs_out$OPT$HESS,
		ROW = vname2,COL = vname2)
	nz = which(diag(bs_out$OPT$HESS) != 0)
	np = nrow(bs_out$OPT$HESS)
	bs_out$OPT$COV = matrix(0,np,np)
	rcon = rcond(bs_out$OPT$HESS[nz,nz])
	if( rcon != 0 ){
		bs_out$OPT$COV[nz,nz] = solve(-bs_out$OPT$HESS[nz,nz],tol = 0.1 * rcon)
	} else {
		bs_out$OPT$COV[nz,nz] = NA
	}
	bs_out$HYPO$PVAL = smart_names(MAT = bs_out$HYPO$PVAL,
		ROW = colnames(RHO),COL = c("TReC","TReCASE","CisTrans"))
	bs_out$HYPO$mat_ETA = smart_names(MAT = bs_out$HYPO$mat_ETA,
		ROW = paste0("allele.",c("A","B")),COL = colnames(RHO))
	bs_out$OPT$eqtl_vars = smart_names(MAT = bs_out$OPT$eqtl_vars,
		ROW = c("KAPPA","ETA","ALPHA"),COL = cell_types)
	
	# Convergence status
	bs_out$OPT$Convergence_Status = "Success"
	if( bs_out$OPT$converge == 0 ){
		bs_out$OPT$Convergence_Status = "Failed"
	}
	if( show ){
		message(sprintf("Convergence_Status = %s\n",bs_out$OPT$Convergence_Status),appendLF = FALSE)
	}
	
	# Output
	bs_out
	
}

#' @title CSeQTL_GS
#' @description Main function that performs eQTL mapping on one 
#'	gene with its associated SNPs.
#' 
#' @param XX A numeric design matrix of baseline covariates 
#'	including the intercept in the first column and centered 
#'	continuous covariates. One of the non-intercept columns
#'	should correspond to centered log-transformed library size.
#'	The \code{rownames(XX)} needs to be specified.
#' @param TREC An integer vector containing one gene's 
#'	total read counts.
#'	The \code{names(TREC)} needs to be specified.
#' @param SNP A nonnegative integer matrix with genotypes. Values
#'	should be coded 0, 1, 2, 3, 5 for genotypes AA, AB, BA, BB, NA, respectively.
#'	Rows correspond to SNPs and columns correspond to subjects.
#'  The \code{rownames(SNP)} needs to be specified.
#' @param hap2 An integer vector of the second haplotype's
#'	counts for the gene. The \code{names(hap2)} needs to be specified
#'	and match \code{names(TREC)}.
#' @param ASREC An integer vector of the total haplotype
#'	counts (ASReC) for the gene. The \code{names(ASREC)} needs to be specified
#'	and match \code{names(TREC)}.
#' @param PHASE A binary 0/1 vector indicating if haplotype counts
#'	for the gene will be used. The \code{names(PHASE)} needs to be specified
#'	and match \code{names(TREC)}.
#' @param RHO A numeric matrix of cell type proportions. Rows 
#'	correspond to subjects and columns correspond to cell types. 
#'	The \code{rownames(RHO)} needs to be specified
#'	and match \code{names(TREC)}. The \code{colnames(RHO)} 
#'	also needs to be specified.
#' @param trim Boolean value set to \code{FALSE} by default to prevent
#'	outcome trimming. If \code{TRUE}, the CSeQTL model will be fitted
#'	without SNP genotype to calculate each subject's Cooks' distance
#'	for the gene.
#' @param thres_TRIM A positive numeric value to perform subject outcome trimming.
#'	Subjects with standardized Cooks' Distances greater than the threshold are trimmed.
#' @param numAS A positive integer to determine if a subject has
#'	enough total haplotype counts.
#' @param numASn A positive integer to determine how many subjects
#'	have at least \code{numAS} to use the haplotype counts.
#' @param numAS_het A positive integer to determine how many subjects
#'	with at least \code{numAS} are heterozygous (AB or BA). If
#'	\code{sum(PHASE == 1 & ASREC >= numAS & (SNP == 1 | SNP == 2)) >= numAS_het},
#'	those subjects haplotype counts will be used for TReCASE and cis/trans 
#'	testing and estimation.
#' @param cistrans A numeric value specifying the cis/trans test p-value cutoff
#'	to determine if the eQTLs from TReC-only or TReCASE model is reported.
#' @param ncores A positive integer specifying the number of threads available
#'	to decrease computational runtime when performing trimming and looping through SNPs.
#' @param show A boolean value to display verbose output.
#' 
#' @export
CSeQTL_GS = function(XX,TREC,SNP,hap2,ASREC,PHASE,RHO,trim = TRUE,
	thres_TRIM = 20,numAS = 5,numASn = 5,numAS_het = 5,cistrans = 0.01,
	ncores = 1,show = TRUE){
	
	# Check input classes are correct
	chk_XX(XX = XX)
	if( !is(TREC,"numeric") ) stop("TREC should be a integer/numeric vector!")
	if( !is(SNP,"matrix") ) 	stop("SNP should be a matrix!")
	if( !is(hap2,"numeric") ) 		stop("hap2 should be a integer/numeric vector!")
	if( !is(ASREC,"numeric") ) 	stop("ASREC should be a integer/numeric vector!")
	if( !is(PHASE,"numeric") ) 	stop("PHASE should be a integer/numeric vector!")
	if( !is(RHO,"matrix") ) 	stop("RHO should be a matrix!")
	if( !is(trim,"logical") ) 				stop("trim should be TRUE/FALSE!")
	
	if( trim ){
		if( !is(thres_TRIM,"numeric") )
			stop("thres_TRIM should be a positive numeric value!")
		if( thres_TRIM <= 0 )
			stop("thres_TRIM should be a positive numeric value!")
	}
	
	# Constants/Names
	NN = nrow(XX)
	QQ = ncol(RHO)
	ID = names(TREC)
	if( is.null(ID) ) 						stop("Specify names(TREC)!")
	if( is.null(names(hap2)) ) 		stop("Specify names(hap2)!")
	if( is.null(names(ASREC)) ) 	stop("Specify names(ASREC)!")
	if( is.null(names(PHASE)) ) 	stop("Specify names(PHASE)!")
	if( is.null(rownames(RHO)) ) 	stop("Specify rownames(RHO)!")
	if( is.null(colnames(RHO)) ) 	stop("Specify colnames(RHO)!")
	if( is.null(rownames(SNP)) ) 	stop("Specify rownames(SNP)!")
	if( is.null(colnames(SNP)) ) 	stop("Specify colnames(SNP)!")
	
	snp_ids 				= rownames(SNP)
	celltypes 			= colnames(RHO)
	SNP_subj_order 	= colnames(SNP)
	
	# Check sample size and subject ordering
	if( nrow(RHO) != NN || ncol(SNP) != NN || length(TREC) != NN 
		|| length(hap2) != NN || length(ASREC) != NN || length(PHASE) != NN )
		stop("Subject count mismatch!")
	if( !all(ID == names(hap2)) || !all(ID == names(ASREC))
		|| !all(ID == names(PHASE)) || !all(ID == rownames(RHO)) )
		stop("Subject order mismatch")
	if( !all(SNP_subj_order %in% ID) ) stop("colnames(SNP) issue!")
	
	# Make GS_index
	GS_index = matrix(NA,1,2)
	num_snps = nrow(SNP)
	GS_index[1,1] = 0
	GS_index[1,2] = num_snps - 1
	
	# Normalize RHO's rows to sum to 1
	RHO = RHO / rowSums(RHO)
	
	# Prep input formats
	TREC_0 	= matrix(TREC,nrow = 1)
	hap2 		= matrix(hap2,nrow = 1)
	ASREC_0 = matrix(ASREC,nrow = 1)
	PHASE_0 = matrix(PHASE,nrow = 1)
	
	# Run analysis
	start_time = Sys.time()
	gs_out = Rcpp_CSeQTL_GS(XX = XX,TREC_0 = TREC_0,SNP = SNP,hap2 = hap2,
		ASREC = ASREC_0,PHASE_0 = PHASE_0,RHO = RHO,GS_index = GS_index,trim = trim,
		swapCT = TRUE,trim_thres = thres_TRIM,numAS = numAS,numASn = numASn,
		numAS_het = numAS_het,cistrans = cistrans,max_iter = 4e3,eps = 1e-10,
		gr_eps = 1e-2,conv_eps = 5e-5,show = show,prompt = show,ncores = ncores)
	end_time = Sys.time()
	
	# Add row/column labels
	colnames(gs_out$iBETA) = colnames(XX)
	gs_out$ETA 					= smart_names(MAT = gs_out$ETA,ROW = snp_ids,COL = celltypes)
	gs_out$LRT_trec 		= smart_names(MAT = gs_out$LRT_trec,ROW = snp_ids,COL = celltypes)
	gs_out$LRT_trecase 	= smart_names(MAT = gs_out$LRT_trecase,ROW = snp_ids,COL = celltypes)
	gs_out$LRT_cistrans = smart_names(MAT = gs_out$LRT_cistrans,ROW = snp_ids,COL = celltypes)
	gs_out$MU_A 				= smart_names(MAT = gs_out$MU_A,ROW = snp_ids,COL = celltypes)
	gs_out$MU_B 				= smart_names(MAT = gs_out$MU_B,ROW = snp_ids,COL = celltypes)
	names(gs_out$LL) 		= snp_ids
	if( QQ == 1 ){
		pars_label = c(colnames(XX),"logPHI","logETA","logPSI","logALPHA")
	} else {
		pars_label = c(colnames(XX),"logPHI",paste0("logKAP",seq(2,QQ)),
			paste0("logETA",seq(QQ)),"logPSI",paste0("logALPHA",seq(QQ)))
	}
	gs_out$PARS 				= smart_names(MAT = gs_out$PARS,ROW = snp_ids,COL = pars_label)
	
	# Store run times
	gs_out$start_time = start_time
	gs_out$end_time = end_time
	
	# Adjust for multiple testing, calculate finalized pvalues
	PVAL_trec 			= 1 - pchisq(gs_out$LRT_trec[,,drop=FALSE],df = 1)
	PVAL_trecase 		= 1 - pchisq(gs_out$LRT_trecase[,,drop=FALSE],df = 1)
	PVAL_cistrans 	= 1 - pchisq(gs_out$LRT_cistrans[,,drop=FALSE],df = 1)
	CIS_eqtl 				= ifelse(PVAL_cistrans > cistrans,1,0)
	LRT_eqtl				= CIS_eqtl * gs_out$LRT_trecase + (1 - CIS_eqtl) * gs_out$LRT_trec
	PVAL_eqtl				= 1 - pchisq(LRT_eqtl,df = 1)
	
	ext_idx 				= apply(LRT_eqtl,2,which.max)
	MU_A 						= diag(gs_out$MU_A[ext_idx,,drop=FALSE])
	MU_B 						= diag(gs_out$MU_B[ext_idx,,drop=FALSE])
	ETA 						= apply(cbind(MU_A,MU_B),1,function(xx) ifelse(xx[1]==xx[2],1,xx[2]/xx[1]))
	ETA 						= as.numeric(ETA)
	
	out_res = smart_df(CELLTYPE = celltypes,SNP = snp_ids[ext_idx],
		MUa_trec = MU_A,ETA_trec = ETA,ETA_final = diag(gs_out$ETA[ext_idx,,drop=FALSE]),
		PVAL_final = diag(PVAL_eqtl[ext_idx,,drop=FALSE]),
		CISTRANS = ifelse(diag(CIS_eqtl[ext_idx,,drop=FALSE])==1,"CIS","TRANS"))
	print(out_res)
	
	gs_out$SNP_subj_order = SNP_subj_order
	
	# Output
	return(gs_out)
	
}

run_OLS = function(TREC,log_depth,XX,SNP,RHO,trim = FALSE,thres_TRIM = 20){
	if(FALSE){
		TREC = inputs$TREC
		log_depth = inputs$dat$log_depth
		XX = inputs$XX
		SNP = inputs$SNP
		RHO = inputs$RHO
		YY = NULL
		
		TREC = TREC; log_depth = dat$log_depth; XX = XX; 
		SNP = SNP; RHO = RHO; trim = trim; thres_TRIM = trim_thres
		
		
	}
	
	chk_XX(XX = XX)
	
	# Constants
	NN = nrow(RHO)
	QQ = ncol(RHO)
	SS = nrow(SNP)
	nX = ncol(XX)
	num_covars = nX + 1 + (QQ - 1) + (QQ - 1)
	# XX + SNP + RHO(-1) + SNP*RHO(-1)
	
	# Contrast matrix
	KK = matrix(0,QQ,num_covars)
	if( QQ > 1 ){
		colnames(KK) = c(colnames(XX),"SNP",colnames(RHO)[-1],
			paste0("SNP_CT.",colnames(RHO)[-1]))
	} else {
		colnames(KK) = c(colnames(XX),"SNP")
	}
	KK[,"SNP"] = 1 # eQTL 1st cell type
	if( QQ > 1 ){
		if( QQ == 2 ){
			KK[-1,grepl("SNP_CT",colnames(KK))] = 1
		} else {
			diag(KK[-1,grepl("SNP_CT",colnames(KK))]) = 1
		}
	}
	
	# SNP code as 0/1/2
	SNP2 = SNP
	SNP2[SNP2 == 2] = 1 	# Code heterozygous
	SNP2[SNP2 == 3] = 2 	# Code homozygous alternate
	SNP2[SNP2 == 5] = NA	# Code missing genotype
	
	# Inverse quantile normalized outcome after adjusting for library size
	YY = c(TREC) / exp(log_depth)
	YY = qnorm(p = rank(YY) / (NN + 1),mean = 0,sd = 1)
	YY_final = YY
	
	if( trim == TRUE ){
		if( QQ > 1 ){
			null_lm = lm(YY ~ .,data = smart_df(XX[,-1,drop=FALSE],RHO[,-1,drop=FALSE]))
		} else {
			null_lm = lm(YY ~ .,data = smart_df(XX[,-1,drop=FALSE]))
		}
		pred 		= predict(null_lm)
		cooksd 	= cooks.distance(null_lm)
		cooksd 	= (cooksd - median(cooksd)) / mad(cooksd)
		idx_high_cooksd = which(cooksd >= thres_TRIM)
		YY_final[idx_high_cooksd] = pred[idx_high_cooksd]
	}
	
	# Loop through snps for eqtl mapping
	PVAL 	= matrix(NA,SS,QQ); 
	PVAL 	= smart_names(PVAL,ROW = rownames(SNP2),COL = colnames(RHO))
	EST 	= PVAL
	for(ss in seq(SS)){
		# ss = 1
		# ss = which(out_res$SNP == rownames(SNP2))
		if( ss %% 5 == 0 ) message(".",appendLF = FALSE)
		if( ss %% 2e2 == 0 || ss == SS ) message(sprintf("%s out of %s\n",ss,SS),appendLF = FALSE)
		
		# Create interaction covariates between proportions and genotype
		if(QQ == 1){
			iXX = smart_df(SNP = SNP2[ss,])
		} else if(QQ >= 2){
			iXX = apply(RHO[,-1,drop=FALSE],2,function(xx) xx * SNP2[ss,])
			iXX = smart_df(iXX)
			names(iXX) = paste0("SNP_CT",seq(2,QQ))
			iXX = smart_df(SNP = SNP2[ss,],RHO[,-1,drop=FALSE],iXX); 
		}
		rownames(iXX) = NULL
		
		# Prepare final covariate dataframe
		notna_idx = as.integer(which(!is.na(SNP2[ss,])))
		fin_XX = smart_df(XX[,-1],iXX)[notna_idx,]; rm(iXX)
		YY_ana = YY_final[notna_idx]
		lm_full = lm(YY_ana ~ .,data = fin_XX)
		
		if(ss == 1){
			BETA = matrix(NA,SS,num_covars)
			BETA = smart_names(BETA,ROW = rownames(SNP2),
				COL = names(lm_full$coefficients))
		}
		BETA[ss,] = as.numeric(lm_full$coefficients)
		
		# Hypothesis testing
		tout = glht(model = lm_full,linfct = KK,alternative = "two.sided")
		# tout2 = summary(tout,test = adjusted("none"))$test$pvalues
		tout2 = pchisq(( summary(tout,test = adjusted("none"))$test$tstat )^2,
			df = 1,lower.tail = FALSE)
		out_df = smart_df(CT = seq(QQ),PVAL = tout2); rm(tout,tout2)
		
		if( any(out_df$PVAL == 0) ) stop("PVAL equals zero")
		
		# Get eqtl effect size estimates
		eqtl_est = c(KK %*% lm_full$coefficients)
		out_df$EST = eqtl_est
		
		PVAL[ss,] = out_df$PVAL
		EST[ss,]	= out_df$EST
		
		rm(out_df,eqtl_est,lm_full,fin_XX)
	}
	
	# Output
	list(EST = EST,BETA = BETA,PVAL = PVAL)
	
}



###

