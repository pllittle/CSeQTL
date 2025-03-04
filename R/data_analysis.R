# ----------
# Real Data Analysis functions
# ----------
#' @title prep_gene_info
#' @description Prepare GTF files for gene and 
#'	exon level information
#' @param work_dir A character string specifying
#'	the working directory to store gene and exon
#'	data.
#' @param gtf_gz_fn A character string full path
#'	of the gz compressed gtf file
#' @return Null. No value returned.
#' @export
prep_gene_info = function(work_dir,gtf_gz_fn = NULL){
	
	if(FALSE){
		work_dir = "C:/Users/Admin/Downloads"
		gtf_gz_fn = file.path(work_dir,
			"gencode.v26.GRCh38.ERCC.genes.gtf.gz")
		# library(GenomicFeatures)
		
	}
	
	smart_mkdir(work_dir)
	
	# Get gtf
	if( is.null(gtf_gz_fn) ){
		tmp_url = "https://personal.broadinstitute.org"
		tmp_url = file.path(tmp_url,"francois/topmed")
		tmp_url = file.path(tmp_url,"gencode.v26.GRCh38.ERCC.genes.gtf.gz")
		message(sprintf("A possible input file can be obtained at \n  %s\n",tmp_url),appendLF = FALSE)
		
	}
	
	# Prep input for TREC/ASREC script
	exon_fn = file.path(work_dir,"exon_by_genes.rds")
	if( !file.exists(exon_fn) ){
		message("Make TxDb...\n",appendLF = FALSE)
		exdb = suppressWarnings(makeTxDbFromGFF(file = gtf_gz_fn,format = "gtf"))
		exons_list_per_gene = exonsBy(x = exdb,by = "gene")
		saveRDS(object = exons_list_per_gene,file = exon_fn)
	}
	
	# Get gene name mapping, lengths, start/end gene body for collecting SNPs and deconvolution
	gene_fn = file.path(work_dir,"gene_info.tsv.gz")
	if( !file.exists(gene_fn) ){
		gtf = fread(input = gtf_gz_fn,header = FALSE,
			sep = "\t",data.table = FALSE)
		gtf = gtf[which(gtf$V3 == "gene"),]
		gtf = name_change(gtf,"V1","Chr")
		gtf = name_change(gtf,"V4","Start")
		gtf = name_change(gtf,"V5","End")
		gtf = name_change(gtf,"V2","SOURCE")
		gtf = name_change(gtf,"V7","Strand")
		gtf$ENSG = sapply(gtf$V9,function(xx) strsplit(gsub("\\\"","",
			strsplit(xx,";")[[1]][1])," ")[[1]][2],USE.NAMES = FALSE)
		gtf$GENE = sapply(gtf$V9,function(xx) strsplit(gsub("\\\"","",
			strsplit(xx,";")[[1]][4])," ")[[1]][3],USE.NAMES = FALSE)
		gtf = smart_rmcols(gtf,c("V3","V6","V8","V9"))
		rownames(gtf) = NULL
		head(gtf)
		
		gene_fn2 = gsub(".gz$","",gene_fn)
		fwrite(x = gtf,file = gene_fn2)
		gzip(gene_fn2)
		
	}
	
	

}

###