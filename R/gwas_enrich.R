# ----------
# GWAS enrichment functions
# ----------
get_GWAS_catalog = function(work_dir){
	
	gwas_rds_fn = file.path(work_dir,"gwas_catalog_v1.0.rds")
	
	if( !file.exists(gwas_rds_fn) ){
	
		gwas_fn = file.path(work_dir,"gwas_catalog_v1.0.tsv")
		if( !file.exists(gwas_fn) ){
			setwd(work_dir)
			url_link = "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
			wget(url = url_link)
			gwas = fread("full",header = TRUE,sep = "\t",
				blank.lines.skip = TRUE,data.table = FALSE)
			dim(gwas); gwas[1:3,]
			fwrite(x = gwas,file = gwas_fn,sep = "\t",
				col.names = TRUE,row.names = FALSE)
			
		}
		
		message(sprintf("%s: Import GWAS catalog ...\n",date()),appendLF = FALSE)
		gwas = fread(gwas_fn,header = TRUE,sep = "\t",
			blank.lines.skip = TRUE,data.table = FALSE)
		
		# Polish gwas
		gwas = name_change(gwas,"DISEASE/TRAIT","PHENO")
		gwas$PHENO = toupper(gwas$PHENO)
		gwas = gwas[which(gwas$CHR_ID != ""),]
		gwas = gwas[which(gwas$CHR_ID %in% as.character(seq(22))),]
		tab = table(gwas$CHR_ID,useNA = "ifany"); tab
		gwas = smart_rmcols(gwas,c("CNV","MERGED","CONTEXT","PLATFORM [SNPS PASSING QC]",
			"STUDY","DATE ADDED TO CATALOG","PUBMEDID","FIRST AUTHOR","LINK","INITIAL SAMPLE SIZE",
			"REPLICATION SAMPLE SIZE"))
		gwas$CHR_POS = as.integer(gwas$CHR_POS)
		
		# Prep Phenotype groupings
		message(sprintf("%s: Make phenotype groupings ...\n",date()),appendLF = FALSE)
		if( TRUE ){
		gwas$myPHENO = NA
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("INSOMNIA",gwas$PHENO) | grepl("ALZHEIMER",gwas$PHENO)
			| grepl("CHRONOTYPE",gwas$PHENO) | grepl("SLEEP",gwas$PHENO)
			| grepl("ADHD",gwas$PHENO) | grepl("MIGRAINE",gwas$PHENO)
			| grepl("SCHIZOPHRENIA",gwas$PHENO) | grepl("SLEEP",gwas$PHENO)
			| grepl("SCLEROSIS",gwas$PHENO) | grepl("FLUID INTELL",gwas$PHENO)
			| grepl("NEUROTIC",gwas$PHENO) | grepl("PSYCHIATRIC",gwas$PHENO)
			| grepl("PYSCHOLOGIC",gwas$PHENO) | grepl("PARKINSON",gwas$PHENO)
			| grepl("DEPRESS",gwas$PHENO) | grepl("EPILEPSY",gwas$PHENO)
			| grepl("BRAIN",gwas$PHENO) | grepl("COGNITIVE",gwas$PHENO)
			| grepl("REACTION TIME",gwas$PHENO) | grepl("WHITE MATTER MICROSTRUCTURE",gwas$PHENO) 
			| grepl("BIPOLAR",gwas$PHENO) | grepl("GENERAL RISK TOLERANCE",gwas$PHENO) 
			| grepl("CORTICAL SURFACE AREA",gwas$PHENO) | grepl("ISCHAEMIC STROKE",gwas$PHENO) 
			| grepl("POSITIVE AFFECT",gwas$PHENO) 
			| grepl("NEUROCITICISM",gwas$PHENO) | grepl("ATTENTION DEFICIT HYPERACTIVITY",gwas$PHENO) 
			| grepl("AMPHETAMINES",gwas$PHENO) | grepl("INTRACRANIAL ANEURYSM",gwas$PHENO)
			| grepl("ISCHEMIC STROKE",gwas$PHENO) | grepl("CEREBROSPINAL FLUID",gwas$PHENO)
			| grepl("MOOD SWINGS",gwas$PHENO) | grepl("GLIOMA",gwas$PHENO)
			| grepl("HIPPOCAMPAL VOLUME",gwas$PHENO) | grepl("FEELING NERVOUS",gwas$PHENO)
			| grepl("IRRITABLE MOOD",gwas$PHENO) | grepl("CEREBROSPINAL FLUID",gwas$PHENO)
			| grepl("FEELING MISERABLE",gwas$PHENO) | grepl("FEELING FED-UP",gwas$PHENO)
			| grepl("CEREBRAL AMYLOID DEPOSITION",gwas$PHENO) | grepl("RESPONSE TO ANTIPSYCHOTIC",gwas$PHENO)
			| grepl("LONELINESS",gwas$PHENO) | grepl("NEUROFIBRILLARY TANGLES",gwas$PHENO)
			))] = "Psychiatric_neurologic"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("EGG BW",gwas$PHENO) | grepl("INTRACRANEAL",gwas$PHENO)
			| grepl("HEIGHT",gwas$PHENO) | grepl("BMI",gwas$PHENO)
			| grepl("BODY FAT",gwas$PHENO) | grepl("FOREARM",gwas$PHENO)
			| grepl("BODY MASS INDEX",gwas$PHENO) | grepl("WAIST-HIP RATIO",gwas$PHENO)
			| grepl("WAIST CIRCUMFERENCE",gwas$PHENO) | grepl("BIRTH WEIGHT",gwas$PHENO)
			| grepl("HIP CIRCUMFERENCE",gwas$PHENO) | grepl("WEIGHT",gwas$PHENO)
			| grepl("LEAN BODY MASS",gwas$PHENO)
			))] = "Anthropometric"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("CARDIO",gwas$PHENO) | grepl("GLUCOSE",gwas$PHENO)
			| grepl("HYPERTENSION",gwas$PHENO) | grepl("CHOLESTEROL",gwas$PHENO)
			| grepl("VASCULAR",gwas$PHENO) | grepl("HDL",gwas$PHENO)
			| grepl("LDL",gwas$PHENO) | grepl("FAT",gwas$PHENO)
			| grepl("INSULIN",gwas$PHENO) | grepl("TRIGLYCERIDE",gwas$PHENO)
			| grepl("OBESITY-RELATED",gwas$PHENO) | grepl("CORONARY ARTERY",gwas$PHENO)
			| grepl("ABDOMINAL AORTIC",gwas$PHENO) | grepl("PR INTERVAL",gwas$PHENO)
			| grepl("ATRIAL FIBRILLATION",gwas$PHENO) | grepl("QT INTERVAL",gwas$PHENO)
			| grepl("MEAN ARTERIAL PRESSURE",gwas$PHENO) | grepl("QRS DURATION",gwas$PHENO)
			| grepl("HEART FAILURE WITH REDUCED EJECTION FRACTION",gwas$PHENO)
			| grepl("VISCERAL ADIPOSE",gwas$PHENO) | grepl("RESTING HEART RATE",gwas$PHENO)
			| grepl("OBESITY",gwas$PHENO) | grepl("CORONARY HEART",gwas$PHENO)
			| grepl("ADIPONECTIN LEVELS",gwas$PHENO) | grepl("METABOLITE LEVELS",gwas$PHENO)
			| grepl("METABOLIC SYNDROME",gwas$PHENO) | grepl("TAKAYASU ARTERITIS",gwas$PHENO)
			| grepl("SUBCUTANEOUS ADIPOSE",gwas$PHENO) | grepl("MYOCARDIAL INFARCTION",gwas$PHENO)
			))] = "Cardiometabolic"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("EOSINOPHIL",gwas$PHENO) | grepl("GRANULOCYTE",gwas$PHENO)
			| grepl("RETICULOCYTE",gwas$PHENO) | grepl("LYMPHOCYTE",gwas$PHENO)
			| grepl("MONOCYTE",gwas$PHENO) | grepl("MYELOID",gwas$PHENO)
			| grepl("NEUTROPHIL",gwas$PHENO) | grepl("PLATELET",gwas$PHENO)
			| grepl("RED BLOOD",gwas$PHENO) | grepl("BASOPHIL",gwas$PHENO)
			| grepl("WHITE BLOOD",gwas$PHENO) | grepl("BLOOD PROTEIN",gwas$PHENO)
			| grepl("SYSTOLIC BLOOD",gwas$PHENO) | grepl("PULSE PRESSURE",gwas$PHENO)
			| grepl("BLOOD PRESSURE",gwas$PHENO) | grepl("CORPUSCULAR HEMOGLOBIN",gwas$PHENO)
			| grepl("RED CELL DISTRIBUTION",gwas$PHENO) | grepl("C-REACTIVE PROTEIN",gwas$PHENO)
			| grepl("SERUM METABOLITE",gwas$PHENO) | grepl("MEAN CORPUSCULAR",gwas$PHENO)
			| grepl("APOLIPOPROTEIN A1",gwas$PHENO) | grepl("HEMOGLOBIN",gwas$PHENO)
			| grepl("BLOOD METABOLITE",gwas$PHENO) | grepl("APOLIPOPROTEIN B",gwas$PHENO)
			| grepl("SERUM URIC ACID",gwas$PHENO) | grepl("HEMATOCRIT",gwas$PHENO)
			| grepl("BLOOD UREA NITROGEN",gwas$PHENO) | grepl("CHRONIC LYMPHOCYTIC LEUKEMIA",gwas$PHENO)
			| grepl("CREATININE LEVEL",gwas$PHENO) | grepl("MACROPHAGE INFLAMMATORY PROTEIN",gwas$PHENO)
			| grepl("VENOUS THROMBOEMBOLISM",gwas$PHENO) | grepl("VWF LEVELS",gwas$PHENO)
			| grepl("PLASMA FREE AMINO ACID",gwas$PHENO) | grepl("ACUTE LYMPHOBLASTIC LEUKEMIA",gwas$PHENO)
			| grepl("SERUM 25-HYDROXYVITAMIN D",gwas$PHENO) | grepl("LEUKOCYTE TELOMERE",gwas$PHENO)
			| grepl("LIPOPROTEIN",gwas$PHENO) | grepl("BILIRUBIN",gwas$PHENO)
			| grepl("FACTOR VIII",gwas$PHENO) | grepl("MULTIPLE MYELOMA",gwas$PHENO)
			| grepl("POLYCHLORINATED BIPHENYL",gwas$PHENO)
			))] = "Blood"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("CROHNS",gwas$PHENO) | grepl("BOWEL DISEASE",gwas$PHENO)
			| grepl("ULCERATIVE",gwas$PHENO) | grepl("IMMUNO",gwas$PHENO)
			| grepl("OKADA",gwas$PHENO) | grepl("DIABETES",gwas$PHENO)
			| grepl("ANKYLOSING",gwas$PHENO) | grepl("PSORIASIS",gwas$PHENO)
			| grepl("RHEUMATOID",gwas$PHENO) | grepl("DOCTOR ASTHMA",gwas$PHENO)
			| grepl("ALLERGIC RHIN",gwas$PHENO) | grepl("SYSTEMIC LUPUS",gwas$PHENO)
			| grepl("CROHN'S DISEASE",gwas$PHENO) | grepl("CELIAC DISEASE",gwas$PHENO)
			| grepl("AUTOIMMUNE",gwas$PHENO) | grepl("REPORTED ASTHMA",gwas$PHENO)
			))] = "Immune"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("AGE AT DEATH",gwas$PHENO) | grepl("MALE-PATTERN BALDNESS",gwas$PHENO)
			| grepl("IGG GLYCOSYLATION",gwas$PHENO)
			| grepl("AGE EFFECT",gwas$PHENO) | grepl("LONGEVITY",gwas$PHENO)
			))] = "Aging"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("IRRITABLE BOWEL",gwas$PHENO) | grepl("GUT MICROBIOTA",gwas$PHENO)
			| grepl("DIVERTICULAR DISEASE",gwas$PHENO)
			))] = "Digestion"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("HYPERTHYROID",gwas$PHENO) | grepl("HYPOTHYROID",gwas$PHENO)
			| grepl("THYROID",gwas$PHENO)
			))] = "Endocrine"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("OSTEOPOROSIS",gwas$PHENO) | grepl("GOUT",gwas$PHENO)
			| grepl("HEEL BONE",gwas$PHENO) | grepl("IDIOPATHIC SCOLIOSIS",gwas$PHENO)
			| grepl("APPENDICULAR LEAN MASS",gwas$PHENO) | grepl("BONE MINERAL",gwas$PHENO)
			| grepl("FEMORAL NECK BONE",gwas$PHENO) | grepl("OSTEOARTHRITIS",gwas$PHENO)
			| grepl("PSORIATIC ARTHRITIS",gwas$PHENO)
			))] = "Skeletal"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("HAIR",gwas$PHENO) 
			| grepl("BALDING",gwas$PHENO)
			))] = "Morphology"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("EDUCATIONAL ATTAINMENT",gwas$PHENO) | grepl("HIGHEST MATH CLASS",gwas$PHENO)
			| grepl("SELF-REPORTED MATH",gwas$PHENO) | grepl("INTELLIGENCE",gwas$PHENO)
			| grepl("HOUSEHOLD INCOME",gwas$PHENO)
			))] = "Education_wealth"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("LUNG FUNCTION",gwas$PHENO) | grepl("BRONCHODILATOR",gwas$PHENO)
			| grepl("ASTHMA",gwas$PHENO) | grepl("SMOKING",gwas$PHENO)
			| grepl("SMOKERS",gwas$PHENO) | grepl("FEV1",gwas$PHENO)
			| grepl("CHRONIC OBSTRUCTIVE PULMONARY",gwas$PHENO)
			| grepl("PEAK EXPIRATORY FLOW",gwas$PHENO) | grepl("TUBERCULOSIS",gwas$PHENO)
			| grepl("LUNG ADENOCARCINOMA",gwas$PHENO) | grepl("RESPIRATORY DISEASES",gwas$PHENO)
			| grepl("SQUAMOUS CELL LUNG",gwas$PHENO) | grepl("CIGARETTES SMOKED",gwas$PHENO)
			| grepl("VELOPHARYNGEAL DYSFUNCTION",gwas$PHENO) | grepl("LUNG CANCER",gwas$PHENO)
			| grepl("PULMONARY FUNCTION",gwas$PHENO) | grepl("DIFFUSING CAPACITY OF CARBON MONOXIDE",gwas$PHENO)
			| grepl("PNEUMONIA",gwas$PHENO)
			))] = "Respiratory"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("GLOMERULAR FILTRATION",gwas$PHENO)
			| grepl("URATE LEVELS",gwas$PHENO) | grepl("URINARY",gwas$PHENO)
			| grepl("CHRONIC KIDNEY",gwas$PHENO) | grepl("DIABETIC KIDNEY",gwas$PHENO)
			| grepl("IGA NEPHROPATHY",gwas$PHENO)
			))] = "Kidney_urinary"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("INTRAOCULAR PRESSURE",gwas$PHENO)
			| grepl("REFRACTIVE ERROR",gwas$PHENO) | grepl("SPHERICAL EQUIVALENT",gwas$PHENO)
			| grepl("CENTRAL CORNEAL THICKNESS",gwas$PHENO) | grepl("MACULAR THICKNESS",gwas$PHENO)
			| grepl("GLAUCOMA",gwas$PHENO) | grepl("OPTIC DISC SIZE",gwas$PHENO)
			| grepl("OPTIC CUP AREA",gwas$PHENO) | grepl("VERTICAL CUP-DISC RATIO",gwas$PHENO)
			| grepl("MYOPIA",gwas$PHENO) | grepl("AGE-RELATED MACULAR",gwas$PHENO)
			| grepl("CORNEAL ASTIGMATISM",gwas$PHENO) | grepl("OPTIC DISC AREA",gwas$PHENO)
			| grepl("OCULAR SARCOIDOSIS",gwas$PHENO) | grepl("REFRACTIVE ASTIGMATISM",gwas$PHENO)
			))] = "Eyes"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("TONSILLECTOMY",gwas$PHENO) | grepl("GINGIVAL CREVICULAR",gwas$PHENO)
			| grepl("EROSIVE TOOTH",gwas$PHENO) | grepl("MOUTH ULCERS",gwas$PHENO)
			| grepl("OROFACIAL CLEFTS",gwas$PHENO) | grepl("DENTAL CARIES",gwas$PHENO)
			| grepl("CLEFT PALATE",gwas$PHENO) | grepl("NONSYNDROMIC CLEFT LIP",gwas$PHENO)
			))] = "Mouth_teeth"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("MOSQUITO BITE",gwas$PHENO) | grepl("BASAL CELL",gwas$PHENO)
			| grepl("ITCH INTENSITY",gwas$PHENO) | grepl("INFLAMMATORY SKIN",gwas$PHENO)
			| grepl("SUNBURNS",gwas$PHENO) | grepl("LOW TAN RESPONSE",gwas$PHENO)
			| grepl("DERMATITIS",gwas$PHENO) | grepl("ECZEMA",gwas$PHENO)
			| grepl("SKIN PIGMENTATION",gwas$PHENO) | grepl("SHINGLES",gwas$PHENO)
			| grepl("MELANOMA",gwas$PHENO) | grepl("PLANTAR WARTS",gwas$PHENO)
			))] = "Skin"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("COLONOSCOPY",gwas$PHENO) | grepl("COLORECTAL",gwas$PHENO)
			))] = "Colon"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("PRIMARY BILIARY",gwas$PHENO) | grepl("GAMMA GLUTAMYL TRANSFERASE",gwas$PHENO)
			))] = "Liver"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("ALCOHOL",gwas$PHENO) | grepl("NICOTINE DEPENDENCE",gwas$PHENO)
			| grepl("COFFEE CONSUMPTION",gwas$PHENO)
			))] = "substance_addiction"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("ENDOMETRIOSIS",gwas$PHENO) | grepl("UTERINE FIBROIDS",gwas$PHENO)
			| grepl("MENARCHE",gwas$PHENO) | grepl("MENOPAUSE",gwas$PHENO)
			| grepl("OVARIAN",gwas$PHENO) | grepl("CERVICAL",gwas$PHENO)
			| grepl("ENDOMETRIAL",gwas$PHENO)
			))] = "Female_reproduction"
		gwas$myPHENO[which(is.na(gwas$myPHENO) & 
			( grepl("EASE OF GETTING UP",gwas$PHENO) | grepl("MORNING PERSON",gwas$PHENO)
			| grepl("MORNINGNESS",gwas$PHENO)
			))] = "Circadian"
		}
		
		saveRDS(gwas,gwas_rds_fn)
		
	}
	gwas = readRDS(gwas_rds_fn)
	
	return(gwas)
}
run_gwasEnrich_analysis = function(DATA,work_dir,
	which_gwas = "gwas_catalog",nBLOCKS = 200,
	verbose = TRUE){
	
	if(FALSE){
		work_dir = "C:/Users/Admin/Downloads"
		which_gwas = "gwas_catalog"
		nBLOCKS = 200
		
		DATA = dat
		
	}
	
	# Check DATA format
	req_nms = c("Chr","POS")
	if( !all(req_nms %in% names(DATA)) ){
		miss_nms = req_nms[!(req_nms %in% names(DATA))]
		stop(sprintf("Missing column(s): %s",paste(miss_nms,collapse = ",")))
	}
	DATA$POS = as.integer(DATA$POS)
	
	all_EE = names(DATA)
	all_EE = all_EE[grepl("^EE_",all_EE)]
	message(sprintf("%s: Detected %s columns with 'EE_' ...\n",date(),length(all_EE)),appendLF = FALSE)
	
	# Get GWAS catalog
	gwas = get_GWAS_catalog(work_dir = work_dir)
	
	# Subset phenotype
	sub_phenotypes = names(table(gwas$myPHENO))
	# print(sub_phenotypes)
	if( which_gwas == "gwas_catalog" ){
		gwas = gwas
	} else if( which_gwas %in% sub_phenotypes ){
		gwas = gwas[which(gwas$myPHENO == which_gwas),]
	} else {
		stop("No code for that phenotype!")
	}
	# dim(gwas); gwas[1:3,]
	
	# Calculate if locus is GWAS hit
	chrs = paste0("chr",1:22)
	DATA$GG = 0
	for(chr in chrs){
		# chr = chrs[1]; chr
		if( verbose ) message(sprintf("%s: Get GG for chr = %s ...\n",date(),chr),appendLF = FALSE)
		gwas_chr = gwas[which(gwas$CHR_ID == gsub("chr","",chr)),]
		idx = DATA$Chr == chr
		eqtl_pos = DATA$POS[idx]
		# length(eqtl_pos); length(unique(eqtl_pos))
		# there are duplicate positions b/c genes may be close together and position was tested multiple times
		u_gwas_pos = intersect(unique(gwas_chr$CHR_POS),unique(eqtl_pos)); length(u_gwas_pos)
		DATA$GG[idx] = ifelse(DATA$POS[idx] %in% u_gwas_pos,1,0)
		rm(idx,u_gwas_pos,gwas_chr,eqtl_pos)
	}
	
	# Calculate overall gwas enrichment
	message(sprintf("%s: Calculate overall enrichment ...\n",date()),appendLF = FALSE)
	num_test_loci = nrow(DATA)
	num_gwas = sum(DATA$GG)
	gwas_enrich = c()
	for(one_EE in all_EE){
		GROUP = gsub("EE_","",one_EE)
		one_TT = sprintf("TT_%s",GROUP)
		num_eqtl = sum(DATA[,one_EE])
		num_eqtl_gwas = sum(DATA[,one_EE] * DATA$GG)
		tmp_df = smart_df(GROUP = GROUP,
			n_test = num_test_loci,n_eqtl = num_eqtl,
			n_gwas = num_gwas,n_gwas_eqtl = num_eqtl_gwas)
		gwas_enrich = rbind(gwas_enrich,tmp_df)
		rm(tmp_df)
	}
	gwas_enrich$frac_gwas = gwas_enrich$n_gwas_eqtl / gwas_enrich$n_gwas
	gwas_enrich$frac_all = gwas_enrich$n_eqtl / gwas_enrich$n_test
	gwas_enrich$enrich = gwas_enrich$frac_gwas / gwas_enrich$frac_all
	gwas_enrich$enrich[is.na(gwas_enrich$enrich) | gwas_enrich$enrich == 0] = 1
	if( verbose ) print(gwas_enrich)
	
	# Sort positions and run block jackknife for enrichment
	message(sprintf("%s: Run block jackknife enrichment ...\n",date()),appendLF = FALSE)
	chrs2 = intersect(chrs,unique(DATA$Chr))
	chrs2 = chrs[chrs %in% chrs2]; chrs2
	DATA$chr2 = factor(DATA$Chr,levels = chrs)
	DATA = DATA[order(DATA$chr2,DATA$POS),]
	BLOCK_size = floor(num_test_loci / nBLOCKS) + 1; BLOCK_size
	DATA$BLOCK = sort(rep(seq(nBLOCKS),BLOCK_size))[1:num_test_loci]
	uBLOCKs = unique(DATA$BLOCK); uBLOCKs
	gwas_enrich_block = c()
	for(BLOCK in uBLOCKs){
		
		if( verbose ){
			if( BLOCK %% 2 == 0 ) message(".",appendLF = FALSE)
			if( BLOCK %% 30 == 0 || BLOCK == nBLOCKS )
				message(sprintf("%s out of %s\n",BLOCK,nBLOCKS),appendLF = FALSE)
		}
		
		DATA_block = DATA[which(DATA$BLOCK != BLOCK),] # remove one block
		num_test_loci = nrow(DATA_block)
		num_gwas = sum(DATA_block$GG)
		
	for(one_EE in all_EE){
		GROUP = gsub("EE_","",one_EE)
		idx = which(gwas_enrich$GROUP == GROUP)
		num_eqtl = sum(DATA_block[,one_EE])
		if( gwas_enrich$n_gwas_eqtl[idx] > 0 ){
			num_eqtl_gwas = sum(DATA_block[,one_EE] * DATA_block$GG)
			tmp_df = smart_df(BLOCK = BLOCK,GROUP = GROUP,
				n_test = num_test_loci,n_eqtl = num_eqtl,
				n_gwas = num_gwas,n_gwas_eqtl = num_eqtl_gwas)
		} else {
			tmp_df = smart_df(BLOCK = BLOCK,GROUP = GROUP,
				n_test = num_test_loci,n_eqtl = num_eqtl,
				n_gwas = num_gwas,n_gwas_eqtl = 0)
		}
		gwas_enrich_block = rbind(gwas_enrich_block,tmp_df)
	}
		rm(DATA_block)
	}
	gwas_enrich_block$frac_gwas = gwas_enrich_block$n_gwas_eqtl / gwas_enrich_block$n_gwas
	gwas_enrich_block$frac_all = gwas_enrich_block$n_eqtl / gwas_enrich_block$n_test
	gwas_enrich_block$enrich = gwas_enrich_block$frac_gwas / gwas_enrich_block$frac_all
	gwas_enrich_block$enrich[is.na(gwas_enrich_block$enrich) | gwas_enrich_block$enrich == 0] = 1
	
	# Get summary statistics
	gwas_enrich$log_enrich = log(gwas_enrich$enrich)
	gwas_enrich$log_enrich_meanJK = NA
	gwas_enrich$log_enrich_lowJK = NA
	gwas_enrich$log_enrich_highJK = NA
	for(one_EE in all_EE){
		# one_EE = all_EE[1]
		GROUP = gsub("EE_","",one_EE)
		idx1 = which(gwas_enrich$GROUP == GROUP)
		idx2 = which(gwas_enrich_block$GROUP == GROUP)
		out_JK = calc_JK(EST = log(gwas_enrich$enrich[idx1]),
			LOO_EST = log(gwas_enrich_block$enrich[idx2]),
			alpha = 0.05)
		gwas_enrich$log_enrich_meanJK[idx1] = out_JK$JK_mean
		gwas_enrich$log_enrich_lowJK[idx1] = out_JK$JK_CI[1]
		gwas_enrich$log_enrich_highJK[idx1] = out_JK$JK_CI[2]
		rm(out_JK,idx1,idx2)
	}
	gwas_enrich = smart_df(wGWAS = which_gwas,gwas_enrich)
	
	return(gwas_enrich)
	
}

# Simulate an input dataset and run GWAS enrichment analysis
sim_EQTL_RES = function(work_dir){ 
	if(FALSE){
		
		work_dir = "C:/Users/Admin/Downloads"
		gwas = get_GWAS_catalog(work_dir = work_dir)
		gwas = gwas[which(gwas$CHR_ID == "1"),]
		# dim(gwas); gwas[1:5,]
		tab = table(gwas$myPHENO); tab
		
	}
	
	num_loci = 4e4
	dat = smart_df(Chr = rep("chr1",num_loci))
	prob = c(1,1); prob = prob / sum(prob)
	dat$CLASS = sample(c(0,1),num_loci,
		replace = TRUE,prob = prob)
	dat$POS = NA
	idx = which(dat$CLASS == 0); length(idx)
	max_pos = max(as.integer(gwas$CHR_POS))
	dat$POS[idx] = sample(max_pos,
		length(idx),replace = FALSE)
	idx = which(dat$CLASS == 1)
	dat$POS[idx] = sample(as.integer(gwas$CHR_POS),
		length(idx),replace = FALSE)
	dat[1:10,]
	
	prob = c(5,1); prob = prob / sum(prob)
	dat$EE_A1 = sample(c(0,1),num_loci,
		replace = TRUE,prob = prob)
	prob = c(5,1); prob = prob / sum(prob)
	dat$EE_A2 = sample(c(0,1),num_loci,
		replace = TRUE,prob = prob)
	
	# Induce enrichment
	message(sprintf("%s: Induce enrichment in simulated dataset ...\n",date()),appendLF = FALSE)
	idx = which(dat$CLASS == 1)
	dat$EE_A2[idx] = sapply(dat$EE_A2[idx],function(xx){
		ifelse(xx == 1,1,rbinom(1,1,0.1))
	},USE.NAMES = FALSE)
	
	table(dat$EE_A1)
	table(dat$EE_A2)
	table(dat$EE_A1,dat$EE_A2)
	# dat[1:5,]
	
	# Run enrichment analyses: by phenotype group and across all phenotypes
	res = run_gwasEnrich_analysis(DATA = dat,work_dir = work_dir,
		which_gwas = "gwas_catalog",nBLOCKS = 200,verbose = FALSE)
	for(wGWAS in names(tab)){
		message(sprintf("%s: wGWAS = %s ...\n",date(),wGWAS),appendLF = FALSE)
		res = rbind(res,run_gwasEnrich_analysis(DATA = dat,
			work_dir = work_dir,which_gwas = wGWAS,
			nBLOCKS = 200,verbose = FALSE))
	}
	
	# Make second group of calculations, aka shuffling some rows
	res2 = res
	idx = sample(nrow(res2))
	res2[,c("wGWAS","GROUP")] = res[idx,c("wGWAS","GROUP")]
	fres = rbind(smart_df(res,MODEL = "CSeQTL"),
		smart_df(res2,MODEL = "OLS"))
	rm(res2)
	
	# Plot
	pd = position_dodge(width = 0.75) # control point spread
	themes = theme(text = element_text(size = 24),
		# axis.text.x = element_text(size = 12),
		# axis.title.y = element_text(size = 20,face = "bold"),
		panel.background = element_blank(),
		panel.spacing = unit(0.5,"lines"),
		panel.border = element_rect(color = "black",fill = NA,size = 1),
		legend.position = c("none","bottom")[2],
		legend.text = element_text(size = 26))
	
	tmp_range = max(abs(fres$log_enrich_meanJK)); tmp_range
	tmp_range = 1.2 * c(-1,1) * ifelse(tmp_range >= 1,tmp_range,1); tmp_range
	fres$tmp_col = ifelse(fres$log_enrich_lowJK > 0,"Significant","Non-significant")
	fres$wGWAS = as.character(fres$wGWAS)
	lev_wGWAS = sort(unique(fres$wGWAS)); # lev_wGWAS
	lev_wGWAS = c("gwas_catalog",lev_wGWAS[lev_wGWAS != "gwas_catalog"])
	lev_wGWAS = rev(lev_wGWAS)
	fres$wGWAS = factor(fres$wGWAS,levels = lev_wGWAS)
	
	log_enrich_meanJK = NULL
	log_enrich_lowJK = NULL
	log_enrich_highJK = NULL
	MODEL = tmp_col = log_enrich = NULL
	
	gg = ggplot(data = fres,aes(x = wGWAS,y = log_enrich_meanJK,
			ymin = log_enrich_lowJK,ymax = log_enrich_highJK,group = MODEL)) +
		# make sure to include 'group' option in above line, otherwise arrangement gets weird
		geom_errorbar(position = pd,size = 1,width = 0.1,aes(color = MODEL)) +
		geom_point(position = pd,size = 5,
			aes(shape = tmp_col,color = MODEL,stroke = c(2,1.5)[2])) + 
		geom_point(data = fres,position = pd,
			mapping = aes(x = wGWAS,y = log_enrich,group = MODEL)) +
		facet_grid(~ GROUP) + labs(shape = "Inference",color = "Cohort") +
		scale_shape_manual(values = c(1,19)) +
		geom_hline(yintercept = 0,linetype = 2) +
		ylab("log(Enrichment)") + xlab("") +
		coord_flip(ylim = tmp_range) +
		themes +
		guides(color = guide_legend(override.aes = list(size = 5)),
			shape = guide_legend(override.aes = list(size = 5)))
	png_fn = file.path(work_dir,"enrich.png")
	ggsave(filename = png_fn,plot = gg,width = 15,height = 13,units = "in")
	
}


###

