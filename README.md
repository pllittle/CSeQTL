<div align="left">
<a href=""><img src="https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink&label=CSeQTL" height="80" /></a>
</div>

<!-- badges: start -->
![C++](https://img.shields.io/badge/C++-%2300599C.svg?style=square&logo=c%2B%2B&logoColor=gold)
![R](https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink)
![CRAN status](https://www.r-pkg.org/badges/version/CSeQTL)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2022.03.31.486605.svg)](https://doi.org/10.1101/2022.03.31.486605)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/github/languages/code-size/pllittle/CSeQTL.svg)](https://github.com/pllittle/CSeQTL)
[![](https://img.shields.io/github/last-commit/pllittle/CSeQTL.svg)](https://github.com/pllittle/CSeQTL/commits/master)
<!-- badges: end -->

## Introduction

Expression quantitative trait loci (eQTL) mapping is the search 
for genomic loci associated with gene expression. Traditional 
methods perform eQTL mapping on bulk tissue gene expression 
without accounting for the underlying sources of gene expression 
per cell type. Attempts to perform ``cell type-aware'' eQTL 
mapping commonly use linear models that transform the gene 
expression while adjusting for the genotype, a proposed measure 
for cell type prevalence (e.g. a marker gene's expression), and 
their interaction and proceed to perform hypothesis testing.

Our R package performs cell type-specific eQTL mapping on bulk 
RNA-sequencing samples by jointly modeling total read count (TReC) 
and allele-specific read count (ASReC). Our method, **CSeQTL**, 
introduces novel features and fully extends **TReCASE** 
[[HTML](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218220/), 
[PDF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218220/pdf/nihms-307768.pdf), 
[Software](https://github.com/Sun-lab/asSeq)] and **pTReCASE** 
[[HTML](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7410098/), 
[PDF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7410098/pdf/nihms-1028292.pdf), 
[Software](https://github.com/Sun-lab/pTReCASE)] methodologies.

<p align="center">
<a href="https://raw.githubusercontent.com/pllittle/CSeQTL/master/images/ex_CSeQTL_reads.png"><img src="images/ex_CSeQTL_reads.png" width="50%" /></a>
<p align="center"><em>Example of total read counts and allele-specific or haplotype read counts.</em></p>
</p>

## Required Input Data

1. TReC data
2. ASReC data (haplotype 1 and 2 read counts)
3. Phased genotypes
4. Sample covariates: 
	 * observed confounders, 
	 * genotype principal components (PCs),
	 * latent batch effects (derived from residual TReC PCs)
5. Cell type proportions

## Installation

<details>
<summary>Click to expand!</summary>

```R
req_packs = c("devtools","Rcpp","RcppArmadillo",
	"smarter","CSeQTL")
all_packs = as.character(installed.packages()[,1])
rerun = 0
build_vign = ifelse(Sys.getenv("RSTUDIO_PANDOC") == "",FALSE,TRUE)

for(pack in req_packs){
	if( pack %in% all_packs ){
		library(package = pack,character.only = TRUE)
		next
	}
	
	bb = NULL
	if( pack %in% c("smartr","CSeQTL") ){
		repo = sprintf("pllittle/%s",pack)
		bb = tryCatch(devtools::install_github(repo = repo,
			build_vignettes = build_vign,
			dependencies = TRUE),
			error = function(ee){"error"})
	} else {
		bb = tryCatch(install.packages(pkgs = pack,
			dependencies = TRUE),
			error = function(ee){"error"})
	}
	
	if( !is.null(bb) && bb == "error" )
		stop(sprintf("Error for package = %s",pack))
	rerun = 1
}

if( rerun == 1 ) stop("Re-run above code")
```

</details>

## Vignette

```R
# An Introduction
vignette(topic = "intro",package = "CSeQTL")
```

## Citation
Little, P., [Zhabotynsky, V.](https://github.com/yaceya), 
[Li, Y.](https://github.com/yunliUNC), 
[Lin, D.Y.](https://sph.unc.edu/adv_profile/danyu-lin-phd/), 
[Sun, W.](https://github.com/sunway1999) (2022). Cell 
type-specific Expression Quantitative Trait Loci. *bioRxiv*. 
[[HTML](https://www.biorxiv.org/content/10.1101/2022.03.31.486605v1), 
[PDF](https://www.biorxiv.org/content/10.1101/2022.03.31.486605v1.full.pdf)]



