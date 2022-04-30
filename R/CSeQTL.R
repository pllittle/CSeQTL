#' @importFrom graphics boxplot text segments plot par 
#'	polygon points abline layout legend mtext hist
#' @importFrom stats rbinom rgamma runif pchisq rnbinom
#'	rbeta rpois rnorm qnorm lm residuals predict 
#'	cooks.distance median mad na.omit quantile
#' @importFrom utils head
#' @importFrom grDevices rgb
#' @importFrom MatrixEQTL modelLINEAR Matrix_eQTL_main SlicedData
#' @importFrom emdbook rbetabinom
#' @importFrom data.table fread
#' @importFrom multcomp glht adjusted
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab labs theme 
#'	element_text ggsave position_dodge
#' @importFrom smarter smart_names smart_df smart_table
#'	smart_progress
#' @importFrom Rcpp sourceCpp
#' @useDynLib CSeQTL
NULL


# Steps to create/check/install package from directory
# bb = strsplit(getwd(),"/")[[1]]; pack_dir = paste(bb[-length(bb)],collapse = "/")
# pack = strsplit(pack_dir,"/")[[1]]; pack = pack[length(pack)]; pack
# if( pack %in% installed.packages()[,1] ){ remove.packages(pack); q("no")}
# Rcpp::compileAttributes(pkgdir = pack_dir)
# devtools::document(pkg = pack_dir); usethis::use_gpl3_license()
# Sys.setenv("RSTUDIO_PANDOC" = "C:/Program Files/RStudio/bin/pandoc")
# check_pandoc = rmarkdown::pandoc_available(); check_pandoc
#### usethis::use_vignette(name = "test",title = "Testing")
# make_vign = check_pandoc && TRUE
# devtools::check(pkg = pack_dir,manual = TRUE,cran = !TRUE,error_on = c("warning","note")[1],vignettes = make_vign)
# devtools::install(pack_dir,build_vignettes = make_vign)

