#' @importFrom graphics boxplot text segments plot par 
#'	polygon points abline layout legend mtext hist
#' @importFrom stats rbinom rgamma runif pchisq rnbinom
#'	rbeta rpois rnorm qnorm lm residuals predict 
#'	cooks.distance median mad na.omit quantile
#' @importFrom utils head
#' @importFrom grDevices rgb
#' @importFrom MatrixEQTL modelLINEAR Matrix_eQTL_main SlicedData
#' @importFrom emdbook rbetabinom
#' @importFrom data.table fread fwrite
#' @importFrom multcomp glht adjusted
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab labs theme 
#'	element_text ggsave position_dodge element_blank unit
#'	element_rect geom_errorbar geom_point facet_grid 
#'	scale_shape_manual geom_hline coord_flip guides guide_legend
#' @importFrom HelpersMG wget
#' @importFrom smarter smart_names smart_df smart_table
#'	smart_progress name_change smart_rmcols calc_JK
#' @importFrom Rcpp sourceCpp
#' @useDynLib CSeQTL
NULL

# Create package
# rm(list=ls()); library(smarter)
# smart_prepPack(pack_dir = "C:/Users/Admin/Desktop/github/CSeQTL",
#		pandoc = "C:/Program Files/RStudio/bin/pandoc",
#		make_vign = TRUE,cran = FALSE,
#		build_dir = NULL)
#		build_dir = "C:/Users/Admin/Desktop")

###
