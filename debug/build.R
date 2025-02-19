# Steps to create/check/install package from directory

rm(list = ls())
git_dir = sprintf("%s/../..",getwd()); git_dir
pack = "CSeQTL"
pack_dir = file.path(git_dir,pack)
if( !file.exists(pack_dir) ) q("no")

chk_pack = tryCatch(find.package(pack),
	error = function(ee){NULL}); chk_pack
if( !is.null(chk_pack) ){
	remove.packages(pack)
	q("no")
}

req_packs = c("usethis","Rcpp","devtools","rmarkdown")
for(pp in req_packs){
	
	chk_pack = tryCatch(find.package(pp),
		warning = function(ww){NULL},
		error = function(ee){NULL})
	
	if( !is.null(chk_pack) ){
		library(pp,character.only = TRUE)
		next
	}
	
	stop(sprintf("Install %s",pp))
}

Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

compileAttributes(pkgdir = pack_dir)
document(pkg = pack_dir)
use_gpl3_license()
check_pandoc = pandoc_available(); check_pandoc
vign_dir = file.path(pack_dir,"vignettes")
chk_vign = length(list.files(vign_dir,pattern = "Rmd")) > 0
override_vign = !FALSE
make_vign = check_pandoc && chk_vign && override_vign; make_vign

# Check: takes some time
perf_chk = TRUE
if( perf_chk ){
	chk = tryCatch(check(pkg = pack_dir,
		manual = TRUE,cran = TRUE,
		error_on = "note",
		vignettes = make_vign),
		error = function(ee){NULL},
		warning = function(ww){NULL})
} else {
	chk = TRUE
}
chk

# Install locally
if( !is.null(chk) ){
	chk = tryCatch(install(pack_dir,
		build_vignettes = make_vign,
		upgrade = FALSE),
		error = function(ee){NULL},
		warning = function(ww){NULL})
}

# Build
if( !is.null(chk) ){
	desk_dir = sprintf("%s/../",git_dir)
	devtools::build(pkg = pack_dir,path = desk_dir)
}

##

