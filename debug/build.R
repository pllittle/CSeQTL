# Create package
rm(list = ls())
git_dir = file.path(getwd(),"../../")
pack = "CSeQTL"
pack_dir = file.path(git_dir,pack)
if( !dir.exists(pack_dir) ) q("no")
chk_pack = tryCatch(find.package(pack),
	warning = function(ww){NULL},
	error = function(ee){NULL})
if( !is.null(chk_pack) ){
	remove.packages(pack)
	q("no")
}

req_packs = c("Rcpp","devtools","rmarkdown")
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

compileAttributes(pkgdir = pack_dir)
document(pkg = pack_dir)
use_gpl3_license()

Sys.getenv("_R_CHECK_SYSTEM_CLOCK_")
Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

make_vign = TRUE
status = find_pandoc()
check_pandoc = pandoc_available(version = status$version)
chk_vign = length(list.files(file.path(pack_dir,"vignettes"),pattern = ".Rmd$")) > 0
make_vign = check_pandoc && chk_vign && make_vign; make_vign

chk = tryCatch(check(pkg = pack_dir,
	manual = TRUE,cran = TRUE,
	vignettes = make_vign,
	error_on = "note"),
	warning = function(ww){NULL},
	error = function(ee){NULL})

# Install locally
if( !is.null(chk) ){
	install(pkg = pack_dir,upgrade = !TRUE,
		build_vignettes = make_vign)
}


##

