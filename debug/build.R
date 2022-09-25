# Create package
rm(list=ls())
library(smarter)
user_dir = gsub("\\\\","/",Sys.getenv("USERPROFILE"))
if( user_dir == "" ){
	git_dir = "~/github"
} else {
	git_dir = file.path(user_dir,"Desktop/github")
}

pack			= "CSeQTL"
pack_dir 	= file.path(git_dir,pack)
all_packs	= as.character(installed.packages()[,1])
if( pack %in% all_packs ){
	remove.packages(pack)
	q("no")
}

if( user_dir == "" ){
	devtools::install_github("pllittle/CSeQTL")
} else {
	smart_prepPack(pack_dir = pack_dir,
		pandoc = "C:/Program Files/RStudio/bin/pandoc",
		make_vign = FALSE,
		cran = FALSE,
		build_dir = NULL)
}

##

