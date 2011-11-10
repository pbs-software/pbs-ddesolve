.onLoad <- function(lib,pkg)
{
	library.dynam("PBSddesolve", pkg, lib);
	pkg_info <- utils::sessionInfo( package="PBSddesolve" )$otherPkgs$PBSddesolve
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- "unkown"
	
	userguide_path <- system.file( "doc/PBSddesolve-UG.pdf", package = "PBSddesolve" )
	
	cat("
-----------------------------------------------------------
PBSddesolve", pkg_info$Version, "-- Copyright (C) 2007-2011 Fisheries and Oceans Canada
(based  on solv95 by Simon Wood)

A complete user guide 'PBSddesolve-UG.pdf' is located at 
", userguide_path, "

Packaged on", pkg_date, "
Pacific Biological Station, Nanaimo
-----------------------------------------------------------


")
}

