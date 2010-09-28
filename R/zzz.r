.First.lib <- function(lib,pkg)
{
	pkg_info <- utils::sessionInfo( package=pkg )$otherPkgs[[ pkg ]]
	pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]

	userguide_path <- system.file( "doc/ddesolveIntro.pdf", package = pkg )
	
	cat("
ddesolve", pkg_info$Version, "-- Copyright (C) 2007-2010 Fisheries & Oceans Canada
(based  on solv95 by Simon Wood)

--------------------------------------------------
*****  DEPRECATED - use package PBSddesolve  *****
--------------------------------------------------

Type 'openVig()' on the command line to see the vignette:
", userguide_path, "

Packaged on", pkg_date, "
Pacific Biological Station, Nanaimo


")
}

