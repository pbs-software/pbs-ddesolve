# Taking cue from Roger Bivand's maptools:
.PBSddeEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onLoad <- function(lib,pkg)
{
	library.dynam("PBSddesolve", pkg, lib)
}

.onAttach <- function(lib, pkg)
{
        # obtain values necessary for the welcome message
	pkg_info <- utils::sessionInfo( package="PBSddesolve" )$otherPkgs$PBSddesolve
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	userguide_path <- system.file( "doc/PBSddesolve-UG.pdf", package = "PBSddesolve" )
	year <- substring(date(),nchar(date())-3,nchar(date()))
	
	packageStartupMessage("
-----------------------------------------------------------
PBSddesolve ", pkg_info$Version, " -- Copyright (C) 2007-",year," Fisheries and Oceans Canada
(based  on solv95 by Simon Wood)

A complete user guide 'PBSddesolve-UG.pdf' is located at 
", userguide_path, "

Demos include 'blowflies', 'cooling', 'icecream', and 'lorenz'
They can be run two ways:

1. Using 'utils' package 'demo' function, run command
> demo(icecream, package='PBSddesolve', ask=FALSE)

2. Using package 'PBSmodelling', run commands
> require(PBSmodelling); runDemos()
and choose PBSddesolve and then one of the four demos.

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSddeEnv)
}

