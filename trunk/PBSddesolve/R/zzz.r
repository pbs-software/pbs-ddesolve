
.First.lib <- function(lib, pkg)
{
	library.dynam("PBSddesolve", pkg, lib)
	cat("
PBSddesolve 1.05 -- based on solv95 by Simon Wood

A complete user guide 'PBSddesolve-UG.pdf' appears 
in the '.../library/PBSddesolve/doc' folder.

Built on Sep 19, 2008
Pacific Biological Station, Nanaimo\n\n")
}
