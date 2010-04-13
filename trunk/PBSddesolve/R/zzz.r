
.First.lib <- function(lib, pkg)
{
	library.dynam("PBSddesolve", pkg, lib)
	cat("
PBSddesolve 1.06 -- based on solv95 by Simon Wood

A complete user guide 'PBSddesolve-UG.pdf' appears 
in the '.../library/PBSddesolve/doc' folder.

Built on Apr 13, 2010
Pacific Biological Station, Nanaimo\n\n")
}
