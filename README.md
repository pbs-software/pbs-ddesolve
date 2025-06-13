## PBSddesolve: Solver for delay differential equations ##
&copy; Fisheries and Oceans Canada (2007-2025)

The R package **PBSddesolve** generates numerical solutions for systems of delay differential equations (DDEs) and ordinary differential equations (ODEs). The numerical routines come from Simon Wood’s program `solv95`, originally written in C for the Microsoft Windows operating systems. With **PBSddesolve**, a user can write the gradient code for a system of DDEs or ODEs in the R language, rather than C. The code will then run on all platforms supported by R, and the results can be inspected using R’s extensive graphics capabilities. Simon has very generously given us permission to publish **PBSddesolve** (including his embedded routines) under the GNU GENERAL PUBLIC LICENSE Version 2. 

For more information, see the User's Guide (file `PBSddesolve-UG.pdf`) featured in this repository. Obtain the package **PBSddesolve** from <a href="https://CRAN.R-project.org/package=PBSddesolve">CRAN</a>, the standard repository for R contributed packages. 

This package replaces an earlier R package **ddesolve**, which is no longer supported.

**PBSddesolve** represents just one of a series of R packages developed at the Pacific Biological Station (<a href="http://www.pac.dfo-mpo.gc.ca/science/facilities-installations/index-eng.html#pbs">PBS</a>) in Nanaimo, British Columbia. A more advanced version of **PBSddesolve** might be available at <a href="https://github.com/pbs-software">pbs-software on GitHub</a>. Any evolving package (Windows binary and source tarball) is built after using CRAN's rigorous `R CMD check --as-cran` routine (on a Windows system). Most of the time, the revision on <a href="https://github.com/pbs-software/pbs-ddesolve">GitHub</a> can be built in R using `devtools::install_github("pbs-software/pbs-ddesolve/PBSddesolve")`; however, not every revision has been checked for CRAN worthiness.

<font color="red"><h3>Disclaimer</h3></font>

"Fisheries and Oceans Canada (DFO) GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. DFO relinquishes control of the information and assumes no responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against DFO stemming from the use of its GitHub project will be governed by all applicable Canadian Federal laws. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favouring by DFO. The Fisheries and Oceans Canada seal and logo, or the seal and logo of a DFO bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DFO or the Canadian Government.”

As with any freely available product, there is no warranty or promise that **PBSddesolve** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 
