\name{PBSddesolve}
\alias{PBSddesolve-package}
\alias{PBSddesolve}
\docType{package}
\title{
  Package: Solver for Delay Differential Equations
}
\description{
  A solver for systems of delay differential equations based on 
  numerical routines from Simon Wood's
  \href{https://webhomes.maths.ed.ac.uk/~swood34/simon/dde.html}{\code{solv95}} programme. 
  This solver is also capable of solving systems of ordinary
  differential equations.
}
\details{
  Please see the user guide \code{PBSddesolve-UG.pdf}, located in R's
  library directory \code{./library/PBSddesolve/doc}, for a 
  comprehensive overview.
}
\author{
  \href{mailto:simon.wood@ed.ac.uk}{Simon Wood} -- Chair of Computational Statistics, School of Mathematics, University of Edinburgh\cr
  \href{mailto:alex@mofo.ca}{Alex Couture-Beil} -- Software Engineer, Earthly Technologies, Victoria BC\cr
  \href{mailto:schnutej@shaw.ca}{Jon T. Schnute} -- Scientist Emeritus, DFO, Nanaimo BC\cr
  \href{mailto:nboers@gmail.com}{Nicholas M. Boers} -- Senior Software Engineer, Jobber, Edmonton AB

  Maintainer: \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-13}
}
\references{
  Wood, S.N. (1999) 
  Solv95: a numerical solver for systems of delay differential equations with switches.
  Mathematical Institute, North Haugh, St. Andrews, Fife KY16 9SS, U.K., 10 p.
}

\seealso{
In package \pkg{PBSddesolve}:\cr
  \code{\link[PBSddesolve]{dde}}\cr
In package \pkg{deSolve}:\cr
  \code{\link[deSolve]{lsoda}}
}
\keyword{package}
\concept{DDE}

