\name{pastvalue}
\alias{pastvalue}
\alias{pastgradient}
\title{
  Retrieve Past Values (lags) During Gradient Calculation
}
\description{
  These routines provides access to variable history at lagged times.
  The lagged time \eqn{t} must not be less than \eqn{t_0}{t0}, nor should it be greater
  than the current time of gradient calculation.
  The routine cannot be directly called by a user, and will only work during the 
  integration process as triggered by the \code{dde} routine.
}
\usage{
pastvalue(t)
pastgradient(t)
}
\arguments{
  \item{t}{\code{numeric} -- time \emph{t} at which history is accessed}
}
\value{
  Vector of variable history at time \emph{t}.
}
\author{
  \href{mailto:alex@mofo.ca}{Alex Couture-Beil} -- Software Engineer, Earthly Technologies, Victoria BC

  Maintainer: \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-13}
}
\seealso{
In package \pkg{PBSddesolve}:\cr
  \code{\link[PBSddesolve]{dde}}
}
\keyword{math}


