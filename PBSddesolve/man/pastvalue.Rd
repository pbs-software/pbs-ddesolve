\name{pastvalue}
\alias{pastvalue}
\alias{pastgradient}
\title{
  Retrieve Past Values (lags) During Gradient Calculation
}
\description{
   These routines provides access to variable history at lagged times.
   The lagged time \eqn{t} must not be less than \eqn{t_0}{t0}, nor
   should it be greater than the current time of gradient calculation.
   The routine cannot be directly called by a user, and will only work
   during the integration process as triggered by the \code{dde} routine.
}
\usage{
pastvalue(t)
pastgradient(t)
}
\arguments{
  \item{t}{Access history at time \code{t}.}
}
\value{
  Vector of variable history at time \code{t}.
}
\author{
  Alex Couture-Beil <alex@mofo.ca>

  Maintainer: \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-03-07}
}
\seealso{
  \code{\link[PBSddesolve]{dde}}
}
\keyword{math}

