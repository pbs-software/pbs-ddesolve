\name{pastvalue}
\alias{pastvalue}
\alias{pastgradient}
\title{Retrieve Past Values (lags) During Gradient Calculation}

\description{
   These routines provides access to variable history at lagged times.
   The lagged time \eqn{t} must not be less than \eqn{t_0}{t0}, nor should 
   it be greater than the current time of gradient calculation. The routine 
   cannot be directly called by a user, and will only work during the
   integration process as triggered by the \code{dde} routine.
}

\usage{ 
pastvalue(t)
pastgradient(t)
}

\arguments{
\item{t}{access history at time \code{t}.}
}

\value{ vector of variable history at time \code{t}. }

\seealso{ \code{\link{dde}} }

\keyword{math}
