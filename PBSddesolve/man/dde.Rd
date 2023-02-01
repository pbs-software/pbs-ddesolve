\name{dde}
\alias{dde}
\title{Solve Delay Differential Equations}

\description{
  A solver for systems of delay differential equations based on numerical routines from
  C source code \code{solv95} by 
  \href{https://www.maths.ed.ac.uk/~swood34/simon/dde.html}{Simon Wood}.
  This solver is also capable of solving systems of ordinary differential equations. 
}
\usage{
dde(y, times, func, parms=NULL, switchfunc=NULL, mapfunc=NULL, 
   tol=1e-08, dt=0.1, hbsize=10000) 
}
\arguments{
  \item{y}{\code{numeric} -- vector of initial values of the DDE system.
    The size of the supplied vector determines the number of variables in the system.}
  \item{times}{\code{numeric} -- vector of specific times to solve. }
  \item{func}{\code{function} -- a user-supplied function that computes the gradients in 
    the DDE system at time \code{t}.
    The function must be defined using the arguments: \code{(t,y)} or \code{(t,y,parms)},
    where \code{t} is the current time in the integration, \code{y} is a vector of the 
    current estimated variables of the DDE system, and \code{parms} is any R object 
    representing additional parameters (optional).\cr
    The argument \code{func} must return one of the two following return types:\cr
    1) a vector containing the calculated gradients for each variable; or \cr
    2) a list with two elements - the first a vector of calculated gradients, the second 
    a vector (possibly named) of values for a variable specified by the user at each 
    point in the integration.}
  \item{parms}{\code{list} -- any constant parameters to pass to \code{func},
    \code{switchfunc}, and \code{mapfunc}.}
  \item{switchfunc}{\code{function} -- an optional function that is used to manipulate
    state values at given times. 
    The switch function takes the arguments \code{(t,y)} or \code{(t,y,parms)} and must 
    return a numeric vector. 
    The size of the vector determines the number of switches used by the model. 
    As values of \code{switchfunc} pass through zero (from positive to negative), a 
    corresponding call to \code{mapfunc} is made, which can then modify any state value.}
  \item{mapfunc}{\code{function} -- if \code{switchfunc} is defined, then a map function
    must also be supplied with arguments \code{(t, y, switch_id)} or 
    \code{t, y, switch_id, parms)}, where \code{t} is the time, \code{y} are the current
    state values, \code{switch_id} is the index of the triggered switch, and \code{parms}
    are additional constant parameters.}
  \item{tol}{\code{numeric} -- maximum error tolerated at each time step (as a proportion
    of the state variable concerned).}
  \item{dt}{\code{numeric} -- maximum initial time step.}
  \item{hbsize}{\code{numeric} -- history buffer size required for solving DDEs.}
}
\details{
  Please see the included demos (\code{'blowflies'}, \code{'cooling'}, \code{'icecream'},
  \code{'lorenz'}) for examples of how to use \code{dde}.\cr\cr
  The demos can be run two ways:
  \enumerate{
    \item Using the package \code{utils}, run the command:\cr
    \code{demo(icecream, package="PBSddesolve", ask=FALSE)}
    \item Using the package \code{PBSmodelling}, run the commands:\cr
    \code{require(PBSmodelling); runDemos()}
  }
  The latter produces a GUI that shows all demos available from locally installed packages.
  Choose \code{PBSddesolve}. 
  Note that the examples are run in the temporary working environment \code{.PBSddeEnv}.

  The user supplied function \code{func} can access past values (lags) of \code{y} by
  calling the \code{\link[PBSddesolve]{pastvalue}} function. 
  Past gradients are accessible by the \code{\link[PBSddesolve]{pastgradient}} function.
  These functions can only be called from \code{func} and can only be passed values of
  \code{t} greater or equal to the start time, but less than the current time of the 
  integration point.
  For example, calling \code{pastvalue(t)} is not allowed, since these values are the
  current values which are passed in as \code{y}.
}
\value{
  A data frame with one column for \code{t}, a column for every variable in the system,
  and a column for every additional value that may (or may not) have been returned by
  \code{func} in the second element of the list.

  If the initial \code{y} values parameter was named, then the solved values column will
  use the same names. Otherwise \code{y1}, \code{y2}, \dots{} will be used.

  If \code{func} returned a list, with a named vector as the second element, then those
  names will be used as the column names. 
  If the vector was not named, then \code{extra1}, \code{extra2}, \dots{} will be used.
}
\author{
  \href{mailto:alex@mofo.ca}{Alex Couture-Beil} -- Software Engineer, Earthly Technologies, Victoria BC\cr

  Maintainer: \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-01-30}
}
\seealso{
  \code{\link[PBSddesolve]{pastvalue}}
}
\examples{
##################################################
## This is just a single example of using dde.
## For more examples see demo(package="PBSddesolve")
## the demos require the package PBSmodelling
##################################################

require(PBSddesolve)
local(env=.PBSddeEnv, expr={
  #create a func to return dde gradient
  yprime <- function(t,y,parms) {
    if (t < parms$tau)
      lag <- parms$initial
    else
      lag <- pastvalue(t - parms$tau)
    y1 <- parms$a * y[1] - (y[1]^3/3) + parms$m * (lag[1] - y[1])
    y2 <- y[1] - y[2]
    return(c(y1,y2))
  }

  #define initial values and parameters
  yinit <- c(1,1)
  parms <- list(tau=3, a=2, m=-10, initial=yinit)

  # solve the dde system
  yout <- dde(y=yinit,times=seq(0,30,0.1),func=yprime,parms=parms)

  # and display the results
  plot(yout$time, yout$y1, type="l", col="red", xlab="t", ylab="y", 
    ylim=c(min(yout$y1, yout$y2), max(yout$y1, yout$y2)))
  lines(yout$time, yout$y2, col="blue")
  legend("topleft", legend = c("y1", "y2"),lwd=2, lty = 1, 
    xjust = 1, yjust = 1, col = c("red","blue"))
})
}
\keyword{math}
\concept{delay differential}

