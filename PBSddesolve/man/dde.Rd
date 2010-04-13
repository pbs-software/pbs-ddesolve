\name{dde}
\alias{dde}
\title{Solve Delay Differential Equations}

\description{
  A solver for systems of delay differential equations based off numerical 
  routines from Simon Wood's \emph{solv95} program. This solver is also 
  capable of solving systems of ordinary differential equations. 
  
  Please see the included demos for examples of how to use \code{dde}.\cr 
  To view available demos run \code{demo(package="PBSddesolve")}.\cr
  The supplied demos require that the package \pkg{PBSmodelling} be installed.
}

\usage{
dde(y, times, func, parms=NULL, switchfunc=NULL, mapfunc=NULL, 
    tol=1e-08, dt=0.1, hbsize=10000) 
}

\arguments{
\item{y}{vector of initial values of the DDE system. The size of the supplied vector
    determines the number of variables in the system.}
\item{times}{numeric vector of specific times to solve. }
\item{func}{a user supplied function that computes the gradients in the DDE 
   system at time \code{t}. The function must be defined using the arguments: 
   \code{(t,y)} or \code{(t,y,parms)}, where \code{t} is the current time
   in the integration, \code{y} is a vector of the current estimated variables 
   of the DDE system, and \code{parms} is any R object representing 
   additional parameters (optional). 

   The argument \code{func} must return one of the two following return types:
   1) a vector containing the calculated gradients for each variable; or 
   2) a list with two elements - the first a vector of calculated gradients, 
   the second a vector (possibly named) of values for a variable specified by 
   the user at each point in the integration.
}
\item{parms}{any constant parameters to pass to \code{func}, \code{switchfunc}, 
   and \code{mapfunc}.
}
\item{switchfunc}{an optional function that is used to manipulate state 
   values at given times. The switch function takes the arguments \code{(t,y)} or 
   \code{(t,y,parms)} and must return a numeric vector. The size of the vector 
   determines the number of switches used by the model. As values of 
   \code{switchfunc} pass through zero (from positive to negative), a corresponding 
   call to \code{mapfunc} is made, which can then modify any state value.
}
\item{mapfunc}{if \code{switchfunc} is defined, then a map function must also 
   be supplied with arguments \cr \code{(t,y,switch_id)} or \code{t,y,switch_id,parms)},
   where \code{t} is the time, \code{y} are the current state values, \code{switch_id} 
   is the index of the triggered switch, and \code{parms} are additional constant parameters.
}
\item{tol}{maximum error tolerated at each time step (as a proportion of 
   the state variable concerned)
}
\item{dt}{maximum initial time step}
\item{hbsize}{history buffer size required for solving DDEs)}
}

\details{
   The user supplied function \code{func} can access past values (lags) of 
   \code{y} by calling the \code{\link{pastvalue}} function. Past gradients 
   are accessible by the \code{\link{pastgradient}} function. These functions 
   can only be called from \code{func} and can only be passed values of 
   \code{t} greater or equal to the start time, but less than the current time 
   of the integration point. For example, calling \code{pastvalue(t)} is not 
   allowed, since these values are the current values which are passed in as 
   \code{y}.
}

\value{
  A data frame with one column for \code{t}, a column for every variable in the system,
  and a column for every additional value that may (or may not) have been returned
  by \code{func} in the second element of the list.
  
  If the initial \code{y} values parameter was named, then the solved values column 
  will use the same names. Otherwise \code{y1}, \code{y2}, ... will be used.
  
  If \code{func} returned a list, with a named vector as the second element, then
  those names will be used as the column names. If the vector was not named, then
  \code{extra1}, \code{extra2}, ... will be used.
}

\seealso{ \code{\link{pastvalue}} }

\examples{
##################################################
# This is just a single example of using dde.
# For more examples see demo(package="PBSddesolve")
# the demos require the package PBSmodelling
##################################################

#create a func to return dde gradient
require(PBSddesolve)
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
plot(yout$t, yout$y1, type="l", col="red", xlab="t", ylab="y", 
     ylim=c(min(yout$y1, yout$y2), max(yout$y1, yout$y2)))
lines(yout$t, yout$y2, col="blue")
legend("topleft", legend = c("y1", "y2"),lwd=2, lty = 1, 
       xjust = 1, yjust = 1, col = c("red","blue"))
}

\keyword{math}
