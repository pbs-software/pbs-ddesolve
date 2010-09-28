\name{openVig}
\alias{openVig}
\title{
  Open a Vignette pkgIntro.pdf
}
\description{
  Open a package vignette called \emph{pkg}\code{Intro.pdf} if it 
  exists in the \code{doc} subdirectory of the installed package.
}
\usage{
openVig(pkg = "ddesolve")
}
\arguments{
  \item{pkg}{character string that denotes an R package.}
}
\value{
  If vignette \code{pdf} is found, the system's shell will attempt to open it.
}
\author{
  Rowan Haigh, Pacific Biological Station, Nanaimo BC
}
\keyword{package}

