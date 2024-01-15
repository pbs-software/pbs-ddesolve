\name{dot-onClosePBSddeExamples}
\alias{.onClosePBSddeExamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  On Close Set Old WD
}
\description{
  A trivial function that sets the user's working directory to an old
  (previous) location before opening the Windows GUI that is now 
  being closed.
}
\usage{
.onClosePBSddeExamples()
}
%- maybe also 'usage' for other objects documented here.
\value{
  \code{setwd(.PBSddeEnv$oldwd)}
}
\author{
  \href{mailto:alex@mofo.ca}{Alex Couture-Beil} -- Software Engineer, Earthly Technologies, Victoria BC\cr

  Maintainer: \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Offsite, Vancouver BC\cr
  Last modified \code{Rd: 2024-01-04}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{iplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{GUI}
% \concept{ ~cpt2 }
% Use only one concept per line.
