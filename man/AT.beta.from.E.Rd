% TODO File path/AT.beta.from.E.Rd
\name{AT.beta.from.E}
\alias{AT.beta.from.E}
\title{AT.beta.from.E}
\description{Returns relativistic speed for many particles
}
\usage{AT.beta.from.E(E.MeV.u)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n) (see also \code{\link{E.MeV.u}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{beta}{ vector of relative particle speed beta = v/c (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Phy
sicsRoutines.c#L34}
}
\examples{
# Compute beta between 1 and 1000 MeV/u:
AT.beta.from.E( E.MeV.u     = c(1,10,100,1000))
}
