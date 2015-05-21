% TODO File path/AT.gamma.from.E.Rd
\name{AT.gamma.from.E}
\alias{AT.gamma.from.E}
\title{AT.gamma.from.E}
\description{Returns relativistic gamma
}
\usage{AT.gamma.from.E(E.MeV.u)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n) (see also \code{\link{E.MeV.u}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{gamma}{ vector of results (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Phy
sicsRoutines.c#L52}
}
\examples{
# Compute mass in MeV/c2 for a proton with 1-1000 MeV kinetic energy:
938.3 * AT.gamma.from.E( E.MeV.u = 10^(0:3))$gamma
}
