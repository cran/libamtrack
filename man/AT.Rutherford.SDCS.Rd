% TODO File path/AT.Rutherford.SDCS.Rd
\name{AT.Rutherford.SDCS}
\alias{AT.Rutherford.SDCS}
\title{AT.Rutherford.SDCS}
\description{Computes the Rutherford single differential cross section
 for the energy spectrum of secondary electrons produced by
 an HCP}
\usage{AT.Rutherford.SDCS(E.MeV.u, particle.no, material.no, T.MeV)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per nucleon (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{T.MeV}{ electron energies (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{dsdT.m2.MeV}{ Rutherford SDCS for given electron energies (array of
 size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Phy
sicsRoutines.c#L805}
}
