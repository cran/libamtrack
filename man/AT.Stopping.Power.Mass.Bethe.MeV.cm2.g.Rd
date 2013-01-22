% TODO File path/AT.Stopping.Power.Mass.Bethe.MeV.cm2.g.Rd
\name{AT.Stopping.Power.Mass.Bethe.MeV.cm2.g}
\alias{AT.Stopping.Power.Mass.Bethe.MeV.cm2.g}
\title{AT.Stopping.Power.Mass.Bethe.MeV.cm2.g}
\description{Computes the mass stopping power using the Bethe formula
 for many particles according to ICRU49, p.6, Eq. 2.1
 BUT WITHOUT shell or density, Bloch or Barkas correction!}
\usage{AT.Stopping.Power.Mass.Bethe.MeV.cm2.g(E.MeV.u, particle.no, material.no, E.restricted.keV)
}
\arguments{
  \item{E.MeV.u}{ energies of particle per nucleon (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle indices (array of size n) (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (single value) (see also \code{\link{material.no}}).}
  \item{E.restricted.keV}{ if positive and smaller than maximally transferable energy, the restricted stopping power will be computed (single value).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{Mass.Stopping.Power.MeV.cm2.g}{ (array of size n)}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataStoppingPower.c#L454}
}


