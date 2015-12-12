% TODO File path/AT.Stopping.Power.with.no.Rd
\name{AT.Stopping.Power.with.no}
\alias{AT.Stopping.Power.with.no}
\title{AT.Stopping.Power.with.no}
\description{Retrieves the electronic stopping power in keV/um for
 the requested energies and particles for a specified
 material and data source. The data source is thereby
 given via its integer id (s. AT_StoppingPowerData.h for
 details). Data that should be read from a file
 cannot be used with this method.
}
\usage{AT.Stopping.Power.with.no(stopping.power.source.no, E.MeV.u,
 particle.no, material.no)
}
\arguments{
  \item{stopping.power.source.no}{ id of the data source (see also
 \code{\link{stopping.power.source.no}}).}
  \item{E.MeV.u}{ kinetic energies in MeV per amu (array of size n) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle numbers (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material number (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{stopping.power.keV.um}{ array to return stopping powers (array of size
 n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Sto
ppingPower.c#L128}
}
