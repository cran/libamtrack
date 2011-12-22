% TODO File path/AT.Stopping.Power.keV.um.Rd
\name{AT.Stopping.Power.keV.um}
\alias{AT.Stopping.Power.keV.um}
\title{AT.Stopping.Power.keV.um}
\description{Main access method to stopping power data - multiple field}
\usage{AT.Stopping.Power.keV.um(stopping.power.source.no, E.MeV.u, particle.no, material.no)
}
\arguments{
  \item{stopping.power.source.no}{ (see also \code{\link{stopping.power.source.no}}).}
  \item{E.MeV.u}{ (array of size number.of.particles) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ (array of size number.of.particles) (see also \code{\link{particle.no}}).}
  \item{material.no}{ (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{Stopping.Power.keV.um}{ (array of size number.of.particles)}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.2/src/AT_DataStoppingPower.c#L211}
}


