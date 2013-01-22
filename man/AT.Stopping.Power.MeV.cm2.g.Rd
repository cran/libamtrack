% TODO File path/AT.Stopping.Power.MeV.cm2.g.Rd
\name{AT.Stopping.Power.MeV.cm2.g}
\alias{AT.Stopping.Power.MeV.cm2.g}
\title{AT.Stopping.Power.MeV.cm2.g}
\description{Main access method to stopping power data - multiple field}
\usage{AT.Stopping.Power.MeV.cm2.g(stopping.power.source.no, E.MeV.u, particle.no, material.no)
}
\arguments{
  \item{stopping.power.source.no}{ (see also \code{\link{stopping.power.source.no}}).}
  \item{E.MeV.u}{ (array of size number.of.particles) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ (array of size number.of.particles) (see also \code{\link{particle.no}}).}
  \item{material.no}{ (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{Stopping.Power.MeV.cm2.g}{ (array of size number.of.particles)}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataStoppingPower.c#L193}
}


