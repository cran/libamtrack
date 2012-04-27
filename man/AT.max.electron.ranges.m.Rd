% TODO File path/AT.max.electron.ranges.m.Rd
\name{AT.max.electron.ranges.m}
\alias{AT.max.electron.ranges.m}
\title{AT.max.electron.ranges.m}
\description{Returns the maximum electron range (track radius) in m
 for a given parametrization
}
\usage{AT.max.electron.ranges.m(E.MeV.u, material.no, er.model)
}
\arguments{
  \item{E.MeV.u}{ kinetic energy for particles in the given field (array of size number.of.particles) (see also \code{\link{E.MeV.u}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{er.model}{ electron-range model index (see also \code{\link{er.model}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{max.electron.range.m}{ electron range (track radius) in m (array of size number.of.particles)}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.3/src/AT_ElectronRange.c#L126}
}
\examples{
# Compute the electron range in PMMA for the Tabata parametrization between 0.3 keV and 30 MeV
AT.max.electron.ranges.m( E.MeV.u     = 0.3 * 10^(-4:2),
                          material.no = 4,
                          er.model    = 7)
}


