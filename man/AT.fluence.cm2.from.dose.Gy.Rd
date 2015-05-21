% TODO File path/AT.fluence.cm2.from.dose.Gy.Rd
\name{AT.fluence.cm2.from.dose.Gy}
\alias{AT.fluence.cm2.from.dose.Gy}
\title{AT.fluence.cm2.from.dose.Gy}
\description{Returns fluence in 1/cm2 for each given particle}
\usage{AT.fluence.cm2.from.dose.Gy(E.MeV.u, particle.no, D.Gy, material.no,
 stopping.power.source.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of
 size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array
 of size n) (see also \code{\link{particle.no}}).}
  \item{D.Gy}{ dose / Gy for each particle type (array of size n).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{fluence.cm2}{ to be allocated by the user which will be used to return
 the results (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Phy
sicsRoutines.c#L333}
}
