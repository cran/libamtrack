% TODO File path/AT.dose.Gy.from.fluence.cm2.Rd
\name{AT.dose.Gy.from.fluence.cm2}
\alias{AT.dose.Gy.from.fluence.cm2}
\title{AT.dose.Gy.from.fluence.cm2}
\description{Returns dose in Gy for each given particle}
\usage{AT.dose.Gy.from.fluence.cm2(E.MeV.u, particle.no, fluence.cm2,
 material.no, stopping.power.source.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of
 size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array
 of size n) (see also \code{\link{particle.no}}).}
  \item{fluence.cm2}{ fluence for each particle type (array of size n).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ stopping power source index (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{dose.Gy}{ be allocated by the user which will be used to return the
 results (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Phy
sicsRoutines.c#L297}
}
\examples{
# Compute dose from protons, He-3, C-12 and O-16 at
# same energy and fluence in air
AT.dose.Gy.from.fluence.cm2( E.MeV.u      = c(10, 10, 10, 10),
                             fluence.cm2  = c(1e7, 1e7, 1e7, 1e7),
                             particle.no  = c(1001, 2003, 6012, 8016),
                             material.no  = 7,
                             stopping.power.source.no = 2)
}
