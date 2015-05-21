% TODO File path/AT.total.D.Gy.Rd
\name{AT.total.D.Gy}
\alias{AT.total.D.Gy}
\title{AT.total.D.Gy}
\description{Computes the total dose of a mixed particle field
}
\usage{AT.total.D.Gy(E.MeV.u, particle.no, fluence.cm2, material.no,
 stopping.power.source.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of
 size \code{\link{number.of.field.components}}) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (array of size
 \code{\link{number.of.field.components}}) (see also
 \code{\link{particle.no}}).}
  \item{fluence.cm2}{ fluences of particles in the mixed particle field (array
 of size \code{\link{number.of.field.components}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{total.dose.Gy}{total.dose.Gy}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Phy
sicsRoutines.c#L474}
}
\examples{
# Compute dose of monoenergtic high energy
# and a mixed (99/1) field of high and low 
# energy protons in water at same fluence
AT.total.D.Gy( E.MeV.u     = 100,
               particle.no = 1001, 
               fluence.cm2 = 100e8, 
               material.no = 1,
               stopping.power.source.no = 2)
AT.total.D.Gy( E.MeV.u     = c(100, 5),
               particle.no = c(1001, 1001), 
               fluence.cm2 = c(99e8, 1e8), 
               material.no = 1,
               stopping.power.source.no = 2)
}
