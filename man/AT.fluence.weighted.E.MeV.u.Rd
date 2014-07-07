% TODO File path/AT.fluence.weighted.E.MeV.u.Rd
\name{AT.fluence.weighted.E.MeV.u}
\alias{AT.fluence.weighted.E.MeV.u}
\title{AT.fluence.weighted.E.MeV.u}
\description{Computes the fluence-weighted average energy of a particle field
 Needed by SuccessiveConvolutions
}
\usage{AT.fluence.weighted.E.MeV.u(E.MeV.u, fluence.cm2)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of
 size \code{\link{number.of.field.components}}) (see also
 \code{\link{E.MeV.u}}).}
  \item{fluence.cm2}{ fluences of particles in the mixed particle field (array
 of size \code{\link{number.of.field.components}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{average.E.MeV.u}{average.E.MeV.u}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_PhysicsR
outines.c#L514}
}
\examples{
# fluence- and dose-weighted mean energy for a simple mixed field
# of high and low (99/1) energy protons
AT.fluence.weighted.E.MeV.u( E.MeV.u       = c(100, 1),
                             fluence.cm2   = c(99e8, 1e8)) 
AT.dose.weighted.E.MeV.u( E.MeV.u       = c(100, 1),
                          particle.no   = c(1001, 1001),
                          fluence.cm2   = c(99e8, 1e8), 
                          material.no   = 1,                   # water
                                   stopping.power.source.no = 0)    
}
