% TODO File path/AT.dose.weighted.E.MeV.u.Rd
\name{AT.dose.weighted.E.MeV.u}
\alias{AT.dose.weighted.E.MeV.u}
\title{AT.dose.weighted.E.MeV.u}
\description{Computes the dose-weighted average energy of a particle field
 Needed by SuccessiveConvolutions
}
\usage{AT.dose.weighted.E.MeV.u(E.MeV.u, particle.no, fluence.cm2,
 material.no, stopping.power.source.no)
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
  \item{dose-weighted}{dose-weighted}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Phy
sicsRoutines.c#L553}
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
                                   stopping.power.source.no = 2)   
}
