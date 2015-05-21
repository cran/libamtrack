% TODO File path/AT.CSDA.energy.after.slab.E.MeV.u.Rd
\name{AT.CSDA.energy.after.slab.E.MeV.u}
\alias{AT.CSDA.energy.after.slab.E.MeV.u}
\title{AT.CSDA.energy.after.slab.E.MeV.u}
\description{Computes the ion energy after transversing a slab of material
 using Bethe stopping power
 and CSDA approach for many energies / particles
}
\usage{AT.CSDA.energy.after.slab.E.MeV.u(E.initial.MeV.u, particle.no,
 material.no, slab.thickness.m)
}
\arguments{
  \item{E.initial.MeV.u}{ initial energy of particle per nucleon (array of
 size n).}
  \item{particle.no}{ particle index (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.m}{ thickness of slab to transversed.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{E.final.MeV.u}{ final energy after slab (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Dat
aRange.c#L160}
}
