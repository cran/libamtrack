% TODO File path/AT.WEPL.Bethe.Rd
\name{AT.WEPL.Bethe}
\alias{AT.WEPL.Bethe}
\title{AT.WEPL.Bethe}
\description{Computes the water equivalent path length (WEPL) using the Bethe
 formula}
\usage{AT.WEPL.Bethe(E.MeV.u, particle.no, material.no, slab.thickness.m)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per nucleon (array of size n) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.m}{ thickness of slab of material different than water,
 in meter.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{WEPL}{ resulting water equivalent path length (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Dat
aRange.c#L202}
}
