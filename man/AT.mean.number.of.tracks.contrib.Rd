% TODO File path/AT.mean.number.of.tracks.contrib.Rd
\name{AT.mean.number.of.tracks.contrib}
\alias{AT.mean.number.of.tracks.contrib}
\title{AT.mean.number.of.tracks.contrib}
\description{Computes the number of track contributing to a representative
 point in a mixed field
}
\usage{AT.mean.number.of.tracks.contrib(E.MeV.u, particle.no, fluence.cm2,
 material.no, er.model, 
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
  \item{er.model}{ chosen electron-range-model (see also
 \code{\link{er.model}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{resulting}{resulting}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_PhysicsR
outines.c#L700}
}
