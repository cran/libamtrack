% TODO File path/AT.Vavilov.energy.loss.distribution.Rd
\name{AT.Vavilov.energy.loss.distribution}
\alias{AT.Vavilov.energy.loss.distribution}
\title{AT.Vavilov.energy.loss.distribution}
\description{Computes the energy loss
 Vavilov distribution acc. to Seltzer & Berger
 No effective projectile charge is considered!}
\usage{AT.Vavilov.energy.loss.distribution(energy.loss.keV, E.MeV.u,
 particle.no, material.no, slab.thickness.um)
}
\arguments{
  \item{energy.loss.keV}{ energy loss (array of size n).}
  \item{E.MeV.u}{ energy of particle per nucleon (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{fDdD}{ (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Ene
rgyLoss.c#L403}
}
