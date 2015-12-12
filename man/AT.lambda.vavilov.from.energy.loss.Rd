% TODO File path/AT.lambda.vavilov.from.energy.loss.Rd
\name{AT.lambda.vavilov.from.energy.loss}
\alias{AT.lambda.vavilov.from.energy.loss}
\title{AT.lambda.vavilov.from.energy.loss}
\description{Computes the lambda parameter for the
 Vavilov distribution acc. to CERN W5013
 No effective projectile charge is considered!}
\usage{AT.lambda.vavilov.from.energy.loss(energy.loss.keV, E.MeV.u,
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
  \item{lambda.vavilov}{ (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Ene
rgyLoss.c#L304}
}
