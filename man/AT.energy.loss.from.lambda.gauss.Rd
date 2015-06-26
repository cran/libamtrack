% TODO File path/AT.energy.loss.from.lambda.gauss.Rd
\name{AT.energy.loss.from.lambda.gauss}
\alias{AT.energy.loss.from.lambda.gauss}
\title{AT.energy.loss.from.lambda.gauss}
\description{Computes the energy loss from the lambda parameter of the
 Gauss distribution for compatibility with CERN W5013
 No effective projectile charge is considered!}
\usage{AT.energy.loss.from.lambda.gauss(lambda.gauss, E.MeV.u, particle.no,
 material.no, slab.thickness.um)
}
\arguments{
  \item{lambda.gauss}{ Gauss lambda (array of size n).}
  \item{E.MeV.u}{ energy of particle per nucleon (array of size n) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{energy.loss.keV}{ (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Ene
rgyLoss.c#L607}
}
