% TODO File path/AT.kappa.Rd
\name{AT.kappa}
\alias{AT.kappa}
\title{AT.kappa}
\description{Computes the kappa criterium for the
 energy loss distribution according to
 Seltzer and Berger, and CERN W5013
 No effective projectile charge is considered!}
\usage{AT.kappa(E.MeV.u, particle.no, material.no, slab.thickness.um)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per amu (array of size n) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{kappa}{ kappa parameter (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Ene
rgyLoss.c#L65}
}
