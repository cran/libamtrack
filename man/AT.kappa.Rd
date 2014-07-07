% TODO File path/AT.kappa.Rd
\name{AT.kappa}
\alias{AT.kappa}
\title{AT.kappa}
\description{Computes the kappa criterium for the
 energy loss distribution according to
 Seltzer & Berger
 No effective projectile charge is considered!}
\usage{AT.kappa(E.MeV.u, particle.no, material.no, slab.thickness.um)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per nucleon (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{result}{result}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_EnergyLo
ss.c#L294}
}
