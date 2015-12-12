% TODO File path/AT.xi.keV.Rd
\name{AT.xi.keV}
\alias{AT.xi.keV}
\title{AT.xi.keV}
\description{Parameter xi - reduced mean energy loss}
\usage{AT.xi.keV(E.MeV.u, particle.no, material.no, slab.thickness.um)
}
\arguments{
  \item{E.MeV.u}{ energies of particle per nucleon (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle indices (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{xi}{xi}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Ene
rgyLoss.c#L26}
}
