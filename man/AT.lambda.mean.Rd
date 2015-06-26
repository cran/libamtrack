% TODO File path/AT.lambda.mean.Rd
\name{AT.lambda.mean}
\alias{AT.lambda.mean}
\title{AT.lambda.mean}
\description{Computes the mean lambda, introduced to enable
 average value for Landau distribution. See Geant3 W5013, p.254
}
\usage{AT.lambda.mean(E.MeV.u, particle.no, material.no, slab.thickness.um)
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
  \item{lambda.mean}{ mean lambda for given particle (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Ene
rgyLoss.c#L93}
}
