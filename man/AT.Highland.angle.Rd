% TODO File path/AT.Highland.angle.Rd
\name{AT.Highland.angle}
\alias{AT.Highland.angle}
\title{AT.Highland.angle}
\description{Returns Highland angles Theta0
}
\usage{AT.Highland.angle(E.MeV.u, particle.charge.e, l.over.lR)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of
 size n).}
  \item{l.over.lR}{ vector of thicknesses of target material in radiation
 lengths (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{Theta0}{ vector of Highland angles in rad (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_Multiple
CoulombScattering.c#L345}
}
