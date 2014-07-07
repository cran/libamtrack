% TODO File path/AT.characteristicple.scattering.angle.Rd
\name{AT.characteristicple.scattering.angle}
\alias{AT.characteristicple.scattering.angle}
\title{AT.characteristicple.scattering.angle}
\description{Returns characteristic multiple scattering angles Theta_M
}
\usage{AT.characteristicple.scattering.angle(E.MeV.u, particle.charge.e,
 target.thickness.cm, element.acronym)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of
 size n).}
  \item{target.thickness.cm}{ vector of thicknesses of target material in cm
 (array of size n).}
  \item{element.acronym}{ vector of elemental symbols of target material
 (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{Theta.M}{ vector of characteristic multiple scattering angles in rad
 (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_Multiple
CoulombScattering.c#L181}
}
