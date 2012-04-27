% TODO File path/AT.characteristic.single.scattering.angle.Rd
\name{AT.characteristic.single.scattering.angle}
\alias{AT.characteristic.single.scattering.angle}
\title{AT.characteristic.single.scattering.angle}
\description{Returns characteristic single scattering angles chi_c
}
\usage{AT.characteristic.single.scattering.angle(E.MeV.u, particle.charge.e, target.thickness.cm, element.acronym)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of size n).}
  \item{target.thickness.cm}{ vector of thicknesses of target material in cm (array of size n).}
  \item{element.acronym}{ vector of elemental symbols of target material (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{chi.c}{ vector of characteristic single scattering angles in rad (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.3/src/AT_MultipleCoulombScattering.c#L34}
}


