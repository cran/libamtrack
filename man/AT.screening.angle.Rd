% TODO File path/AT.screening.angle.Rd
\name{AT.screening.angle}
\alias{AT.screening.angle}
\title{AT.screening.angle}
\description{Returns screening angles chi_a
}
\usage{AT.screening.angle(E.MeV.u, particle.charge.e, element.acronym)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of size n).}
  \item{element.acronym}{ vector of elemental symbols of target material (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{chi.a}{ vector of screening angles in rad (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_MultipleCoulombScattering.c#L70}
}


