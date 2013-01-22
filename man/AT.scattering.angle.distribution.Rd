% TODO File path/AT.scattering.angle.distribution.Rd
\name{AT.scattering.angle.distribution}
\alias{AT.scattering.angle.distribution}
\title{AT.scattering.angle.distribution}
\description{Returns scattering angle distribution f(Theta)
 The distribution is not normalized because the energy loss in the target is not considered.
}
\usage{AT.scattering.angle.distribution(E.MeV.u, particle.charge.e, target.thickness.cm, element.acronym, Theta)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per nucleon [MeV] (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ charge number of particle.}
  \item{target.thickness.cm}{ thickness of the target material in cm.}
  \item{element.acronym}{ elemental symbol of target material (array of size PARTICLE.NAME.NCHAR).}
  \item{Theta}{ vector of polar angles in rad (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{distribution}{ vector of scattering angle distribution values (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_MultipleCoulombScattering.c#L285}
}


