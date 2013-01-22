% TODO File path/AT.reduced.target.thickness.Rd
\name{AT.reduced.target.thickness}
\alias{AT.reduced.target.thickness}
\title{AT.reduced.target.thickness}
\description{Returns reduced target thicknesses B
}
\usage{AT.reduced.target.thickness(E.MeV.u, particle.charge.e, target.thickness.cm, element.acronym)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of size n).}
  \item{target.thickness.cm}{ vector of thicknesses of target material in cm (array of size n).}
  \item{element.acronym}{ vector of elemental symbols of target material (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{B}{ vector of reduced target thicknesses (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_MultipleCoulombScattering.c#L147}
}


