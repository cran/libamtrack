% TODO File path/AT.effective.collision.number.Rd
\name{AT.effective.collision.number}
\alias{AT.effective.collision.number}
\title{AT.effective.collision.number}
\description{Returns effective collision numbers exp_b
}
\usage{AT.effective.collision.number(E.MeV.u, particle.charge.e, target.thickness.cm, element.acronym)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.charge.e}{ vector of charge numbers of particle (array of size n).}
  \item{target.thickness.cm}{ vector of thicknesses of target material in cm (array of size n).}
  \item{element.acronym}{ vector of elemental symbols of target material (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{exp.b}{ vector of effective collision numbers (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.3/src/AT_MultipleCoulombScattering.c#L101}
}


