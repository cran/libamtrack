% TODO File path/AT.particle.no.from.Z.and.A.Rd
\name{AT.particle.no.from.Z.and.A}
\alias{AT.particle.no.from.Z.and.A}
\title{AT.particle.no.from.Z.and.A}
\description{Returns particle index number from given A and Z
}
\usage{AT.particle.no.from.Z.and.A(Z, A)
}
\arguments{
  \item{Z}{ atomic numbers (array of size n).}
  \item{A}{ mass number (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{particle.no}{ corresponding particle index numbers (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Dat
aParticle.c#L34}
}
