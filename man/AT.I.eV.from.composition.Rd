% TODO File path/AT.I.eV.from.composition.Rd
\name{AT.I.eV.from.composition}
\alias{AT.I.eV.from.composition}
\title{AT.I.eV.from.composition}
\description{Computes the I value for a given material composition
}
\usage{AT.I.eV.from.composition(Z, A, weight.fraction)
}
\arguments{
  \item{Z}{ atomic numbers of constituents (array of size n).}
  \item{A}{ mass numbers of constituents (array of size n).}
  \item{weight.fraction}{ relative fractions of weight of constituents (array
 of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{I.eV}{ I value in eV}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Dat
aMaterial.c#L452}
}
