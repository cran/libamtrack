% TODO File path/AT.average.Z.from.composition.Rd
\name{AT.average.Z.from.composition}
\alias{AT.average.Z.from.composition}
\title{AT.average.Z.from.composition}
\description{Computes the average atomic number for a given material composition
}
\usage{AT.average.Z.from.composition(Z, weight.fraction)
}
\arguments{
  \item{Z}{ atomic numbers of constituents (array of size n).}
  \item{weight.fraction}{ relative fractions of weight of constituents (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{average.Z}{ average Z}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.2/src/AT_DataMaterial.c#L403}
}


