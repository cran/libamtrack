% TODO File path/AT.electron.density.m3.from.composition.Rd
\name{AT.electron.density.m3.from.composition}
\alias{AT.electron.density.m3.from.composition}
\title{AT.electron.density.m3.from.composition}
\description{Computes the electron density for a given material composition
}
\usage{AT.electron.density.m3.from.composition(density.g.cm3, Z, A, weight.fraction)
}
\arguments{
  \item{density.g.cm3}{ physical density (in g per cm3) of material.}
  \item{Z}{ atomic numbers of constituents (array of size n).}
  \item{A}{ mass numbers of constituents (array of size n).}
  \item{weight.fraction}{ relative fractions of weight of constituents (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{electron.density.m3}{ electron density per m3}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataMaterial.c#L364}
}


