% TODO File path/AT.effective.Z.from.composition.Rd
\name{AT.effective.Z.from.composition}
\alias{AT.effective.Z.from.composition}
\title{AT.effective.Z.from.composition}
\description{Computes the effective atomic number for a given material composition
}
\usage{AT.effective.Z.from.composition(Z, weight.fraction, electron.densities.cm3, exponent)
}
\arguments{
  \item{Z}{ atomic numbers of constituents (array of size n).}
  \item{weight.fraction}{ relative fractions of weight of constituents (array of size n).}
  \item{electron.densities.cm3}{ if not zero, weight fractions will additionally include electron densities per volume (array of size n).}
  \item{exponent}{ exponent for additivity rule reflecting the photon energy regime (usually 3.5 at ~ 100 kV).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{effective.Z}{ effective Z}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.2/src/AT_DataMaterial.c#L420}
}


