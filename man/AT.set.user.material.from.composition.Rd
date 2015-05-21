% TODO File path/AT.set.user.material.from.composition.Rd
\name{AT.set.user.material.from.composition}
\alias{AT.set.user.material.from.composition}
\title{AT.set.user.material.from.composition}
\description{Initializes user defined material from composition data. The
 material can then be used with material index number 0. !Be aware! that is
 definition is only valid during run-time. When the library is reloaded, the
 default settings are restored and the
 material is removed.
}
\usage{AT.set.user.material.from.composition(density.g.cm3, Z, A,
 weight.fraction)
}
\arguments{
  \item{density.g.cm3}{ physical density (in g per cm3) of material.}
  \item{Z}{ atomic numbers of constituents (array of size n).}
  \item{A}{ mass numbers of constituents (array of size n).}
  \item{weight.fraction}{ relative fractions of weight of constituents (array
 of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{status}{ material defined successfully if zero}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Dat
aMaterial.c#L494}
}
