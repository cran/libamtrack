% TODO File path/AT.set.user.material.Rd
\name{AT.set.user.material}
\alias{AT.set.user.material}
\title{AT.set.user.material}
\description{Initializes user defined material. The material can then be used
 with material index number 0. !Be aware! that is
 definition is only valid during run-time. When the library is reloaded, the
 default settings are restored and the
 material is removed.}
\usage{AT.set.user.material(density.g.cm3, I.eV, average.A, average.Z)
}
\arguments{
  \item{density.g.cm3}{ physical density in g per cm3.}
  \item{I.eV}{ I value in eV.}
  \item{average.A}{ average mass number.}
  \item{average.Z}{ average atomic number.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{status}{ material defined successfully if zero}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Dat
aMaterial.c#L479}
}
