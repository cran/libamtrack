% TODO File path/AT.electron.density.m3.Rd
\name{AT.electron.density.m3}
\alias{AT.electron.density.m3}
\title{AT.electron.density.m3}
\description{Returns electron density from average A and Z}
\usage{AT.electron.density.m3(density.g.cm3, average.A, average.Z)
}
\arguments{
  \item{density.g.cm3}{ material density in g/cm3 (array of size n).}
  \item{average.A}{ average mass number (array of size n).}
  \item{average.Z}{ average atomic number (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{electron.density.m3}{ electron density in 1/m3 (array of size n)}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataMaterial.c#L350}
}


