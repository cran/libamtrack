% TODO File path/AT.get.materials.data.Rd
\name{AT.get.materials.data}
\alias{AT.get.materials.data}
\title{AT.get.materials.data}
\description{Returns material data for list of materials}
\usage{AT.get.materials.data(material.no)
}
\arguments{
  \item{material.no}{ material indices (array of size number.of.materials) (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{density.g.cm3}{ material density in g/cm3 (array of size number.of.materials)}
  \item{I.eV}{ mean ionization potential in eV (array of size number.of.materials)}
  \item{alpha.g.cm2.MeV}{ fit parameter for power-law representation of stp.power/range/E-dependence (array of size number.of.materials)}
  \item{p.MeV}{ fit parameter for power-law representation of stp.power/range/E-dependence (array of size number.of.materials)}
  \item{m.g.cm2}{ fit parameter for the linear representation of fluence changes due to nuclear interactions based on data from Janni 1982 (array of size number.of.materials)}
  \item{average.A}{ average mass number (array of size number.of.materials)}
  \item{average.Z}{ average atomic number (array of size number.of.materials)}
}


