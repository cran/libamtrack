% TODO File path/AT.Stopping.Power.Rd
\name{AT.Stopping.Power}
\alias{AT.Stopping.Power}
\title{AT.Stopping.Power}
\description{Retrieves the electronic stopping power in keV/um for
 the requested energies and particles for a specified
 material and data source. The data source is thereby
 given via its name (s. AT_StoppingPowerData.h for
 details), except for data that should be read from
 a file, in this case the (path and) filename has to be
 provided. In this case, the user has to make sure that
 energy and stopping power units are correct and that
 the data match the given material (use material.no = 0
 for custom-defined material) for density scaling.
 The file has to be plain
 ASCII with three columns (separated by space)
 charge, energy, and stopping power
 and sorted in ascending order by first charge than energy
 any alphanumeric comment can be inserted (in separate
 lines)
}
\usage{AT.Stopping.Power(stopping.power.source, E.MeV.u, particle.no,
 material.no)
}
\arguments{
  \item{stopping.power.source}{ name of the data source (see also
 \code{\link{stopping.power.source}}).}
  \item{E.MeV.u}{ kinetic energies in MeV per amu (array of size n) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle numbers (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material number (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{stopping.power.keV.um}{ array to return stopping powers (array of size
 n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Sto
ppingPower.c#L86}
}
