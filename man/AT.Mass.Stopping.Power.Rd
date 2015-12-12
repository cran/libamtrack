% TODO File path/AT.Mass.Stopping.Power.Rd
\name{AT.Mass.Stopping.Power}
\alias{AT.Mass.Stopping.Power}
\title{AT.Mass.Stopping.Power}
\description{Retrieves the electronic mass stopping power in MeV*cm2/g
 for the requested energies and particles for a specified
 material and data source. The data source is thereby
 given via its name (s. AT_StoppingPowerData.h from
 details), except for data that should be read for
 a file, in this case the (path and) filename has to be
 provided. In this case, the user has to make sure that
 energy and stopping power units are correct and that
 the data match the given material (use material.no = 0
 for custom-defined material).
 The file has to be plain
 ASCII with three columns (separated by space)
 charge, energy, and stopping power
 and sorted in ascending order by first charge than energy
 any alphanumeric comment can be inserted (in separate
 lines)
}
\usage{AT.Mass.Stopping.Power(stopping.power.source, E.MeV.u, particle.no,
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
  \item{stopping.power.MeV.cm2.g}{ array to return stopping powers (array of
 size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Sto
ppingPower.c#L39}
}
\examples{
# Use stopping power data from Bethe equation
AT.Mass.Stopping.Power( E.MeV.u               = 270.55, 
                        particle.no           = 6012, 
                        material.no           =
 AT.material.no.from.material.name("Water, Liquid"), 
                        stopping.power.source = "Bethe")
# Use stopping power data from PSTAR
AT.Mass.Stopping.Power( E.MeV.u               = 270.55, 
                        particle.no           = 6012, 
                        material.no           =
 AT.material.no.from.material.name("Water, Liquid"), 
                        stopping.power.source = "PSTAR")
# Use stopping power data from ICRU
AT.Mass.Stopping.Power( E.MeV.u               = 270.55, 
                        particle.no           = 6012, 
                        material.no           =
 AT.material.no.from.material.name("Water, Liquid"), 
                        stopping.power.source = "ICRU")
# Use stopping power data from file shipped with libamtrack
# For water:
AT.Mass.Stopping.Power( E.MeV.u               = 270.55, 
                        particle.no           = 6012, 
                        material.no           =
 AT.material.no.from.material.name("Water, Liquid"), 
                        stopping.power.source = system.file("extdata",
 "FLUKA_DEDX_WATER_76.8eV.txt", package = "libamtrack"))

# And alumina (CAVE: the data file contains stopping power in keV/um, which is
# return, despite the name of the function):
AT.Mass.Stopping.Power( E.MeV.u               = 270.55, 
                        particle.no           = 6012, 
                        material.no           =
 AT.material.no.from.material.name("Aluminum Oxide"), 
                        stopping.power.source = system.file("extdata",
 "FLUKA_DEDX_ALOX.txt", package = "libamtrack"))

# When using density scaling, the value becomes clearily wrong! Thus, do
# preferably use mass stopping power data in the files
AT.Stopping.Power( E.MeV.u               = 270.55, 
                   particle.no           = 6012, 
                   material.no           =
 AT.material.no.from.material.name("Aluminum Oxide"), 
                   stopping.power.source = system.file("extdata",
 "FLUKA_DEDX_ALOX.txt", package = "libamtrack"))
}
