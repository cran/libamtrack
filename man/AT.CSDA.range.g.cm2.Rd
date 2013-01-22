% TODO File path/AT.CSDA.range.g.cm2.Rd
\name{AT.CSDA.range.g.cm2}
\alias{AT.CSDA.range.g.cm2}
\title{AT.CSDA.range.g.cm2}
\description{Returns CSDA range (in g/cm2) from pstar tables for given energy.
 In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)}
\usage{AT.CSDA.range.g.cm2(E.MeV.u, particle.no, material.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of size number.of.particles) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array of size number.of.particles) (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{CSDA.range.g.cm2}{ (array of size number.of.particles) to be allocated by the user which will be used to return the results}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataLET.c#L81}
}
\examples{
# Range of 270 MeV/u carbon ions and 142 MeV/u protons in water
AT.CSDA.range.g.cm2(    E.MeV.u               = c(270, 142),
                        particle.no           = AT.particle.no.from.particle.name(c("12C", "1H")),
                        material.no           = AT.material.no.from.material.name("Water, Liquid"))
}


