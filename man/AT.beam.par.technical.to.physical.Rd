% TODO File path/AT.beam.par.technical.to.physical.Rd
\name{AT.beam.par.technical.to.physical}
\alias{AT.beam.par.technical.to.physical}
\title{AT.beam.par.technical.to.physical}
\description{Converts technical, accelerator parameters of a symmetric, double lateral Gaussian shape beam, i.e.
 total number of particles and FWHM to
 physical beam parameters, i.e.
 central (=peak) fluence and width (= 1 standard deviation)
}
\usage{AT.beam.par.technical.to.physical(N, FWHM.mm)
}
\arguments{
  \item{N}{ absolute particle numbers (array of size n).}
  \item{FWHM.mm}{ FWHMs (in mm) (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{fluence.cm2}{ resulting fluence in beam center (array of size n)}
  \item{sigma.cm}{ resulting beam width stdev (array of size n)}
}
\seealso{
View the C source code here:

\url{ http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.4.1/src/AT_PhysicsRoutines.c#L320
}

}
\examples{
# Get peak dose of a 142.66 MeV protons in Alox
# from technical beam parameters
peak.fluence.cm2 <- AT.beam.par.technical.to.physical(  N        = 3.2e8, 
                                                        FWHM.mm  = 15.2)[1]
AT.dose.Gy.from.fluence.cm2( E.MeV.u     = 142.66,
                             particle.no = AT.particle.no.from.particle.name("1H"),
                             material.no = AT.material.no.from.material.name("Aluminum Oxide"),
                             fluence.cm2 = peak.fluence.cm2,
                             stopping.power.source.no = 0)
}


