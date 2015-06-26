% TODO File path/AT.beam.par.physical.to.technical.Rd
\name{AT.beam.par.physical.to.technical}
\alias{AT.beam.par.physical.to.technical}
\title{AT.beam.par.physical.to.technical}
\description{Converts physical beam parameters of a symmetric, double lateral
 Gaussian shape beam, i.e.
 central (=peak) fluence and width (= 1 standard deviation)
 to technical, accelerator parameters, i.e.
 total number of particles and FWHM
}
\usage{AT.beam.par.physical.to.technical(fluence.cm2, sigma.cm)
}
\arguments{
  \item{fluence.cm2}{ fluence in beam center (array of size n).}
  \item{sigma.cm}{ beam width stdev (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{N}{ resulting absolute particle numbers (array of size n)}
  \item{FWHM.mm}{ resulting FWHMs (in mm) (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Phy
sicsRoutines.c#L430}
}
\examples{
# Technical parameters for a double Gaussian beam
# of Carbon ions at 270 MeV/u in water with 
# 0.5 cm sigma width and 1 Gy in the peak
AT.beam.par.physical.to.technical(  fluence.cm2 = AT.fluence.cm2.from.dose.Gy(
 E.MeV.u      = 270, 
                                                                              
 D.Gy         = 1.0, 
                                                                              
 particle.no  = AT.particle.no.from.particle.name("12C"), 
                                                                              
 material.no  = AT.material.no.from.material.name("Water, Liquid"),
                                                                              
 stopping.power.source.no = 2), 
                                    sigma.cm    = 0.5)
}
