% TODO File path/AT.run.IGK.method.Rd
\name{AT.run.IGK.method}
\alias{AT.run.IGK.method}
\title{AT.run.IGK.method}
\description{Computes HCP response and relative efficiency/RBE using Katz' Ion-Gamma-Kill approach
 according to Waligorski, 1988
}
\usage{AT.run.IGK.method(E.MeV.u, particle.no, fluence.cm2.or.dose.Gy, material.no, stopping.power.source.no, 
	rdd.model, rdd.parameters, er.model, gamma.model, gamma.parameters, 
	saturation.cross.section.factor, write.output)
}
\arguments{
  \item{E.MeV.u}{ particle energy for each component in the mixed particle field [MeV/u] (array of size \code{\link{number.of.field.components}}) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle type for each component in the mixed particle field (array of size \code{\link{number.of.field.components}}) (see also \code{\link{particle.no}}).}
  \item{fluence.cm2.or.dose.Gy}{ if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size \code{\link{number.of.field.components}}) (see also \code{\link{fluence.cm2.or.dose.Gy}}).}
  \item{material.no}{ index number for detector material (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ stopping power source number (PSTAR,...) (see also \code{\link{stopping.power.source.no}}).}
  \item{rdd.model}{ index number for chosen radial dose distribution (see also \code{\link{rdd.model}}).}
  \item{rdd.parameters}{ parameters for chosen radial dose distribution (array of size 4).}
  \item{er.model}{ index number for chosen electron-range model (see also \code{\link{er.model}}).}
  \item{gamma.model}{ index number for chosen gamma response.}
  \item{gamma.parameters}{ parameters for chosen gamma response (array of size 9).}
  \item{saturation.cross.section.factor}{ scaling factor for the saturation cross section.}
  \item{write.output}{ if true, a protocol is written to a file in the working directory.}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{relative.efficiency}{ particle response at dose D / gamma response at dose D}
  \item{S.HCP}{ absolute particle response}
  \item{S.gamma}{ absolute gamma response}
  \item{sI.cm2}{ resulting ion saturation cross section in cm2}
  \item{gamma.dose.Gy}{ dose contribution from gamma kills}
  \item{P.I}{ ion kill probability}
  \item{P.g}{ gamma kill probability}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_Algorithms_IGK.c#L34}
}
\examples{
# Compute the relative efficiency of an Alanine detector in a mixed
# carbon / proton field
AT.run.IGK.method( particle.no                          = 1001,                        # namely protons with
                   E.MeV.u                              = 10,                          # 10 MeV/u
                   fluence.cm2.or.dose.Gy               = c(-1.0),                     # delivering 1 Gy
                   material.no                          = 5,                           # i.e. Alanine
                   rdd.model                            = 4,                           # Katz parametrization of radial dose distribution with simplified extended targets
                   rdd.parameter                        = c(5e-8,1e-10),               # with 50 nm target size and 1e-10 dose minimum
                   er.model                             = 2,                           # Butts&Katz parametrization of track radius
                   gamma.model                          = 2,                           # Use general target/hit model but here...
                   gamma.parameters                     = c(1,500,1,1,0),              # ...as exponential saturation with characteristic dose 500 Gy
                   saturation.cross.section.factor      = 1.4,                         # factor to take 'brush' around track into account
                   write.output                         = TRUE,                        # write a log file
                   stopping.power.source.no             = 0)                      
}


