% TODO File path/AT.run.CPPSC.method.Rd
\name{AT.run.CPPSC.method}
\alias{AT.run.CPPSC.method}
\title{AT.run.CPPSC.method}
\description{Computes HCP response and relative efficiency/RBE using compound
 Poison processes and
 successive convolutions (CPP_SC, the SPIFF algorithm)
}
\usage{AT.run.CPPSC.method(E.MeV.u, particle.no, fluence.cm2.or.dose.Gy,
 material.no, stopping.power.source.no, 
	rdd.model, rdd.parameters, er.model, gamma.model, gamma.parameters, 
	N2, fluence.factor, write.output, shrink.tails, shrink.tails.under, 
	adjust.N2, lethal.events.mode)
}
\arguments{
  \item{E.MeV.u}{ particle energy for each component in the mixed particle
 field [MeV/u] (array of size \code{\link{number.of.field.components}}) (see
 also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle type for each component in the mixed particle
 field (array of size \code{\link{number.of.field.components}}) (see also
 \code{\link{particle.no}}).}
  \item{fluence.cm2.or.dose.Gy}{ if positive, particle fluence for each
 component in the mixed particle field [1/cm2]; if negative, particle dose for
 each component in the mixed particle field [Gy] (array of size
 \code{\link{number.of.field.components}}) (see also
 \code{\link{fluence.cm2.or.dose.Gy}}).}
  \item{material.no}{ index number for detector material (see also
 \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
  \item{rdd.model}{ index number for chosen radial dose distribution (see also
 \code{\link{rdd.model}}).}
  \item{rdd.parameters}{ parameters for chosen radial dose distribution (array
 of size 4).}
  \item{er.model}{ index number for chosen electron-range model (see also
 \code{\link{er.model}}).}
  \item{gamma.model}{ index number for chosen gamma response.}
  \item{gamma.parameters}{ parameters for chosen gamma response (array of size
 9).}
  \item{N2}{ number of bins per factor of two for the dose scale of local dose
 histogram.}
  \item{fluence.factor}{ factor to scale the fluences / doses given in
 \code{\link{fluence.cm2.or.dose.Gy}} with.}
  \item{write.output}{ if true, a log-file is written to
 SuccessiveConvolutions.txt in the working directory.}
  \item{shrink.tails}{ if true, tails of the local dose distribution,
 contributing less than shrink.tails.under are cut.}
  \item{shrink.tails.under}{ limit for tail cutting in local dose
 distribution.}
  \item{adjust.N2}{ if true, N2 will be increase if necessary at high fluence
 to ensure sufficient local dose histogram resolution.}
  \item{lethal.events.mode}{ if true, computations are done for dependent
 subtargets.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{N2}{ number of bins per factor of two for the dose scale of local dose
 histogram}
  \item{relative.efficiency}{ particle response at dose D / gamma response at
 dose D}
  \item{d.check}{ sanity check: total dose (in Gy) as returned by the
 algorithm}
  \item{S.HCP}{ absolute particle response}
  \item{S.gamma}{ absolute gamma response}
  \item{mean.number.of.tracks.contrib}{ mean number of tracks contributing to
 representative point}
  \item{start.number.of.tracks.contrib}{ low fluence approximation for mean
 number of tracks contributing to representative point (start value for
 successive convolutions)}
  \item{n.convolutions}{ number of convolutions performed to reach requested
 dose/fluence}
  \item{lower.Jensen.bound}{ lower bound for Jensen's inequity}
  \item{upper.Jensen.bound}{ upper bound for Jensen's inequity}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Alg
orithms_CPP.c#L34}
}
\examples{
# Compute the relative efficiency of an Alanine detector in a mixed
# carbon / proton field
AT.run.CPPSC.method( particle.no                          = c(6012, 1001,
 1001),         # namely carbon, protons, and protons with
                     E.MeV.u                              = c(270, 270, 5),   
           
# 270 MeV/u (primary Carbon, 270 MeV/u and 5 MeV/u (fast and slow proton
# component) 
                     fluence.cm2.or.dose.Gy               = c(1e8, 1e9, 1e7), 
           # and their corresponding fluences
                     material.no                          = 5,                
           # i.e. Alanine
                     rdd.model                            = 3,                
           # simple 'Geiss' parametrization of radial dose distribution
                     rdd.parameter                        = 50e-9,            
           # with 50 nm core radius
                     er.model                             = 4,                
           # M. Scholz' parametrization of track radius
                     gamma.model                          = 2,                
           # General hit/target X ray response, but
                     gamma.parameters                     = c(1,500,1,1,0),   
           
# as simple single exponential saturation (one hit, one target), saturation
# dose 500 Gy
                     N2                                   = 10,               
           # ten bins per factor 2 for internal local dose histogramming
                     fluence.factor                       = 1.0,              
           # can be used to easily scale total fluence (historical)
                     write.output                         = TRUE,             
           # write a log file
                     shrink.tails                         = TRUE,             
           # cut tails of local dose distribution, if...
                     shrink.tails.under                   = 1e-30,            
           # ... they contribute less then 1e-30 to first moment of histogram
                     adjust.N2                            = TRUE,             
           # perform rebinning if local dose distribution becomes too narrow
                     lethal.events.mode                   = FALSE,            
           # use independent subtargets
                     stopping.power.source.no             = 2)                
      
}
