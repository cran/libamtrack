% TODO File path/AT.run.GSM.method.Rd
\name{AT.run.GSM.method}
\alias{AT.run.GSM.method}
\title{AT.run.GSM.method}
\description{Computes HCP response and relative efficiency/RBE using summation of tracks
 an a Cartesian grid (the GSM algorithm).

 Be aware that this routine can take considerable time to compute depending on
 the arguments, esp. for higher energy (>10 MeV/u) particles. It is therefore
 advantageous to test your settings with a low number of runs first.
}
\usage{AT.run.GSM.method(E.MeV.u, particle.no, fluence.cm2.or.dose.Gy, material.no, stopping.power.source.no, 
	rdd.model, rdd.parameters, er.model, gamma.model, gamma.parameters, 
	N.runs, write.output, nX, voxel.size.m, lethal.events.mode)
}
\arguments{
  \item{E.MeV.u}{ particle energy for each component in the mixed particle field [MeV/u] (array of size \code{\link{number.of.field.components}}) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle type for each component in the mixed particle field (array of size \code{\link{number.of.field.components}}) (see also \code{\link{particle.no}}).}
  \item{fluence.cm2.or.dose.Gy}{ if positive, particle fluence for each component in the mixed particle field [1/cm2]; if negative, particle dose for each component in the mixed particle field [Gy] (array of size \code{\link{number.of.field.components}}) (see also \code{\link{fluence.cm2.or.dose.Gy}}).}
  \item{material.no}{ index number for detector material (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ TODO (see also \code{\link{stopping.power.source.no}}).}
  \item{rdd.model}{ index number for chosen radial dose distribution (see also \code{\link{rdd.model}}).}
  \item{rdd.parameters}{ parameters for chosen radial dose distribution (array of size 4).}
  \item{er.model}{ index number for chosen electron-range model (see also \code{\link{er.model}}).}
  \item{gamma.model}{ index number for chosen gamma response.}
  \item{gamma.parameters}{ parameters for chosen gamma response (array of size 9).}
  \item{N.runs}{ number of runs within which track positions will be resampled.}
  \item{write.output}{ if true, a protocol is written to SuccessiveConvolutions.txt in the working directory.}
  \item{nX}{ number of voxels of the grid in x (and y as the grid is quadratic).}
  \item{voxel.size.m}{ side length of a voxel in m.}
  \item{lethal.events.mode}{ if true, allows to do calculations for cell survival.}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{relative.efficiency}{ particle response at dose D / gamma response at dose D}
  \item{d.check}{ sanity check: total dose (in Gy) as returned by the algorithm}
  \item{S.HCP}{ absolute particle response}
  \item{S.gamma}{ absolute gamma response}
  \item{n.particles}{ average number of particle tracks on the detector grid}
  \item{sd.relative.efficiency}{ standard deviation for relative.efficiency}
  \item{sd.d.check}{ standard deviation for d.check}
  \item{sd.S.HCP}{ standard deviation for S.HCP}
  \item{sd.S.gamma}{ standard deviation for S.gamma}
  \item{sd.n.particles}{ standard deviation for n.particles}
}
\examples{
# Compute the relative efficiency of an Alanine detector
#  in a proton field
AT.run.GSM.method( particle.no                          = 1001,                        # protons with
                   E.MeV.u                              = 10,                          # 10 MeV/u
                   fluence.cm2.or.dose.Gy               = c(-1.0),                     # delivering 1 Gy
                   material.no                          = 5,                           # i.e. Alanine
                   rdd.model                            = 3,                           # simple 'Geiss' parametrization of radial dose distribution
                   rdd.parameter                        = 50e-9,                       # with 50 nm core radius
                   er.model                             = 4,                           # M. Scholz' parametrization of track radius
                   gamma.model                          = 4,                           # Use exponential saturation
                   gamma.parameters                     = c(1,500),                    # max. response normalized to 1, saturation dose 500 Gy
                   N.runs                               = 1000,                        # resample 1000 times
                   write.output                         = TRUE,                        # write a log file
                   nX                                   = 10,                          # use a 10x10 grid
                   voxel.size.m                         = 5e-9,                        # with 5 nm voxel size
                   lethal.events.mode                   = FALSE,                       # use independent subtargets
                   stopping.power.source.no             = 0)                      
}


