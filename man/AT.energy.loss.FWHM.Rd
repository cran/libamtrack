% TODO File path/AT.energy.loss.FWHM.Rd
\name{AT.energy.loss.FWHM}
\alias{AT.energy.loss.FWHM}
\title{AT.energy.loss.FWHM}
\description{Computes the width of the energy loss distribution
 Uses Landau, Vavilov or Gauss depending on
 kappa parameter
 No effective projectile charge is considered!}
\usage{AT.energy.loss.FWHM(E.MeV.u, particle.no, material.no,
 slab.thickness.um)
}
\arguments{
  \item{E.MeV.u}{ energy of particle per nucleon (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{slab.thickness.um}{ slab thickness in um.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{result}{result}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/trunk/src/AT_EnergyLo
ss.c#L683}
}
