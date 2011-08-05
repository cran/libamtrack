% TODO File path/AT.energy.straggling.MeV2.cm2.g.Rd
\name{AT.energy.straggling.MeV2.cm2.g}
\alias{AT.energy.straggling.MeV2.cm2.g}
\title{AT.energy.straggling.MeV2.cm2.g}
\description{Get energy spread with depth according to Bohr's classical theory
 Bohr, N. (1915), Phil. Mag. 30, 581ff, see also Evans, R.D. (1955), The atomic nucleus, McGraw Hill, New York, p. 661
 In the literature dsE2dz is often given in units ergs2/cm. Here we report it mass-normalized MeV2*cm2/g
 Since the effective charge of the particle enters the equation, particle types and energies have to be given
 The equation is however limited to energies > 10 MeV/u and not too heavy ions
 TODO: add William extension for relativistic effects (Williams, E.J. (1945), Revs. Mod. Phys. 17, 217ff)}
\usage{AT.energy.straggling.MeV2.cm2.g(E.MeV.u, particle.no, material.no)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array of size n) (see also \code{\link{particle.no}}).}
  \item{material.no}{ index number for slab material (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{dsE2dz.MeV2.cm2.g}{ Increase of energy straggling variance sigma.E^2 per unit length of material (array of size n)}
}


