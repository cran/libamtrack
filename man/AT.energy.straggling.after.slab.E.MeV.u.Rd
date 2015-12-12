% TODO File path/AT.energy.straggling.after.slab.E.MeV.u.Rd
\name{AT.energy.straggling.after.slab.E.MeV.u}
\alias{AT.energy.straggling.after.slab.E.MeV.u}
\title{AT.energy.straggling.after.slab.E.MeV.u}
\description{Get energy spread of an ion beam after traversing
 a material slab according to Bohr's classical theory.
 Bohr, N. (1915), Phil. Mag. 30, 581ff, see also Evans, R.D. (1955), The
 atomic nucleus, McGraw Hill, New York, p. 661
 Please note that the effective charge is assumed to be constant over the
 material slab
 If this is not the case you should apply this routine multiple times to
 subslices}
\usage{AT.energy.straggling.after.slab.E.MeV.u(E.MeV.u, particle.no,
 material.no, slab.thickness.m, initial.sigma.E.MeV.u)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array
 of size n) (see also \code{\link{particle.no}}).}
  \item{material.no}{ index number for slab material (see also
 \code{\link{material.no}}).}
  \item{slab.thickness.m}{ thickness of slab in m.}
  \item{initial.sigma.E.MeV.u}{ energy spread - 1 sigma - before traversing
 the slab - can be 0 (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{sigma.E.MeV.u}{ energy spread - 1 sigma - after traversing the slab
 (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Phy
sicsRoutines.c#L318}
}
