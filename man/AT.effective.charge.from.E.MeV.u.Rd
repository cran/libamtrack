% TODO File path/AT.effective.charge.from.E.MeV.u.Rd
\name{AT.effective.charge.from.E.MeV.u}
\alias{AT.effective.charge.from.E.MeV.u}
\title{AT.effective.charge.from.E.MeV.u}
\description{Effective charge according to Barkas-Bethe-approximation:
 for particles with given kinetic energy per nucleon
}
\usage{AT.effective.charge.from.E.MeV.u(E.MeV.u, particle.no)
}
\arguments{
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array of size n) (see also \code{\link{particle.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{effective.charge}{ Effective charge according to Barkas-Bethe-approximation (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.2/src/AT_PhysicsRoutines.c#L149}
}
\examples{
# Charge pick-up of several nuclids depending on ion energy
require(lattice)
df          <- data.frame( E.MeV.u          = 10^seq(-1, 2, length.out = 50),
                           particle.name    = c("1H", "3He", "6Li", "12C", "16O"),
                           effective.charge = 0)
for(i in 1:nrow(df)){
   df$effective.charge[i] <- AT.effective.charge.from.E.MeV.u(  E.MeV.u     = df$E.MeV.u[i], 
                                                                particle.no = AT.particle.no.from.particle.name(df$particle.name[i])[1])
}
xyplot( effective.charge ~ log10(E.MeV.u),
        df,
        type     = 'l',
        groups   = particle.name,
        auto.key = list( space = 'right', points = FALSE, lines = TRUE),
        xlab     = 'kinetic energy / (MeV/u)',
        ylab     = 'effective charge',
        panel    = function(...){
                       panel.grid(-1,-1)
                       panel.xyplot(...)},
        scales     = list( x = list( at        = -1:2,
                                      labels = c("0.1", "1", "10", "100"))))
}


