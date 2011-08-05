% TODO File path/AT.max.E.transfer.MeV.Rd
\name{AT.max.E.transfer.MeV}
\alias{AT.max.E.transfer.MeV}
\title{AT.max.E.transfer.MeV}
\description{Kinetic energy maximally transferred from an ion to an electron
 in a collision - relativistic or non-relativistic
}
\usage{AT.max.E.transfer.MeV(E.MeV.u)
}
\arguments{
  \item{E.MeV.u}{ energies of particle per nucleon [MeV/u]; if positive, the computation will be relativistic; if negative, the classic formular will be used (array of size n) (see also \code{\link{E.MeV.u}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{max.E.transfer.MeV}{ maximal energies transferred (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{ http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.4.1/src/AT_PhysicsRoutines.c#L167
}

}
\examples{
# Plot maximum energy transferred in a collision in the range from 1 to 1000 MeV
# (this is independent of the ion type!) and compare the classical with
# the relativistic approach
require(lattice)
E.MeV.u     <- 10^seq(0, 3, length.out = 50)
df          <- data.frame( E.MeV                  = E.MeV.u,
                           max.E.keV.classical    = AT.max.E.transfer.MeV(-1.0 * E.MeV.u)$max.E.transfer.MeV * 1000,
                           max.E.keV.relativistic = AT.max.E.transfer.MeV(E.MeV.u)$max.E.transfer.MeV * 1000)
plot1 <- xyplot( log10(max.E.keV.relativistic) ~ log10(E.MeV),
                 df,
                 type     = 'o',
                 xlab     = list('ion kinetic energy', font = 4),
                 ylab     = list('maximum electron energy (relativistic)', font = 4),
                 scales   = list( x = list( at = 0:3, labels = c("1 MeV/u", "10 MeV/u", "100 MeV/u", "1 GeV/u")),
                                  y = list( at = 1:3, labels = c("10 keV", "100 keV", "1 MeV"))))
plot2 <- xyplot( log10(max.E.keV.classical) ~ log10(max.E.keV.relativistic),
                 df,
                 type     = 'l',
                 col      = 'red',
                 lwd      = 2,
                 xlab     = 'maximum electron energy (relativistic) / keV',
                 ylab     = 'maximum electron energy (classical) / keV',
                 xlim     = c(0, 4),
                 ylim     = c(0, 4),
                 scales   = list( at = 1:3, labels = c("10", "100", "1000")),
                 panel    = function(...){
                                panel.grid(-1, -1)
                                panel.abline(a = 0, b = 1, lty = 2, col = 'grey')
                                panel.xyplot(...)})
plot(plot1, split = c(1, 1, 2, 1))
plot(plot2, split = c(2, 1, 2, 1), newpage = FALSE)
}


