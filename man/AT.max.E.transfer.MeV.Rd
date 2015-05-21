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
  \item{E.MeV.u}{ energies of particle per nucleon [MeV/u]; if positive, the
 computation will be relativistic; if negative, the classic formular will be
 used (array of size n) (see also \code{\link{E.MeV.u}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{max.E.transfer.MeV}{ maximal energies transferred (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Phy
sicsRoutines.c#L209}
}
\examples{

# Plot maximum energy transferred in a collision in the range from 1 to 1000
# MeV
# (this is independent of the ion type!) and compare the classical with
# the relativistic approach
E.MeV.u     <- 10^seq(0, 3, length.out = 50)
df          <- data.frame( E.MeV                  = E.MeV.u,
                           max.E.keV.classical    = AT.max.E.transfer.MeV(-1.0
 * E.MeV.u)$max.E.transfer.MeV * 1000,
                           max.E.keV.relativistic =
 AT.max.E.transfer.MeV(E.MeV.u)$max.E.transfer.MeV * 1000)
}
