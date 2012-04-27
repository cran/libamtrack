% TODO File path/AT.E.MeV.u.from.momentum.MeV.c.u.Rd
\name{AT.E.MeV.u.from.momentum.MeV.c.u}
\alias{AT.E.MeV.u.from.momentum.MeV.c.u}
\title{AT.E.MeV.u.from.momentum.MeV.c.u}
\description{Returns energy per nucleon for particles with given momentum per nucleon
}
\usage{AT.E.MeV.u.from.momentum.MeV.c.u(momentum.MeV.c.u)
}
\arguments{
  \item{momentum.MeV.c.u}{ vector of particle momenta per nucleon [MeV/c], (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV], (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.3/src/AT_PhysicsRoutines.c#L112}
}
\examples{
# Get kinetic energy for a 502 MeV/c antiproton beam (CERN AD/ACE experiment):
AT.E.MeV.u.from.momentum.MeV.c.u( momentum.MeV.c.u = 502)
}


