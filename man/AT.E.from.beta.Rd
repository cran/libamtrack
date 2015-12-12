% TODO File path/AT.E.from.beta.Rd
\name{AT.E.from.beta}
\alias{AT.E.from.beta}
\title{AT.E.from.beta}
\description{Returns energy per nucleon of particle with relativistic speed
 beta
}
\usage{AT.E.from.beta(beta)
}
\arguments{
  \item{beta}{ vector of relative particle speed beta = v/c (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{E.MeV.u}{ vector of energies of particle per nucleon [MeV] (array of
 size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Phy
sicsRoutines.c#L71}
}
\examples{
# Compute energies for betas between 0.1 and 0.99:
AT.E.from.beta( beta = c(0.1, 0.1*(1:9), 0.99))
}
