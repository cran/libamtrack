% TODO File path/AT.max.E.transfer.MeV.new.Rd
\name{AT.max.E.transfer.MeV.new}
\alias{AT.max.E.transfer.MeV.new}
\title{AT.max.E.transfer.MeV.new}
\description{Kinetic energy maximally transferred from an ion to an electron
 in a collision - relativistic or non-relativistic
}
\usage{AT.max.E.transfer.MeV.new(E.MeV.u, A)
}
\arguments{
  \item{E.MeV.u}{ energies of particle per nucleon [MeV/u]; if positive, the
 computation will be relativistic; if negative, the classic formular will be
 used (array of size n) (see also \code{\link{E.MeV.u}}).}
  \item{A}{ atomic mass (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{max.E.transfer.MeV}{ maximal energies transferred (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Phy
sicsRoutines.c#L209}
}
