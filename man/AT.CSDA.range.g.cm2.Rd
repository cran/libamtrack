% TODO File path/AT.CSDA.range.g.cm2.Rd
\name{AT.CSDA.range.g.cm2}
\alias{AT.CSDA.range.g.cm2}
\title{AT.CSDA.range.g.cm2}
\description{Computes the CSDA range using the PSTAR data}
\usage{AT.CSDA.range.g.cm2(E.initial.MeV.u, E.final.MeV.u, particle.no,
 material.no)
}
\arguments{
  \item{E.initial.MeV.u}{ initial energy of particle per nucleon (array of
 size n).}
  \item{E.final.MeV.u}{ final energy of particle per nucleon (array of size
 n).}
  \item{particle.no}{ particle index (array of size n) (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{CSDA.range.cm2.g}{ resulting range (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Dat
aRange.c#L101}
}
\examples{
# Range of 270 MeV/u carbon ions and 142 MeV/u protons in water
#AT.CSDA.range.g.cm2(    E.MeV.u               = c(270, 142),

#                        particle.no           =
# AT.particle.no.from.particle.name(c("12C", "1H")),

#                        material.no           =
# AT.material.no.from.material.name("Water, Liquid"))
}
