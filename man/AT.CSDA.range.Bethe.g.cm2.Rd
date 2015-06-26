% TODO File path/AT.CSDA.range.Bethe.g.cm2.Rd
\name{AT.CSDA.range.Bethe.g.cm2}
\alias{AT.CSDA.range.Bethe.g.cm2}
\title{AT.CSDA.range.Bethe.g.cm2}
\description{Computes the CSDA range using the Bethe formula
 (AT_Stopping_Power_Bethe_Number)
 according to ICRU49, p.6, Eq. 2.1
 BUT WITHOUT shell or density, Bloch or Barkas correction!}
\usage{AT.CSDA.range.Bethe.g.cm2(E.initial.MeV.u, E.final.MeV.u, particle.no,
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
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Dat
aRange.c#L100}
}
