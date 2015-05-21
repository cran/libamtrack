% TODO File path/AT.total.fluence.cm2.Rd
\name{AT.total.fluence.cm2}
\alias{AT.total.fluence.cm2}
\title{AT.total.fluence.cm2}
\description{Computes the total fluence of a mixed particle field
}
\usage{AT.total.fluence.cm2(E.MeV.u, particle.no, D.Gy, material.no,
 stopping.power.source.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of
 size \code{\link{number.of.field.components}}) (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle index (array of size
 \code{\link{number.of.field.components}}) (see also
 \code{\link{particle.no}}).}
  \item{D.Gy}{ doses of particles in the mixed particle field (array of size
 \code{\link{number.of.field.components}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{total.fluence.cm}{total.fluence.cm}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.0/src/AT_Phy
sicsRoutines.c#L502}
}
\examples{
# Compute total fluence in a 
# and a mixed field of
# high and low energy protons
# that deliver the same dose to water
AT.total.fluence.cm2( E.MeV.u     = c(100, 5),
                      particle.no = c(1001, 1001), 
                      D.Gy        = c(1, 1), 
                      material.no = 1,
                      stopping.power.source.no = 2)
}
