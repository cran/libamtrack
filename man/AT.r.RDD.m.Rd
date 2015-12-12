% TODO File path/AT.r.RDD.m.Rd
\name{AT.r.RDD.m}
\alias{AT.r.RDD.m}
\title{AT.r.RDD.m}
\description{Returns distance as a function of dose
}
\usage{AT.r.RDD.m(D.RDD.Gy, E.MeV.u, particle.no, material.no, rdd.model, 
	rdd.parameter, er.model, stopping.power.source.no)
}
\arguments{
  \item{D.RDD.Gy}{ dose [Gy] (array of size n).}
  \item{E.MeV.u}{ particle (ion) energy per nucleon [MeV/u] (see also
 \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle code number (see also
 \code{\link{particle.no}}).}
  \item{material.no}{ material code number (see also
 \code{\link{material.no}}).}
  \item{rdd.model}{ Radial Dose Distribution model code number (see also
 \code{\link{rdd.model}}).}
  \item{rdd.parameter}{ Radial Dose Distribution model parameters vector
 (array of size 4).}
  \item{er.model}{ delta electron range model code number (see also
 \code{\link{er.model}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{r.RDD.m}{ distance [m] (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_RDD
.c#L662}
}
