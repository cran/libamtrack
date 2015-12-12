% TODO File path/AT.translate.dose.into.DSB.distribution.Rd
\name{AT.translate.dose.into.DSB.distribution}
\alias{AT.translate.dose.into.DSB.distribution}
\title{AT.translate.dose.into.DSB.distribution}
\description{Converts a local dose into a DSB distribution assuming
 Poissonian rule for creation.
}
\usage{AT.translate.dose.into.DSB.distribution(f.d.Gy, f.dd.Gy, f,
 enhancement.factor, DSB.per.Gy.per.domain, 
	domains.per.nucleus, write.output)
}
\arguments{
  \item{f.d.Gy}{ bin midpoints for f (array of size n.bins.f).}
  \item{f.dd.Gy}{ bin widths for f (array of size n.bins.f).}
  \item{f}{ dose frequency (array of size n.bins.f).}
  \item{enhancement.factor}{ dose enhancement factor (array of size n.bins.f).}
  \item{DSB.per.Gy.per.domain}{ number of DSBs per domain per Gy.}
  \item{domains.per.nucleus}{ number of domains in nucleus.}
  \item{write.output}{ if true, a log file will be written
 ("dose.to.DSBs.log") containing the DBS distribution.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{total.pDSBs}{ probability sum of DSB probability (quality check, has
 to be ~1)}
  \item{total.nDSBs}{ number of DSBs in nucleus}
  \item{number.of.iDSBs}{ number of isolated DSBs in nucleus}
  \item{number.of.cDSBs}{ number of complex DSBs in nucleus}
  \item{avg.number.of.DSBs.in.cDSBs}{ average number of DSBs in complex DSBs}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Suc
cessiveConvolutions.c#L1130}
}
