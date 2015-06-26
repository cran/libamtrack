% TODO File path/AT.Vavilov.IDF.Rd
\name{AT.Vavilov.IDF}
\alias{AT.Vavilov.IDF}
\title{AT.Vavilov.IDF}
\description{Computes the Vavilov probability density function using CERNLIB
 (G116)
}
\usage{AT.Vavilov.IDF(rnd, kappa, beta)
}
\arguments{
  \item{rnd}{ random number from uniform distribution between 0 and 1 (array
 of size n).}
  \item{kappa}{ straggling parameter (array of size n).}
  \item{beta}{ relativistic speed, between 0 and 1 (array of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{lambda.vavilov}{ resulting Vavilov lambda (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.2/src/AT_Ene
rgyLoss.c#L279}
}
