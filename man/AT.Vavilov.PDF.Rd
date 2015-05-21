% TODO File path/AT.Vavilov.PDF.Rd
\name{AT.Vavilov.PDF}
\alias{AT.Vavilov.PDF}
\title{AT.Vavilov.PDF}
\description{Compute Vavilov distribution using CERNLIB (G116)
}
\usage{AT.Vavilov.PDF(lambda.V, kappa, beta)
}
\arguments{
  \item{lambda.V}{ lambda (array of size n).}
  \item{kappa}{ straggling parameter.}
  \item{beta}{ relativistic speed, between 0 and 1.}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{density}{ resulting density (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_Ene
rgyLoss.c#L32}
}
