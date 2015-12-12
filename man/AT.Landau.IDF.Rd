% TODO File path/AT.Landau.IDF.Rd
\name{AT.Landau.IDF}
\alias{AT.Landau.IDF}
\title{AT.Landau.IDF}
\description{Computes the Landau inverse distribution function using CERNLIB
 (G115)
}
\usage{AT.Landau.IDF(rnd)
}
\arguments{
  \item{rnd}{ random number from uniform distribution between 0 and 1 (array
 of size n).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{lambda.landau}{ resulting Landau lambda (array of size n)}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.3/src/AT_Ene
rgyLoss.c#L85}
}
