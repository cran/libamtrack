% TODO File path/AT.D.RDD.Gy.Rd
\name{AT.D.RDD.Gy}
\alias{AT.D.RDD.Gy}
\title{AT.D.RDD.Gy}
\description{Returns local dose as a function of distance r_m for a given
 radial dose distribution model
}
\usage{AT.D.RDD.Gy(r.m, E.MeV.u, particle.no, material.no, rdd.model, 
	rdd.parameter, er.model, stopping.power.source.no)
}
\arguments{
  \item{r.m}{ distance [m] (array of size n).}
  \item{E.MeV.u}{ particle (ion) energy per nucleon [MeV/u] (single number, no
 mixed fields) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ particle code number (single number, no mixed fields)
 (see also \code{\link{particle.no}}).}
  \item{material.no}{ material code number (single number, no mixed fields)
 (see also \code{\link{material.no}}).}
  \item{rdd.model}{ radial dose distribution model index (see also
 \code{\link{rdd.model}}).}
  \item{rdd.parameter}{ radial dose distribution model parameters (array of
 size 4).}
  \item{er.model}{ electron range / track with model index (see also
 \code{\link{er.model}}).}
  \item{stopping.power.source.no}{ TODO (see also
 \code{\link{stopping.power.source.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD
% NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{D.RDD.Gy}{ dose [Gy] (array of size n)}
  \item{status}{status}
}
\seealso{
View the C source code here:
\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.6.1/src/AT_RDD
.c#L485}
}
\examples{

# Compute dose in several distances of an 100 MeV/u neon ion in water
# according to 'Site' parametrization
AT.D.RDD.Gy(    r.m              = 10^(-9:-4),
                E.MeV.u          = 100,
                particle.no      = 10020,
                material.no      = 1,
                rdd.model        = 4,
                rdd.parameter    = c(5e-8, 1e-10),
                er.model         = 2,
                stopping.power.source.no = 2)

# Compare the Geiss parametrization of RDD for protons and Carbon ions at
# different energies:
df    <-    expand.grid( E.MeV.u        = 10^seq(0, 3, length.out = 4),      
 # from 1 to 1000 MeV/u in 4 steps
                         particle.no    = c(1001,6012),                      
 # protons and carbons
                         r.m            = 10^seq(-9, -2, length.out = 100),  
 # from 1 nm to 1 cm in 100 steps
                         material.no    = 2,                                 
 # Aluminium Oxide
                         rdd.model      = 3,                                 
 # Geiss parametrization
                         rdd.parameter  = 5e-8,                              
 # Fixed core size of 50 nm
                         er.model       = 4,                                 
 # Geiss track width parametrization
                         D.Gy           = 0)                                 
 # For later use
ii                   <-  df$particle.no == 1001                              
 # Add particle names
df$particle.name     <-  "Carbon-12"
df$particle.name[ii] <-  "Protons"
for (i in 1:nrow(df)){                                                       
 # Loop through particles/energies
    df$D.Gy[i]    <-    AT.D.RDD.Gy( r.m              = df$r.m[i],
                                     E.MeV.u          = df$E.MeV.u[i],
                                     particle.no      = df$particle.no[i],
                                     material.no      = df$material.no[i],
                                     rdd.model        = df$rdd.model[i],
                                     rdd.parameter    = df$rdd.parameter[i],
                                     er.model         = df$er.model[i],
                                     stopping.power.source.no = 2)$D.RDD.Gy
}
}
