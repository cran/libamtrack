% TODO File path/AT.CSDA.range.m.Rd
\name{AT.CSDA.range.m}
\alias{AT.CSDA.range.m}
\title{AT.CSDA.range.m}
\description{Returns CSDA range (in m) from pstar tables for given energy.
 In case of ions a simple scaling procedure (A/Z^2) will be used (even effective charge will be neglected)}
\usage{AT.CSDA.range.m(E.MeV.u, particle.no, material.no)
}
\arguments{
  \item{E.MeV.u}{ energy of particles in the mixed particle field (array of size number.of.particles) (see also \code{\link{E.MeV.u}}).}
  \item{particle.no}{ type of the particles in the mixed particle field (array of size number.of.particles) (see also \code{\link{particle.no}}).}
  \item{material.no}{ material index (see also \code{\link{material.no}}).}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{CSDA.range.m}{ (array of size number.of.particles) to be allocated by the user which will be used to return the results}
}
\seealso{
View the C source code here:

\url{http://sourceforge.net/apps/trac/libamtrack/browser/tags/0.5.4/src/AT_DataLET.c#L127}
}
\examples{
require(lattice)
df                <-    expand.grid( particle.name         = c("1H", "3He", "12C", "16O"),              # Define parameter space:
                                     particle.energy.MeV.u = 10^seq(-1, 3, length.out = 500),           # 1 nuclid, energy between 0.1 and 1000 MeV/u,
                                     CSDA.range.m          = 0)
for (i in 1:nrow(df)){
    df$CSDA.range.m[i]        <-    AT.CSDA.range.m( E.MeV.u               = df$particle.energy.MeV.u[i],
                                                     particle.no           = AT.particle.no.from.particle.name(df$particle.name[i]),
                                                     material.no           = AT.material.no.from.material.name("Water, Liquid"))$CSDA.range.m
}
xyplot( log10(CSDA.range.m) ~ log10(particle.energy.MeV.u)|"CSDA ranges in water",
        df,
        groups      =    particle.name,
        type        =    "l",
        aspect      =    1,
        as.table    =    TRUE,
        auto.key    =    list(    columns = 4, points = FALSE, lines = TRUE),
        xlab        =    list(    "energy / (MeV/u)", cex = 1.25),
        ylab        =    list( "CSDA range", cex = 1.25),
        panel       =    function(...){
                             panel.grid(h = -1, v = -1)
                             panel.xyplot(...)
                         },
        scales      =    list( y = list(    at        =    -6:0,
                                            labels    =    c("1 um", "10 um", "0.1 mm", "1 mm", "1 cm", "10 cm", "1 m"),
                                            cex       =    1.0),
                               x = list(    at        =    log10(c(.1, .3, 1, 3, 10, 30, 100, 300, 1000)),
                                            labels    =    c("0.1", "0.3", "1", "3", "10", "30", "100", "300", "1000"),
                                            cex       =    1.0))
    )
}


