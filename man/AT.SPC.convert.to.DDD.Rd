\name{AT.SPC.convert.to.DDD}
\alias{AT.SPC.convert.to.DDD}
\title{AT.SPC.convert.to.DDD}
\description{Converts spectral fluence data differential in depth and energy
 (spc) into a depth-dose. It also allows for the export of TRiP98 compatible
 DDD files. The routine calls AT.SPC.read so if you experience problems,
 please check if additional arguments for AT.SPC.read are needed (via '...').
 For the DDD format, see M. Kr\"amer's TRiP98 documentation:
 \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtddd.html}.}
\usage{AT.SPC.convert.to.DDD(file.name.spc, file.name.ddd = NULL, endian =
 'little', plot = TRUE, write = TRUE, ...)}
\arguments{
  \item{file.name.spc}{ path and file name for spc data}
  \item{file.name.ddd}{ path and file name for ddd data, if not given
 file.name.spc will be used but with *.ddd extension (only applicable if
 \code{write = TRUE}}
  \item{endian}{ endianess of spc data, see also \link{AT.SPC.read}}
  \item{plot}{ if true, the resulting depth-dose curve will be plotted, e.g.
 for error-checking}
  \item{write}{ if true, the resulting depth-dose curve will be exported in
 TRiP98 format}
  \item{...}{ additional arguments for the \link{AT.SPC.read} routine}
}
\value{
  A data frame containing the depth-dose curve
}
