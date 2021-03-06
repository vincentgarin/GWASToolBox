% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QTL_select.R
\name{QTL_select}
\alias{QTL_select}
\title{QTL candidates selection}
\usage{
QTL_select(Qprof, threshold = 3, window = 20, verbose = TRUE)
}
\arguments{
\item{Qprof}{Object of class \code{QTLprof} returned by the function
\code{\link{mpp_SIM}}, \code{\link{mpp_CIM}} or \code{\link{MQE_CIM}}.}

\item{threshold}{\code{Numeric} value representing -log10(p-value) threshold
above which a position can be considered as a QTL candidate. Default = 3.}

\item{window}{\code{Numeric} value in centi-Morgan representing the minimum
distance between two selected positions. Default = 20.}

\item{verbose}{\code{Logical} value specifying if the detection of no QTL
should be printed. Default = TRUE.}
}
\value{
Return:

\item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
1) QTL marker or in between position names; 2) chromosomes;
3) interger position indicators on the chromosome;
4) positions in centi-Morgan; and 5) -log10(p-values).}
}
\description{
Selection of QTL candidate positions.
}
\details{
The function select QTL positions that are above the given
\code{threshold} per chromosome. Once a position has been selected, and
exlusion \code{window} is set around that position. Positions falling into
that region will not be candidate anymore. The search continue until there
is no more candidate position above the \code{threshold}.
}
\references{
This function is a modification of the QTL.reduce function
coming from the Biometris pipepline.

RAP (R Analytical Pipeline) (V0.9.1) May 2011

Authors: Paul Eilers (1), Gerrit Gort (1), Sabine Schnabel (1), Lucia
Gutierrez(1, 2), Marcos Malosetti(1), Joost van Heerwaarden, and Fred van
Eeuwijk(1)

(1) Wageningen University and Research Center, Netherlands (2) Facultad de
Agronomia, UDELAR, Uruguay
}
\seealso{
\code{\link{mpp_SIM}}, \code{\link{mpp_CIM}}, \code{\link{mpp_perm}}
}
