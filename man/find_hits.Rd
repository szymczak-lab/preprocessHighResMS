% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_features.R
\name{find_hits}
\alias{find_hits}
\title{Find putative theoretical features for each measured feature}
\usage{
find_hits(info.features, dat, ppm)
}
\arguments{
\item{info.features}{data.frame with information about theoretical possible
features (e.g. as generated by the function
\code{\link{chem_formula_2_adducts}}).}

\item{dat}{data.frame with information about measured features with columns
id (= unique identifier), mz (= measured mz value), intensity (= measured
intensity).}

\item{ppm}{Numeric. Resolution of mass spectrometer.}
}
\value{
matrix with measured features in rows and theoretical possible features in
columns (1 = measured feature within ppm distance to theoretical feature,
0 = otherwise).
}
\description{
Identifies all possible theoretical features within an user defined ppm
window around the measured mz value of each measured feature. The result
includes only those theoretical features that are potential candidates for
at least one measured feature and orphan isotopes are removed if no hit is
identified for the corresponding monoisotope.
}
\examples{
data("se.example")
data("info.features")

dat = prepare_data_for_annotation(se = se.example)
hits.m = find_hits(info.features = info.features,
                   dat = dat,
                   ppm = 20)
}
