% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.R
\name{bridgeNormalise}
\alias{bridgeNormalise}
\title{Adjust quantification values in MSnSet}
\usage{
bridgeNormalise(
  obj,
  plex1_ix = 1:10,
  plex2_ix = 12:21,
  bridge1_ix = 11,
  bridge2_ix = 22,
  keep_bridge = 0
)
}
\arguments{
\item{plex1_ix}{index for quantification values in multiplex 1}

\item{plex2_ix}{index for quantification values in multiplex 2}

\item{bridge1_ix}{index for bridge sample(s) in multiplex 1}

\item{bridge2_ix}{index for bridge sample(s) in multiplex 1}

\item{keep_bridge}{How many values to retain for each bridge sample (0, 1, or 2)}

\item{values}{vector of quantification values}
}
\value{
bridge normalised values
}
\description{
Adjust quantification values in MSnSet using bridging channels to
normalise across isobaric tag multiplexes.
Note that the bridge channels are not normalised.
}
