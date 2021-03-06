% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute.R
\name{SnKnn}
\alias{SnKnn}
\title{Impute missing values using KNN imputation with prior sum normalisation}
\usage{
SnKnn(obj, k = 5, verbose = FALSE)
}
\arguments{
\item{obj}{Input}

\item{k}{nearest neighbours}

\item{verbose}{Log the start and end time for imputation}
}
\value{
Input data with missing values imputed
}
\description{
Impute missing values using KNN imputation with prior sum normalisation
}
