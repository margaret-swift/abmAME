% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_options.R
\name{sample_options}
\alias{sample_options}
\title{Sample from weighted options}
\usage{
sample_options(W)
}
\arguments{
\item{W}{A vector of weights.}
}
\value{
An integer corresponding to the chosen weight in the initially
provided vector. \strong{NOTE: indexing begins at 0 matching C++ convention}.
}
\description{
The R function to call the C++ \emph{cpp_sample_options} function.
}
\details{
The function will pre-normalise provided values before passing them
to the Rcpp sample function and returning the index randomly selected.
}
\examples{
set.seed(2022)

sampleOut <- NULL
for(i in 1:10000){
  sampleOut[i] <- sample_options(c(0.25, 0.15, 0.5, 0.05, 0.05))
}
hist(sampleOut)
table(sampleOut) / 10000

}
