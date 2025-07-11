% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataResampling.R
\name{resamp.2x}
\alias{resamp.2x}
\title{randomizing a community phylogenetics matrix: '2x' of Hardy (2008)}
\usage{
resamp.2x(obj, level = 0.1)
}
\arguments{
\item{obj}{a community dataset in spacodiR format (see as.spacodi)}

\item{level}{numeric vector of length one giving a proportion specifying the extent of data shuffling}
}
\value{
A shuffled dataset
}
\description{
Used for resampling data within a community dataset
}
\details{
A resampling procedure for a species-by-plots matrix, based on Gotelli swapping. 
Shuffles abundances within a pair of plots and for a pair of species. 
The level defines the degree of sampling, with larger values dictating a higher level of reshuffling. 
For instance, if level = 0.4 and the dataset involves 5 species and 10 plots, a total of 20 (0.4x5x10) Gotelli swaps are performed
}
\examples{
data(sp.example)
attach(sp.example)
spl


# shuffle dataset
resamp.2x(obj=spl, level=0.2) 

}
\references{
HARDY OJ. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
GOTELLI NJ. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621.
}
\author{
Timothy Paine and Jonathan Eastman
}
