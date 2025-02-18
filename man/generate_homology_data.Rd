% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spagi2_master.R
\name{generate_homology_data}
\alias{generate_homology_data}
\title{generate_homology_data}
\usage{
generate_homology_data(species1 = "hsapiens", species2)
}
\arguments{
\item{species1}{Species 1 name, in the form 'hsapiens'. Default is 'hsapiens'.}

\item{species2}{Species 2 name, in the form 'hsapiens'}
}
\value{
This function returns a list containing the homology table and matrix of the two species.
}
\description{
This function generates homology table and matrix for all ensembl genes of two species.
}
\details{
This function generates homology table and matrix for all ensembl genes of two species.
}
\examples{
human.zebrafish.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "drerio")
summary(human.zebrafish.homology.data)
}
