% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spagi2_master.R
\docType{package}
\name{spagi2-package}
\alias{spagi2-package}
\alias{spagi2}
\alias{SPAGI2}
\title{SPAGI2: Identification of active signalling pathways for cross-species by integrating gene expression and human protein interaction network}
\description{
SPAGI2 is an R package for cross species active signalling pathway identification for a gene expression profile. This package contains the neccessary R code to run SPAGI2 as described in "SPAGI2: Identification of active signalling pathways for cross-species by integrating gene expression and human protein interaction network". SPAGI2 is implemented within a helpful framework to identify active pathways from gene expression profiles.
}
\examples{
## Do a sample analysis using mouse embryonic dental epitheliam cell at E13.5 microarray gene expression data.

## Here we will use "pathway.path.new" and "rp.median.exp.4" as background data from the SPAGI2 repository.
## We will also use "housekeeping.gene" as background data from the SPAGI repository.
## Also we will use "tooth.epi.E13.5" as query microarray gene expression data. This data is for mouse embryonic dental epitheliam cell at E13.5.
## These data sets are loaded automatically with the packages SPAGI2 and SPAGI.

## Pre-process the query data (tooth.epi.E13.5), The data has already been made in normalized and log2 scale format. Here, we will use the expression cutoff as 5.0 of the query data. Also we have already made the replicate names same for the data.
tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
## Generate the homology data for the two species - 'hsapiens' and 'mmusculus'. It will take little bit more time to access the biomaRt ensembl server.
mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
## Generate the mouse homology pathway path data. Here we will use the background human pathway.path.new data to get the mouse homology pathway path data.
mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
## Identify active pathway paths of the processed query data. Here, we will use the mouse.homology.pathway.path data to get the active pathway paths of the processed query data.
tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
## Get pathway activity score (i.e., pathway name gene expression and pathway specific gene count proportion) of the processed query data. This function uses automatically loaded housekeeping gene data from SPAGI package to get the cell/tissue expressed specific pathway gene count proportion.
tooth.epi.E13.5.active.pathway.score<-get_pathway_activity_score_new(active.pathway.path = tooth.epi.E13.5.active.pathway, processed.query.data = tooth.epi.E13.5.processed.data, homology.data = mouse.homology.data)
## Plot the pathway activity score result (i.e., pathway name gene expression versus pathway specific gene count proportion) in a 2D plane (black=specifically active, gray=generally active). This function uses generally highly expressed receptor data from the package to separate the cell/tissue specific active pathways and generally active pathways in many cells.
display_pathway_activity_score(pathway.activity.score = tooth.epi.E13.5.active.pathway.score, homology.data = mouse.homology.data)
## To separate the top ranked pathways we can do this
abline(v=0.2, h=10, lty=2, col="black")

}
\author{
Md Humayun Kabir <humayun.mkabir@gmail.com>
}
