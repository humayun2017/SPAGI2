% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spagi2_master.R
\name{identify_active_pathway_path_new}
\alias{identify_active_pathway_path_new}
\title{identify_active_pathway_path_new}
\usage{
identify_active_pathway_path_new(pathway.path, processed.query.data)
}
\arguments{
\item{pathway.path}{A list with pathway path data where each sublist denotes a path of the pathway. This is used as background data.}

\item{processed.query.data}{A list with expressed query data where each sublist corresponds for each cell/tissue type.}
}
\value{
This function returns a list of pathways where each sublist denote a path of the pathway.
}
\description{
This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
}
\details{
This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
}
\examples{
#Pre-process the 'tooth.epi.E13.5' data
tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#Generate the homology data for the two species - 'hsapiens' and 'mmusculus'
mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
#Generate the mouse homology pathway path data
mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
#Identify active pathway paths of the processed query data
tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
head(summary(tooth.epi.E13.5.active.pathway$E13.5_Epi))

}
