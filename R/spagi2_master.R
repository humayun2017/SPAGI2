#########################################################################################
#' @title SPAGI2: Identification of active signalling pathways for cross-species by integrating gene expression and human protein interaction network
#'
#' @description SPAGI2 is an R package for cross species active signalling pathway identification for a gene expression profile. This package contains the neccessary R code to run SPAGI2 as described in "SPAGI2: Identification of active signalling pathways for cross-species by integrating gene expression and human protein interaction network". SPAGI2 is implemented within a helpful framework to identify active pathways from gene expression profiles.
#'
#'
#' @author Md Humayun Kabir <humayun.mkabir@gmail.com>
#'
#' @docType package
#' @name spagi2-package
#' @aliases spagi2 SPAGI2
#'
#' @examples
#' ## Do a sample analysis using mouse embryonic dental epitheliam cell at E13.5 microarray gene expression data.
#'
#' ## Here we will use "pathway.path.new" and "rp.median.exp.4" as background data from the SPAGI2 repository.
#' ## We will also use "housekeeping.gene" as background data from the SPAGI repository.
#' ## Also we will use "tooth.epi.E13.5" as query microarray gene expression data. This data is for mouse embryonic dental epitheliam cell at E13.5.
#' ## These data sets are loaded automatically with the packages SPAGI2 and SPAGI.
#'
#' ## Pre-process the query data (tooth.epi.E13.5), The data has already been made in normalized and log2 scale format. Here, we will use the expression cutoff as 5.0 of the query data. Also we have already made the replicate names same for the data.
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' ## Generate the homology data for the two species - 'hsapiens' and 'mmusculus'. It will take little bit more time to access the biomaRt ensembl server.
#' mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
#' ## Generate the mouse homology pathway path data. Here we will use the background human pathway.path.new data to get the mouse homology pathway path data.
#' mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
#' ## Identify active pathway paths of the processed query data. Here, we will use the mouse.homology.pathway.path data to get the active pathway paths of the processed query data.
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' ## Get pathway activity score (i.e., pathway name gene expression and pathway specific gene count proportion) of the processed query data. This function uses automatically loaded housekeeping gene data from SPAGI package to get the cell/tissue expressed specific pathway gene count proportion.
#' tooth.epi.E13.5.active.pathway.score<-get_pathway_activity_score_new(active.pathway.path = tooth.epi.E13.5.active.pathway, processed.query.data = tooth.epi.E13.5.processed.data, homology.data = mouse.homology.data)
#' ## Plot the pathway activity score result (i.e., pathway name gene expression versus pathway specific gene count proportion) in a 2D plane (black=specifically active, gray=generally active). This function uses generally highly expressed receptor data from the package to separate the cell/tissue specific active pathways and generally active pathways in many cells.
#' display_pathway_activity_score(pathway.activity.score = tooth.epi.E13.5.active.pathway.score, homology.data = mouse.homology.data)
#' ## To separate the top ranked pathways we can do this
#' abline(v=0.2, h=10, lty=2, col="black")
#'
NULL
####################






##################################################################################
# Copyright University of Rajshahi 2020
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################






######################################################################################
# Users should also check the names of functions loaded by this script
#
#
######################################
# Please note that it is necessary to pre-process the RNA-seq query data to calculate RPKM/FPKM/CPM of raw read count data or microarray query data to normalized form and then make the data into log2 scale format before utilizing these data with the SPAGI2 package.
# Also note that the background pathway path and housekeeping gene data are in official gene symbol format. So please make your gene ids as valid gene ids supported by ensembl before using with SPAGI2 package.
# The SPAGI2 or SPAGI package does not perform any normalizations or log2 scale format. It assumes that all RNA-seq query data are in RPKM/FPKM/CPM normalized and log2 scale format. For micro-array data sets be aware that extra steps like background subtraction and quantile normalization are recommended before log2 scale format.
# To utilize the SPAGI2 package, you have to provide an expression cut-off threshold and a high expression threshold (i.e., an expression value that is high enough) for your query data. To look at the distribution of the query data you can unilize the helping function 'show_distribution()' of the package and then fix the expression cut-off threshold and a high expression threshold values.
######################################
#
#
######################################
# It is required to provide the valid gene ids and species name supported by ensembl. For invalid gene ids and species name one can not get the result.
# To see the valid gene Ids of a species that have mapping to the external_gene_name, one can look at the find_valid_geneID() function of the package.
# Also, to see the currently ensembl supported species , one can look at the find_supported_datasets() function of the package.
######################################
#
#
#
##############################
# Sometimes you might get temporarily disconnected from the BioMart web server with an error message like:
#
#  Error in bmRequest(request = request, verbose = verbose) :
#    Internal Server Error (HTTP 500).
#
# If this happens just try to re assign the variable 'host.ID' with an archived version.
# For example -
# host.ID <- "http://sep2019.archive.ensembl.org"
#
# Alternatively the BioMart web service is temporarily down. Just wait for the BioMart web service up again.
###############################
#
#
#
###########################################
# SPAGI2 depends on seven libraries - biomaRt, slam, fastmatch, data.table, doParallel, igraph and spagi
# Also to generate the background pathway data it depends on STRINGdb. Make sure you have installed all the packages.
# To install these libraries execute the following commands:
#
# If 'BiocManager' is not already installed, first install 'BiocManager' to insall the bioconductor packages. Type the following in an R command window:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
# BiocManager::install("STRINGdb")
# install.packages("slam")
# install.packages("fastmatch")
# install.packages("data.table")
# install.packages("doParallel")
# install.packages("igraph")
#
# Make sure you have *devtools* installed, then install the *SPAGI* and *SPAGI2* package from github repository:
# install.packages('devtools')
# devtools::install_github('VCCRI/SPAGI')
# devtools::install_github('humayun2017/SPAGI2')
#
# Finally load all the libraries with
# library(biomaRt)
# library(STRINGdb)
# library(slam)
# library(fastmatch)
# library(data.table)
# library(doParallel)
# library(igraph)
# library(spagi)
# library(spagi2)
#############################################
#
#
#
####################
##!!When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed,
##!!When used in package, these variables need to be global..
##!!So, assign these variables before running the package.
#biomart.ID <- "ENSEMBL_MART_ENSEMBL"
#host.ID <-  "www.ensembl.org"
##or,
#host.ID <- "asia.ensembl.org"
##or, for example to an archived version
#host.ID <- "http://sep2019.archive.ensembl.org"
####################
#######################################################################################






#########################################################################################
# For valid gene Ids of a species that have mapping to the external_gene_name
# This function returns the valid gene ids of the species that have mapping to the ensembl_gene_id

#' @title find_valid_geneID
#'
#' @description
#' This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species. For supported species that have ensembl id mapping run 'supportedSpecies' (this command is from XGSA package), it will return all the names of the supported species. For each species there exists separate ensembl datasets.
#'
#' @rdname find_valid_geneID
#' @name find_valid_geneID
#'
#' @details This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species. For supported species that have ensembl id mapping run 'supportedSpecies' (this command is from XGSA package), it will return all the names of the supported species. For each species there exists separate ensembl dataset.
#' @return This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species
#'
#' @param species Species name in the form "hsapiens"
#'
#' @importFrom biomaRt useMart listAttributes
#'
#' @export
#'
#' @examples
#' human.valid.gene.IDs<-find_valid_geneID(species = "hsapiens")
#' head(human.valid.gene.IDs)

find_valid_geneID <- function(species = "hsapiens"){
  dataset.name <- paste(species, "_gene_ensembl", sep="")
  if(interactive()){
    ensembl = useMart("ensembl", dataset=dataset.name)
    return(listAttributes(ensembl))
  }
}
####################






#############################################################################################
# This function determines which species are currently supported by Ensembl - taken from xgsa

#' @title find_supported_datasets
#'
#' @description
#' This function determines which species are currently supported by Ensembl
#'
#' @rdname find_supported_datasets
#' @name find_supported_datasets
#'
#' @details
#' This function determines which species are currently supported by Ensembl
#'
#' @return This function returns a vector or organism names in the form 'hsapiens'
#'
#'
#' @importFrom biomaRt useMart listDatasets
#'
#' @export
#'
#' @examples
#' supported_species <- find_supported_datasets()
#' head(supported_species)

find_supported_datasets <- function(default=TRUE){
  ensembl <- useMart(biomart.ID, dataset='hsapiens_gene_ensembl', host=host.ID)
  sets <- listDatasets(ensembl)
  organisms <- gsub("([a-z]*)_.*","\\1",sets[,1],"_")
  return(organisms)
}
####################






###########################################################################################################
## look for the distribution for all cells / tissues based on their replicates' average value
## this function uses the 'tooth.epi.E13.5' data to test the function that is automatically loaded with the package

#' @title show_distribution
#'
#' @description
#' This function helps to look for the distribution for all cells / tissues based on their replicates' average value.
#'
#' @rdname show_distribution
#' @name show_distribution
#'
#' @details
#' This function helps to look for the distribution for all cells / tissues based on their replicates' average value.
#'
#' @return This function displays the distribution of each cell / tissue based on their replicates' average value. It returns null value.
#'
#' @param cell.tissue.data An expression matrix with replicated column headers per-replicate. It is assumed that all query data are in RPKM/FPKM/CPM and log normalized form. Also assume that gene ids are official gene symbols. For the matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param data.format Format of cell.tissue.data. Default is "matrix".
#' @param experiment.descriptor A vector corresponding to the matrix column names of cell.tissue.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#'
#' @importFrom spagi format_matrix_data compute_mean
#'
#' @export
#'
#' @examples
#' show_distribution(cell.tissue.data = tooth.epi.E13.5)
#' abline(v=5,col="red")

show_distribution<-function(cell.tissue.data, data.format = "matrix", experiment.descriptor = NULL){
  if(data.format=="matrix"){
    #calculate average value for each gene based on replicates for each cell type /tissue
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(cell.tissue.data)
    cell.data.formatted<-format_matrix_data(cell.tissue.data, experiment.descriptor)

    #show the distribution for each cell/tissue based on average value of the replicates
    cell.tissue.names<-colnames(cell.data.formatted)
    for(i in 1:ncol(cell.data.formatted)){
      title.text<-paste("Distribution of", cell.tissue.names[i], sep = " ")
      devAskNewPage(ask = FALSE) #to turn off the "Hit <Return> to see next plot" prompt
      hist(cell.data.formatted[,i], breaks=100, main=title.text, xlim=c(0,15), ylim=c(0,2000),
           xlab = "Expression value", ylab = "# of genes",  col="grey")
      #Sys.sleep(1)
    }
    return(NULL)
  }else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }
}
####################






###############################################################################################
# Returns a data frame with two columns containing all gene symbols and their corresponding gene IDs (eg. Ensembl gene IDs) for the given species
# Input parameter is the species string like 'hsapiens'
# gene mapping to external_gene_name with other gene id
# Here, geneID.forMapping is used for the existing gene ID (that has 'external_gene_name' mapping in the biomaRt database), such as "ensembl_gene_id" for ensembl gene id.
# For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the package.

#' @title make_ENSEMBL_symbol_map_new
#'
#' @description
#' This helper function returns a data frame with two columns containing all gene symbols and their corresponding gene IDs for the given species
#'
#' @rdname make_ENSEMBL_symbol_map_new
#' @name make_ENSEMBL_symbol_map_new
#'
#' @details This helper function returns a data frame with two columns containing all gene symbols and their corresponding gene IDs for the given species
#' @return This helper function returns a data frame with two columns containing all gene symbols and their corresponding gene IDs for the given species
#'
#' @param species Species name in the form 'hsapiens'
#' @param geneSymbol.toMapping The valid gene ID names to which the user cell data will be converted, default is "external_gene_name"
#' @param geneID.forMapping The valid gene ID names of the user cell data, default is "ensembl_gene_id"
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human.ensembl.symbol.map<-make_ENSEMBL_symbol_map_new('hsapiens', geneSymbol.toMapping="external_gene_name", geneID.forMapping="ensembl_gene_id")
#' head(human.ensembl.symbol.map)
#'

make_ENSEMBL_symbol_map_new <- function(species, geneSymbol.toMapping, geneID.forMapping){
  dataset.name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart = biomart.ID, dataset = dataset.name, host = host.ID)
  ENSEMBL.symbol.map <- getBM(attributes=c(geneSymbol.toMapping, geneID.forMapping), mart = ensembl)
  return(ENSEMBL.symbol.map)
}
####################






############################################################################################################
# This function converts gene IDs to gene symbols for cells/tissues gene expression profiles
# for example, used external_gene_name for geneSymbol.toMapping and ensembl_gene_id for geneID.forMapping
# the example gene symbols should be - rownames(query.data)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")

#' @title convert_geneID_to_geneSymbol
#'
#' @description
#' This function converts either a vector of gene IDs, a matrix of gene or a named cell / tissue expression profile or matrix to ensembl IDs
#'
#' @rdname convert_geneID_to_geneSymbol
#' @name convert_geneID_to_geneSymbol
#'
#' @details
#' This function converts gene IDs of a named cell/tissue expression profile to gene symbols
#'
#' @return This function returns cell/tissue profiles with converted gene symbols
#'
#' @param user.cell.data Cell/tissue profile with valid gene IDs
#' @param species The species of the given gene IDs or cell/tissue profile, default is "hsapiens"
#' @param geneSymbol.toMapping The valid gene ID names (that has ensembl id mapping in the biomaRt database) to which the user.cell.data will be converted, default is "external_gene_name". For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the package.
#' @param geneID.forMapping The valid gene ID names (that has ensembl id mapping in the biomaRt database) of the user.cell.data, default is "ensembl_gene_id". For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the package.
#' @param collapse.method Used when one external_gene_name has more then one probes, usually two options are used, maximum probe value (i.e., "max") or average of probe values (i.e., "mean"), default is "max".
#'
#' @export
#'
#' @examples
#' query.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(query.data)<-c("ENSG00000160202", "ENSG00000109846", "ENSG00000100053", "ENSG00000007372", "ENSG00000181449", "ENSG00000117707", "ENSG00000138083", "ENSG00000150938", "ENSG00000244752", "ENSG00000101144")
#' colnames(query.data)<-c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' convert_geneID_to_geneSymbol(query.data)
#'

convert_geneID_to_geneSymbol <- function(user.cell.data, species="hsapiens", geneSymbol.toMapping="external_gene_name", geneID.forMapping = "ensembl_gene_id", collapse.method = "max"){
  ensembl.symbol.map<-make_ENSEMBL_symbol_map_new(species, geneSymbol.toMapping, geneID.forMapping)

  #for cell/tissue gene expression data sets
  cell.tissue.names<-colnames(user.cell.data)
  list.gene.names<-rownames(user.cell.data)
  ensembl.genes.exist.in.user.cell.data <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% list.gene.names,]

  #duplicate probes / symbol IDS - multiple ENS genes match to the same probe or symbol - probably mostly NA or ""
  duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.cell.data[[2]])
  ensembl.genes.exist.in.user.cell.data.without.duplication <- ensembl.genes.exist.in.user.cell.data[!duplicate.ensembl.gene.list,]

  idlist <- split(ensembl.genes.exist.in.user.cell.data.without.duplication[[2]], ensembl.genes.exist.in.user.cell.data.without.duplication[[1]])

  if(collapse.method == "mean"){
    #for mean
    collapsed <- lapply(idlist,function(X){
      return(colMeans(user.cell.data[unique(X),,drop=F]))
    })
  }
  else if(collapse.method == "max"){
    #for max
    collapsed <- lapply(idlist,function(X){
      return(user.cell.data[unique(X),,drop=F][which.max(rowMeans(user.cell.data[unique(X),,drop=F])),])
    })
  }
  else {
    print("ERROR: Other collapse methods not supported!!")
    return(NULL)
  }

  user.cell.data.clean <- do.call(rbind, collapsed)
  colnames(user.cell.data.clean)<-cell.tissue.names

  return(user.cell.data.clean)
}
####################






###########################################################################################################
# For gene ID conversion to gene symbol of the query data and also -
# For average value calculation of the replicates and then get expressed genes based on expression cut-off;
# return a list of expressed genes for each cell type or tissue

#' @title preprocess_querydata_new
#'
#' @description
#' This function preprocesses the query data to convert gene IDs to gene symbols and to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @rdname preprocess_querydata_new
#' @name preprocess_querydata_new
#'
#' @details
#' This function preprocesses the query data to convert gene IDs to gene symbols and to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @return This function returns a list with specifically expressed genes for each cell type / tissue with gene IDs as gene symbols
#'
#' @param cell.tissue.data An expression matrix with replicated column headers per-replicate. It is assumed that all query data are in RPKM/FPKM/CPM and log normalized form. Also assume that gene ids are official gene symbols. For the matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param exp.cutoff.th An expression cut-off threshold value for the query data.
#' @param species The species abbreviation of the query data (cell.tissue.data). Default is "hsapiens".
#' @param data.format Format of cell.tissue.data. Default is "matrix".
#' @param geneID The code for the type of gene IDs used by cell.tissue.data, as used by the biomaRt database. To find the valid codes for gene IDs for a species, please see the find_valid_geneID() function of the package. Default is "external_gene_name".
#' @param experiment.descriptor A vector corresponding to the matrix column names of cell.tissue.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#' @param collapse.method How to summarise values when one ensembl_gene_id has more then one value (multiple microarray probes or transcripts to one gene for example). Currently two options are implemented, 'max' or 'mean'. Default is "max".
#'
#' @importFrom spagi format_matrix_data compute_mean
#'
#' @export
#'
#' @examples
#' query.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(query.data)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(query.data)<-c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' preprocess_querydata_new(cell.tissue.data=query.data, exp.cutoff.th=5)
#'

preprocess_querydata_new<-function(cell.tissue.data, exp.cutoff.th, species = "hsapiens", data.format = "matrix", geneID = "external_gene_name", experiment.descriptor = NULL, collapse.method = "max"){
  if(data.format=="matrix"){
    ##first convert the gene ids to gene symbols
    if(geneID=="external_gene_name"){
      cell.data.genesymbol<-cell.tissue.data
    }
    else{
      cell.data.genesymbol<-convert_geneID_to_geneSymbol(user.cell.data = cell.tissue.data, species = species, geneID.forMapping = geneID, collapse.method = collapse.method)
    }
    ##

    ##process the query data
    #calculate average value for each gene based on replicates for each cell type
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(cell.data.genesymbol)
    cell.data.formatted<-format_matrix_data(cell.data.genesymbol, experiment.descriptor)

    #get expressed genes for each cell based on cut-off threshold
    expressed.cell.data<-list()
    cell.names<-colnames(cell.data.formatted)
    for(i in 1:ncol(cell.data.formatted)){
      each.cell.data<-cell.data.formatted[,i]
      names(each.cell.data)<-rownames(cell.data.formatted)
      #now get expressed genes based on cut-off value
      each.cell.exp.data<-each.cell.data[each.cell.data>=exp.cutoff.th]
      expressed.cell.data[[cell.names[i]]]<-each.cell.exp.data
    }

    #return the list with expressed genes data for each cell
    return(expressed.cell.data)
    ##
  }
  else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }
}
####################






#######################################################################################################
# This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional)
# It uses the simple_triplet_matrix structure from the package 'slam'
# This function is taken from the XGSA package

#' @title generate_sparse_matrix
#'
#' @description
#' This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
#' It uses the simple_triplet_matrix structure from the package 'slam'.
#'
#' @rdname generate_sparse_matrix
#' @name generate_sparse_matrix
#'
#' @details
#' This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
#' It uses the simple_triplet_matrix structure from the package 'slam'.
#' @return This function return a sparse matrix
#'
#' @param homology_table Table of homology mapping between two species, as retrieved by the function get_homology_matrix()
#'
#' @importFrom slam simple_triplet_matrix
#'
#' @export
#'
#' @examples
#' Not used trivially
#'

generate_sparse_matrix <- function(homology_table){
  #require(slam)
  # Compare only genes with homology mapping
  homology_table <- unique(homology_table)[!unique(homology_table)[,2]=="",]

  if(ncol(homology_table)==3){
    #hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), x = homology_table[,3])
    hms <- simple_triplet_matrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), v = homology_table[,3])
  } else if(ncol(homology_table)==2){
    hms <- simple_triplet_matrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), v = rep(1,nrow(homology_table)))
    #hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])))
  } else {
    print("ERROR: incorrect number of columns in homology table ")
    return(NULL)
  }
  rownames(hms) = unique(homology_table[,1])
  colnames(hms) = unique(homology_table[,2])
  return(hms)
}
####################






#######################################################################################################
# use ENSEMBL to generate homology table
# input parameteres are strings like 'hsapiens' and 'drerio' and a number indicating which species to get sequence identity wrt.
# This function is taken from the XGSA package and is used with some modifications

#' @title get_homology_table_new
#'
#' @description
#' This function uses information from Ensembl to generate a homology table between two species.
#'
#' @rdname get_homology_table_new
#' @name get_homology_table_new
#'
#' @details
#' This function uses information from Ensembl to generate a homology table between two species.
#' @return This function return a data frame with 2 or 3 columns, representing the Ensembl gene IDs in species 1 and 2, and the sequence identity score if requested.
#'
#' @param species1 Species 1 name, in the form 'hsapiens'
#' @param species2 Species 2 name, in the form 'hsapiens'
#' @param sequence_identity_reference Flag of which sequence identity value to retrieve. 0 (default) returns no value); 1 returns the sequence of identity of genes in species 2 compared to species 1; 2 returns the sequence of identity of genes in species 1 compared to species 2.
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_zebrafish_homology_table <- get_homology_table_new('hsapiens', 'drerio')
#' head(human_zebrafish_homology_table)
#'

get_homology_table_new <- function(species1, species2, sequence_identity_reference = 0){
  #require(biomaRt)
  dataset_name <- paste(species1, "_gene_ensembl", sep="")
  homolog_attribute <- paste(species2, "_homolog_associated_gene_name", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)

  if(sequence_identity_reference == 1){
    homology_table <- getBM(attributes=c("external_gene_name", homolog_attribute, paste(species2, "_homolog_perc_id", sep="")), mart = ensembl)
  } else if (sequence_identity_reference == 2){
    homology_table <- getBM(attributes=c("external_gene_name", homolog_attribute, paste(species2, "_homolog_perc_id_r1", sep="")), mart = ensembl)
  } else {
    homology_table <- getBM(attributes=c("external_gene_name", homolog_attribute), mart = ensembl)
  }
  return(homology_table)
}
####################






#############################################################################################
# Returns a data frame with one column containing all Ensembl gene IDs for the given species
# Input parameter is the species string (eg. 'hsapiens')
# This function is taken from the XGSA package with some modifications

#' @title get_ENSEMBL_gene_list_new
#'
#' @description
#' This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#'
#' @rdname get_ENSEMBL_gene_list_new
#' @name get_ENSEMBL_gene_list_new
#'
#' @details This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#' @return This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_ensembl_genes <- get_ENSEMBL_gene_list_new('hsapiens')
#' head(human_ensembl_genes)
#'

get_ENSEMBL_gene_list_new <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  gene.list <- getBM(attributes=c("external_gene_name"), mart = ensembl)
  return(gene.list)
}
####################






######################################################################################################
# Creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix
# Returns a sparse matrix
# Input parameters are the two species strings (eg. 'hsapiens', 'drerio') and a boolean of whether to return the sequence identity values or not
# This function is taken from the XGSA package and used with some modifications

#' @title get_homology_matrix_new
#'
#' @description
#' This function creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix.
#'
#' @rdname get_homology_matrix_new
#' @name get_homology_matrix_new
#'
#' @details
#' This function creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix.
#' @return This function returns a sparse matrix where rows represent the genes in species 1 and columns the genes in species 2.
#'
#' @param species1 Species 1 name, in the form 'hsapiens'
#' @param species2 Species 2 name, in the form 'hsapiens'
#' @param homology.table Homology table for two species
#' @param seq.identity Boolean of whether or not  to retrieve sequence identity information, defaults to FALSE
#'
#' @export
#'
#' @examples
#' human.zebrafish.homology.table<-get_homology_table_new('hsapiens','drerio')
#' human.zebrafish.homology.matix <- get_homology_matrix_new('hsapiens','drerio', human.zebrafish.homology.table)
#' print(data.matrix(human.zebrafish.homology.matix[1:5,1:5]))

get_homology_matrix_new <- function(species1, species2, homology.table, seq.identity = FALSE){
  if(seq.identity == FALSE){
    if(!exists("homology.matrix.list", envir = .GlobalEnv)){
      assign("homology.matrix.list", list(), envir = .GlobalEnv)
    }
    direct <- !is.null(homology.matrix.list[[species1]][[species2]])
    indirect <- !is.null(homology.matrix.list[[species2]][[species1]])
    if(direct){
      # homology matrix already genearted and stored
      return(homology.matrix.list[[species1]][[species2]])
    } else if(indirect){
      # inverse homology matrix already generated and stored
      return(t(homology.matrix.list[[species2]][[species1]]))
    } else {
      if(species1 == species2){
        # do a same species analysis
        gene.list <- get_ENSEMBL_gene_list_new(species1)
        hm <- generate_sparse_matrix(cbind(gene.list, gene.list))
      } else {
        # create a new homology matrix and return it, here 'homology.table' is taken as a function parameter
        hm <- generate_sparse_matrix(homology.table)
      }
      if(species1 %in% names(homology.matrix.list)){
        homology.matrix.list[[species1]][[species2]] <<- hm
      } else {
        homology.matrix.list[[species1]] <<- list()
        homology.matrix.list[[species1]][[species2]] <<- hm
      }
      return(homology.matrix.list[[species1]][[species2]])
    }
  } else {
    print("ERROR: can not generate homology matrix with sequence identity values at this moment")
    return(NULL)
  }
}
####################






#######################################################################################################
##Generate homology table and matrix for all ensembl genes of two species

#' @title generate_homology_data
#'
#' @description
#' This function generates homology table and matrix for all ensembl genes of two species.
#'
#' @rdname generate_homology_data
#' @name generate_homology_data
#'
#' @details
#' This function generates homology table and matrix for all ensembl genes of two species.
#' @return This function returns a list containing the homology table and matrix of the two species.
#'
#' @param species1 Species 1 name, in the form 'hsapiens'. Default is 'hsapiens'.
#' @param species2 Species 2 name, in the form 'hsapiens'
#'
#' @export
#'
#' @examples
#' human.zebrafish.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "drerio")
#' summary(human.zebrafish.homology.data)

generate_homology_data<-function(species1 = "hsapiens", species2){
  #create a list for homology data
  species.homology.data<-list()
  #check for same species and if same species then simply return the 'pathway.path' as 'homology.pathway'
  if(species1==species2){
    species.homology.data[["species.homology.table"]]<-NULL
    species.homology.data[["species.homology.matix"]]<-NULL
    species.homology.data[["species"]]<-species1
  }else{
    #get the homology table for the species
    species.homology.table<-get_homology_table_new(species1, species2)
    #get the homology matrix for the species
    species.homology.matix <- get_homology_matrix_new(species1, species2, species.homology.table)
    species.homology.data[["species.homology.table"]]<-species.homology.table
    species.homology.data[["species.homology.matix"]]<-species.homology.matix
    species.homology.data[["species"]]<-species2
  }
  return(species.homology.data)
}
####################






#######################################################################################################
##Generate the pathway molecule homology table for the two species with respect to the list of molecules of species1

#' @title generate_pathway_molecule_homology_table
#'
#' @description
#' This function creates a pathway molecule homology table for the two species with respect to the list of molecules of species1 by retrieving the homology table.
#'
#' @rdname generate_pathway_molecule_homology_table
#' @name generate_pathway_molecule_homology_table
#'
#' @details
#' This function creates a pathway molecule homology table for the two species with respect to the list of molecules of species1 by retrieving the homology table.
#' @return This function returns a homology table for the list of molecules of species1 of the two species.
#'
#' @param species.homology.data A list containing the homology data (i.e., homology table, homology matrix and species2 name) for the two species.
#' @param species1.pathway.molecule A vector of pathway molecules for species1
#'
#' @export
#'
#' @examples
#' human.zebrafish.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "drerio")
#' human.genes<-c("PAX6", "SOX2", "PROX1", "FGFR1", "BMP7", "MAPK1", "ACE2")
#' human.zebrafish.selected.genes.homology.table<-generate_pathway_molecule_homology_table(species.homology.data = human.zebrafish.homology.data, species1.pathway.molecule = human.genes)
#' human.zebrafish.selected.genes.homology.table

generate_pathway_molecule_homology_table<-function(species.homology.data, species1.pathway.molecule){
  #get the homology table for the species
  species.homology.table<-species.homology.data[[1]]
  #get the homology matrix for the species
  species.homology.matix<-species.homology.data[[2]]
  #get only the species1 pathway molecules that are exist in the homology matrix
  species1.pathway.molecule.filter<-unique(as.character(species1.pathway.molecule[species1.pathway.molecule%in%rownames(species.homology.matix)]))
  #get the pathway homology molecules for the species2 from the homology matrix
  species2.pathway.molecule<-colnames(species.homology.matix)[which(col_sums(species.homology.matix[species1.pathway.molecule.filter,,drop=FALSE])>0)]
  #now filter the homology table with only the hsapiens and other species pathway molecules
  species.homology.table.filtered<-species.homology.table[species.homology.table[[1]] %in% species1.pathway.molecule.filter &
                                                            species.homology.table[[2]] %in% species2.pathway.molecule,]
  rownames(species.homology.table.filtered)<-NULL
  return(species.homology.table.filtered)
}
####################






#######################################################################################
##change existing pathway names by homology pathway molecules with multiple replacement

#' @title alter_pathway_names
#'
#' @description
#' This function alters existing pathway names by homology pathway molecules with multiple replacement.
#'
#' @rdname alter_pathway_names
#' @name alter_pathway_names
#'
#' @details
#' This function alters existing pathway names by homology pathway molecules with multiple replacement.
#' @return This function alters the pathway data with the existing pathways' name change by the homology molecules.
#'
#' @param homology.table Homology table for the pathway molecules of the two species - species1 and species2
#' @param human.pathway A list with human pathway path data where each sublist denotes a path of the pathway. This is used as background data.
#'
#' @export
#'
#' @examples
#' P<-list(c('P','G','A','N'), c('P','C','X','T'), c('P','M','N','D'))
#' Q<-list(c('Q','M','O','T','L'), c('Q','U','T','F'), c('Q','M','N','I','B'))
#' R<-list(c('R','X','B','K'), c('R','G','H','O'), c('R','U','C','E'))
#' QQ<-list(c('Q','M','O','T','L'), c('Q','U','T','F'), c('Q','M','N','I','B'))
#' S<-list(c('S','M','NN'), c('S','N','O','B1'))
#' U<-list(c('SS','OO','nn'), c('ss','oo','hh'))
#' V<-list(c('V','X','B','D','J','E','K'), c('V','M','O','F','D','G'), c('V','B','C'))
#' human.pathway<-list("P"=P, "Q"=Q, "R"=R, "QQ"=QQ, "S"=S, "U"=U, "V"=V)
#' homology.table<-data.frame('human_mol'=c('B1','C','X','R',LETTERS), 'other_mol'=c('b','c1','x1','r1',letters))
#' human.pathway.update<-alter_pathway_names(homology.table, human.pathway)
#' human.pathway.update

alter_pathway_names<-function(homology.table, human.pathway){
  pathway.update<-list()
  pathway.name<-names(human.pathway)
  for(i in 1:length(pathway.name)){
    tmp.elm<-as.vector(homology.table[[2]][which(pathway.name[i]==as.vector(homology.table[[1]]))])
    if(length(tmp.elm)==0){#1st check for non-existing pathways
      NULL
    }else if(length(tmp.elm)==1){#2nd check for non-duplicate pathways - most of them - hence differently
      pathway.update[[tmp.elm]]<-human.pathway[[pathway.name[i]]]
    }else{
      for(j in 1:length(tmp.elm)){#finally for multiple replacement pathways
        pathway.update[[tmp.elm[j]]]<-human.pathway[[pathway.name[i]]]
      }
    }
  }
  return(pathway.update)
}
####################






#############################################################################################
##Alter human pathway elements by homology molecules
#First check whether all elements of a path are exist in the data frame human genes
#then alter each element
#also check for duble genes mapped to single genes, e.g., PAX6=pax6a and PAX6=pax6b for drerio

#' @title alter_pathway_elements
#'
#' @description
#' This function alters human pathway elements by homology molecules of other species.
#'
#' @rdname alter_pathway_elements
#' @name alter_pathway_elements
#'
#' @details
#' This function alters human pathway elements by homology molecules of other species.
#' @return This function returns the pathway data by altering the pathway molecules by the homology molecules.
#'
#' @param homology.table Homology table for the pathway molecules of the two species - species1 and species2
#' @param human.pathway A list with human pathway path data where each sublist denotes a path of the pathway. This is used as background data.
#'
#' @importFrom fastmatch fin
#'
#' @export
#'
#' @examples
#' P<-list(c('P','G','A','N'), c('P','C','X','T'), c('P','M','N','D'))
#' Q<-list(c('Q','M','O','T','L'), c('Q','U','T','F'), c('Q','M','N','I','B'))
#' R<-list(c('R','X','B','K'), c('R','G','H','O'), c('R','U','C','E'))
#' QQ<-list(c('Q','M','O','T','L'), c('Q','U','T','F'), c('Q','M','N','I','B'))
#' S<-list(c('S','M','NN'), c('S','N','O','B1'))
#' U<-list(c('SS','OO','nn'), c('ss','oo','hh'))
#' V<-list(c('V','X','B','D','J','E','K'), c('V','M','O','F','D','G'), c('V','B','C'))
#' human.pathway<-list("P"=P, "Q"=Q, "R"=R, "QQ"=QQ, "S"=S, "U"=U, "V"=V)
#' homology.table<-data.frame('human_mol'=c('B1','C','X','R',LETTERS), 'other_mol'=c('b','c1','x1','r1',letters))
#' homology.pathway<-alter_pathway_elements(homology.table, human.pathway)
#' homology.pathway

alter_pathway_elements<-function(homology.table, human.pathway){
  #first get exiting pathway names change with multiple replacement
  human.pathway.update<-alter_pathway_names(homology.table, human.pathway)

  #next check and replace each pathway path elements
  homology.pathway<-list()
  #homology.table.human.molecule<-sort(unique(as.vector(homology.table[[1]])))
  homology.table.human.molecule<-unique(as.vector(homology.table[[1]]))
  human.pathway.update.name<-names(human.pathway.update)
  for(i in 1:length(human.pathway.update)){
    each.pathway<-human.pathway.update[[i]]
    individual.path.update<-lapply(each.pathway, function(x){
      #for taking only the path where all molecules are present in the data frame human molecules
      #if(!(anyNA(chmatch(x, homology.table.human.molecule)))){
      if(all(x %fin% homology.table.human.molecule)){ #when x has small number of elements it is fast enough
        individual.path.elem.update<-lapply(x, function(z){
          z<-as.vector(homology.table[[2]][which(z==as.vector(homology.table[[1]]))])
          return(z)
        })
      }else{#if any element of a path is not present at human molecule list at homology table
        individual.path.elem.update<-NULL
      }

      #check for non existing path, i.e., if any element of a path is not present at human molecule list
      if(!(is.null(individual.path.elem.update))){
        tmp.ll<-list()
        r<-1
        #then check for number of molecules for each path, in our case all path lengths are 3/4/5/6/7
        if(length(individual.path.elem.update)==3){
          elem1<-human.pathway.update.name[i]
          elem2<-unlist(individual.path.elem.update[[2]])
          elem3<-unlist(individual.path.elem.update[[3]])
          for(j in 1:length(elem2)){
            for(k in 1:length(elem3)){
              tmp.ll[[r]]<-c(elem1,elem2[j],elem3[k])
              r<-r+1
            }
          }
        }else if(length(individual.path.elem.update)==4){
          elem1<-human.pathway.update.name[i]
          elem2<-unlist(individual.path.elem.update[[2]])
          elem3<-unlist(individual.path.elem.update[[3]])
          elem4<-unlist(individual.path.elem.update[[4]])
          for(j in 1:length(elem2)){
            for(k in 1:length(elem3)){
              for(l in 1:length(elem4)){
                tmp.ll[[r]]<-c(elem1,elem2[j],elem3[k],elem4[l])
                r<-r+1
              }
            }
          }
        }else if(length(individual.path.elem.update)==5){
          elem1<-human.pathway.update.name[i]
          elem2<-unlist(individual.path.elem.update[[2]])
          elem3<-unlist(individual.path.elem.update[[3]])
          elem4<-unlist(individual.path.elem.update[[4]])
          elem5<-unlist(individual.path.elem.update[[5]])
          for(j in 1:length(elem2)){
            for(k in 1:length(elem3)){
              for(l in 1:length(elem4)){
                for(m in 1:length(elem5)){
                  tmp.ll[[r]]<-c(elem1,elem2[j],elem3[k],elem4[l],elem5[m])
                  r<-r+1
                }
              }
            }
          }
        }else if(length(individual.path.elem.update)==6){
          elem1<-human.pathway.update.name[i]
          elem2<-unlist(individual.path.elem.update[[2]])
          elem3<-unlist(individual.path.elem.update[[3]])
          elem4<-unlist(individual.path.elem.update[[4]])
          elem5<-unlist(individual.path.elem.update[[5]])
          elem6<-unlist(individual.path.elem.update[[6]])
          for(j in 1:length(elem2)){
            for(k in 1:length(elem3)){
              for(l in 1:length(elem4)){
                for(m in 1:length(elem5)){
                  for(n in 1:length(elem6)){
                    tmp.ll[[r]]<-c(elem1,elem2[j],elem3[k],elem4[l],elem5[m],elem6[n])
                    r<-r+1
                  }
                }
              }
            }
          }
        }else{#for length 7 paths i.e., (length(individual.path.elem.update)==7)
          elem1<-human.pathway.update.name[i]
          elem2<-unlist(individual.path.elem.update[[2]])
          elem3<-unlist(individual.path.elem.update[[3]])
          elem4<-unlist(individual.path.elem.update[[4]])
          elem5<-unlist(individual.path.elem.update[[5]])
          elem6<-unlist(individual.path.elem.update[[6]])
          elem7<-unlist(individual.path.elem.update[[7]])
          for(j in 1:length(elem2)){
            for(k in 1:length(elem3)){
              for(l in 1:length(elem4)){
                for(m in 1:length(elem5)){
                  for(n in 1:length(elem6)){
                    for(o in 1:length(elem7)){
                      tmp.ll[[r]]<-c(elem1,elem2[j],elem3[k],elem4[l],elem5[m],elem6[n],elem7[o])
                      r<-r+1
                    }
                  }
                }
              }
            }
          }
        }
        return(tmp.ll)
      }else{
        return(NULL)
      }
    })

    #check for null paths
    if(!(is.null(individual.path.update))){
      homology.pathway[[human.pathway.update.name[i]]]<-unlist(individual.path.update, recursive = F)
    }
  }
  return(homology.pathway)
}
####################






#############################################################################################
##The main function to replace the molecules of human pathway data by other species molecules

#' @title generate_homology_pathways
#'
#' @description
#' This function generates homology pathways by altering the molecules of human pathway paths.
#'
#' @rdname generate_homology_pathways
#' @name generate_homology_pathways
#'
#' @details
#' This function generates homology pathways by altering the molecules of human pathway paths. The human pathway path data is automatically loaded with the package.
#' @return This function returns the homology pathways by altering the pathway molecules by the homology molecules.
#'
#' @param species.homology.data A list containing the homology data (i.e., homology table, homology matrix and species2 name) for the two species.
#' @param pathway.path A list with human pathway path data where each sublist denotes a path of the pathway. This is used as a background data.
#'
#' @export
#'
#' @examples
#' human.zebrafish.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "drerio")
#' zebrafish.homology.pathway.path<-generate_homology_pathways(species.homology.data = human.zebrafish.homology.data, pathway.path = pathway.path.new)
#' head(zebrafish.homology.pathway.path[[1]])

generate_homology_pathways<-function(species.homology.data, pathway.path){
  #check for same species and if same species then simply return the 'pathway.path' as 'homology.pathway'
  if(species.homology.data$species=="hsapiens"){
    homology.pathway<-pathway.path
    return(homology.pathway)
  }else{
    #get the unique molecules from the human 'pathway.path' background data
    human.pathway.molecules<-unique(as.vector(unlist(pathway.path)))

    #generate the pathway molecule homology table for the two species with respect to the list of molecules of species1
    pathway.molecules.homology.table<-generate_pathway_molecule_homology_table(species.homology.data = species.homology.data, species1.pathway.molecule = human.pathway.molecules)

    #generate human pathway elements by homology molecules
    homology.pathway<-alter_pathway_elements(homology.table = pathway.molecules.homology.table, human.pathway = pathway.path)

    return(homology.pathway)
  }
}
####################






###########################################################################################
# identify active pathway path for RNA-seq gene expression profile
# generate a list of pathways where each sublist denote a path of the pathway

#' @title identify_active_pathway_path_new
#'
#' @description
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @rdname identify_active_pathway_path_new
#' @name identify_active_pathway_path_new
#'
#' @details
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @return This function returns a list of pathways where each sublist denote a path of the pathway.
#'
#' @param pathway.path A list with pathway path data where each sublist denotes a path of the pathway. This is used as background data.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type.
#'
#' @importFrom fastmatch fin
#'
#' @export
#'
#' @examples
#' #Pre-process the 'tooth.epi.E13.5' data
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' #Generate the homology data for the two species - 'hsapiens' and 'mmusculus'
#' mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
#' #Generate the mouse homology pathway path data
#' mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
#' #Identify active pathway paths of the processed query data
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' head(summary(tooth.epi.E13.5.active.pathway$E13.5_Epi))
#'

identify_active_pathway_path_new<-function(pathway.path, processed.query.data){
  #get all the unique pathway path genes
  pathway.path.gene<-unique(unlist(pathway.path, use.names = FALSE))


  ##for parallel processing
  #library(doParallel)
  no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  ##


  ##process separately each cell or tissue type in each iteration
  active_pathway_path<-foreach(i=1:length(processed.query.data), .packages = "fastmatch", .final = function(xx) setNames(xx, names(processed.query.data))) %dopar% {
    #get the pathway path genes that are expressed in the cell or tissue
    pathway.gene.expressed.in.cell<-intersect(pathway.path.gene, names(processed.query.data[[i]]))


    ##get the pathway paths in which all elements are expressed in query input data (i.e., in pathway.gene.expressed.in.cell)
    pathway.path.exist<-lapply(pathway.path, function(y){
      tmp.path.exist<-lapply(y, function(x){
        #for taking only the path where all molecules are expressed in gene expression data
        if(all(x %fin% pathway.gene.expressed.in.cell)){ #when x has small number of elements it is fast enough
          return(x)
        }
      })
      return(tmp.path.exist)
    })
    ##


    ##take only the existing pathway paths without null paths
    pathway.path.exist.clean<-lapply(pathway.path.exist, function(x){
      return(x[!(sapply(x, is.null))])
    })
    ##


    #filter and then return only the pathways that have at least one path
    return(pathway.path.exist.clean[lengths(pathway.path.exist.clean) > 0])
  }
  ##


  #now stop the cluster
  stopCluster(cl)


  #Finally return the whole result for active pathway path of every cell/tissue
  #by removing the cells / tissues that have no active pathways
  return(active_pathway_path[lengths(active_pathway_path) > 0])
}
####################






#######################################################################################################
## Generate pathway PPI data frame from the active pathway data to draw pathway figures using cytoscape

#' @title generate_pathway_ppi_data_frame
#'
#' @description
#' This function generates pathway PPI data frame from the active pathway data to draw pathway figures using cytoscape.
#'
#' @rdname generate_pathway_ppi_data_frame
#' @name generate_pathway_ppi_data_frame
#'
#' @details
#' This function generates pathway PPI data frame from the active pathway data to draw pathway figures using cytoscape.
#' @return This function returns pathway PPI data frame from the active pathway data to draw pathway figures using cytoscape.
#'
#' @param active.pathway.path A list of sublist containing active pathway path data for each cell / tissue.
#'
#' @export
#'
#' @examples
#' #Pre-process the 'tooth.epi.E13.5' data
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' #Generate the mouse homology pathway path data
#' mouse.homology.pathway.path<-generate_homology_pathways(species1 = "hsapiens", species2 = "mmusculus", pathway.path = pathway.path.new)
#' #Identify active pathway paths of the processed query data
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' #Generate the active pathway paths data frame
#' tooth.epi.E13.5.active.pathway.df<-generate_pathway_ppi_data_frame(active.pathway.path = tooth.epi.E13.5.active.pathway)
#' tooth.epi.E13.5.active.pathway.df[[1]][1]

generate_pathway_ppi_data_frame<-function(active.pathway.path){
  ###Each cell / tissue active pathway path data will be processed separately
  active.pathway.path.df<-lapply(active.pathway.path, function(each.cell.active.pathway.path){
    ##Here, we will get the interactions from the pathway paths
    pathway.ppi<-list()
    for(i in 1:length(each.cell.active.pathway.path)){
      tmp.individual.pathway.ppi<-NULL
      tmp.path.ppi.1<-NULL
      tmp.path.ppi.2<-NULL
      #for individual pathway at a time -  "each.cell.active.pathway.path[[i]]" denotes an individual pathway
      individual.pathway.ppi<-lapply(each.cell.active.pathway.path[[i]], function(x){
        #here, each x denotes a single path of a pathway
        tmp.path.vec<-unlist(x)
        #this loop get interactions for each path of a pathway
        for(k in 1:(length(tmp.path.vec)-1)){
          tmp.path.ppi.1<-c(tmp.path.ppi.1, tmp.path.vec[k])
          tmp.path.ppi.2<-c(tmp.path.ppi.2, tmp.path.vec[k+1])
        }
        #combine interactions for each path of a pathway
        tmp.path.ppi<-data.frame("from"=tmp.path.ppi.1, "to"=tmp.path.ppi.2)
        #combine all interactions of the paths of a pathway
        tmp.individual.pathway.ppi<-rbind(tmp.individual.pathway.ppi, tmp.path.ppi)
        return(tmp.individual.pathway.ppi)
      })
      #now combine all interactions for a pathway
      pathway.ppi[[names(each.cell.active.pathway.path)[i]]]<-do.call(rbind, lapply(individual.pathway.ppi, data.frame))
    }
    ##


    ##Now, get unique interactions from each pathway PPI data frame
    pathway.unique.ppi<-list()
    for(i in 1:length(pathway.ppi)){
      #getting the data frame for each pathway
      tmp.pathway.ppi.df<-pathway.ppi[[i]]
      #combine the neighboring factors to treat as a single vector
      tmp.pathway.ppi.df.combined<-list()
      for(j in 1:nrow(tmp.pathway.ppi.df)){
        tmp.pathway.ppi.df.combined[[j]]<-paste(tmp.pathway.ppi.df[j,1], tmp.pathway.ppi.df[j,2], sep="*")
      }
      #taking only the unique links between factors
      tmp.pathway.ppi.df.combined.unique<-unique(unlist(tmp.pathway.ppi.df.combined))
      #now separating the unique links using a list that contains all the links as vectors
      tmp.pathway.ppi.df.combined.unique.split<-lapply(tmp.pathway.ppi.df.combined.unique, function(x) {return(unlist(strsplit(x, split = "[*]")))})
      #making data frame from the unique split lists
      tmp.clean.pathway.ppi.df<-do.call(rbind, lapply(tmp.pathway.ppi.df.combined.unique.split, rbind))
      #set the column names of the data frame
      colnames(tmp.clean.pathway.ppi.df)<-c("from", "to")
      #now assign the clean ppi data frame for each pathway
      pathway.unique.ppi[[names(pathway.ppi)[i]]]<-as.data.frame(tmp.clean.pathway.ppi.df)
    }
    ##

    #it returns each cell / tissue pathway path data frame
    return(pathway.unique.ppi)
  })

  #it returns all cell / tissue pathway path data frame
  return(active.pathway.path.df)
  ###
}
####################






########################################################################################
## This function draws all the pathways and stores in one pdf file for each cell / tissue
## Draws the figures as tree layout

#' @title draw_active_pathways
#'
#' @description
#' This function draws all the pathways as tree layout and stores in one pdf file for each cell / tissue.
#'
#' @rdname draw_active_pathways
#' @name draw_active_pathways
#'
#' @details
#' This function draws all the pathways as tree layout and stores in one pdf file for each cell / tissue.
#' @return This function draws all the pathways as tree layout and stores in one pdf file for each cell / tissue and finally returns null.
#'
#' @param active.pathway.ppi.data.frame A list of sublist containing active pathway data frame for each pathway of each cell / tissue.
#'
#' @export
#'
#' @importFrom igraph graph.data.frame plot.igraph
#'
#' @examples
#' #Pre-process the 'tooth.epi.E13.5' data
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' #Generate the mouse homology pathway path data
#' mouse.homology.pathway.path<-generate_homology_pathways(species1 = "hsapiens", species2 = "mmusculus", pathway.path = pathway.path.new)
#' #Identify active pathway paths of the processed query data
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' #Generate the active pathway paths data frame
#' tooth.epi.E13.5.active.pathway.df<-generate_pathway_ppi_data_frame(active.pathway.path = tooth.epi.E13.5.active.pathway)
#' #Now draw the pathways as tree layout in one pdf file
#' draw_active_pathways(tooth.epi.E13.5.active.pathway.df)

draw_active_pathways<-function(active.pathway.ppi.data.frame){
  ##Each cell / tissue active pathway path data will be processed separately
  for(i in 1:length(active.pathway.ppi.data.frame)){
    each.cell.name<-names(active.pathway.ppi.data.frame)[i]
    each.cell.print.msg<-paste("Drawing for", each.cell.name, sep = " ")
    print(each.cell.print.msg)
    each.cell.file.name<-paste(each.cell.name, "pdf", sep = ".")
    each.cell.active.pathway.ppi.df<-active.pathway.ppi.data.frame[[i]]
    #pdf(file = each.cell.file.name, onefile = T, paper = "a4")
    pdf(file = each.cell.file.name, onefile = T)
    for(j in 1:length(each.cell.active.pathway.ppi.df)){
      tmp.g = graph.data.frame(each.cell.active.pathway.ppi.df[[j]], directed=TRUE)
      #for tree structure
      plot.igraph(tmp.g, vertex.size=5, vertex.label.cex=0.5, layout=layout_as_tree(tmp.g),
                  edge.arrow.size=0.15, edge.width=E(tmp.g)$rels/100,
                  main = names(each.cell.active.pathway.ppi.df)[j])

      #Or, for reingold structure
      #plot.igraph(tmp.g,vertex.size=5,vertex.label.cex=0.5,layout=layout.fruchterman.reingold(tmp.g,niter=10000),
      #            edge.arrow.size=0.15,edge.width=E(tmp.g)$rels/100,
      #            main = names(each.cell.active.pathway.ppi.df)[j])

      print(j)
    }
    dev.off()
  }
  return(NULL)
}
####################






#######################################################################################################
##Generate the homology gene set for species2 with respect to the list of genes of species1

#' @title generate_homology_gene_set
#'
#' @description
#' This function generates the homology gene set for species2 with respect to the list of genes of species1.
#'
#' @rdname generate_homology_gene_set
#' @name generate_homology_gene_set
#'
#' @details
#' This function generates the homology gene set for species2 with respect to the list of genes of species1.
#' @return This function returns a vector of homology gene set for species2 with respect to the list of genes of species1.
#'
#' @param species.homology.matix A sparse matrix containing the homology genes information for the two species.
#' @param species1.genes A vector containing set of genes for species1.
#'
#' @export
#'
#' @examples
#' human.zebrafish.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "drerio")
#' human.genes<-c("PAX6", "SOX2", "PROX1", "FGFR1", "BMP7", "MAPK1", "ACE2")
#' zebrafish.selected.homology.genes<-generate_homology_gene_set(species.homology.matix = human.zebrafish.homology.data[[2]], species1.genes = human.genes)
#' zebrafish.selected.homology.genes

generate_homology_gene_set<-function(species.homology.matix, species1.genes){
  #get only the species1 pathway molecules that are exist in the homology matrix
  species1.genes.filter<-unique(as.character(species1.genes[species1.genes %in% rownames(species.homology.matix)]))
  #get the homology genes for the species2 from the homology matrix
  species2.genes<-colnames(species.homology.matix)[which(col_sums(species.homology.matix[species1.genes.filter,,drop=FALSE])>0)]
  return(species2.genes)
}
####################






##########################################################################################################
# rank the active pathways based on their pathway name gene expression versus pathway expressed specific gene count proportion
# generate a list of pathways' activity score where each sublist corresponds for each cell/tissue type

#' @title get_pathway_activity_score_new
#'
#' @description
#' This function generates pathway activity score (i.e., pathway name gene expression versus pathway expressed specific gene count proportion) of the active pathways for each cell/tissue type. It uses active pathway path and processed query data with the homology information to generate the activity score. It also uses housekeeping gene data from SPAGI package to get the cell/tissue expressed specific pathway gene count proportion.
#'
#' @rdname get_pathway_activity_score_new
#' @name get_pathway_activity_score_new
#'
#' @details
#' This function generates pathway activity score (i.e., pathway name gene expression versus pathway expressed specific gene count proportion) of the active pathways for each cell/tissue type. It uses active pathway path and processed query data with the homology information to generate the activity score. It also uses housekeeping gene data from SPAGI package to get the cell/tissue expressed specific pathway gene count proportion.
#'
#' @return This function returns a list of pathway activity score for each cell/tissue type.
#'
#' @param active.pathway.path A list of active pathway path data for each cell/tissue type as returned by the function 'identify_active_pathway_path'.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type as returned by the function 'preprocess_querydata'.
#' @param homology.data A list containing the homology data (i.e., homology table, homology matrix and species2 name) for the two species.
#' @param hk.gene A vector consisting of housekeeping genes. Default is housekeeping.gene. This data is loaded automatically with the SPAGI package. This data was generated using the gene expression profiles of different cell types and/or tissues from the ENCODE human and mouse project.
#'
#' @export
#'
#' @examples
#' #Pre-process the 'tooth.epi.E13.5' data
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' #Generate the homology data for the two species - 'hsapiens' and 'mmusculus'
#' mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
#' #Generate the mouse homology pathway path data
#' mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
#' #Identify active pathway paths of the processed query data
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' #Get pathway activity score (i.e., pathway name gene expression and pathway specific gene count proportion) of the processed query data
#' tooth.epi.E13.5.active.pathway.score<-get_pathway_activity_score_new(active.pathway.path = tooth.epi.E13.5.active.pathway, processed.query.data = tooth.epi.E13.5.processed.data, homology.data = mouse.homology.data)
#' head(summary(tooth.epi.E13.5.active.pathway.score$E13.5_Epi))
#'

get_pathway_activity_score_new<-function(active.pathway.path, processed.query.data, homology.data, hk.gene = housekeeping.gene){
  ##first check for species, if species='hsapiens', keep the housekeeping genes as it was, otherwise get homology housekeeping genes
  if(homology.data$species=="hsapiens"){
    homology.hk.gene<-hk.gene
  }else{
    homology.hk.gene<-generate_homology_gene_set(homology.data[[2]], hk.gene)
  }
  ##


  ##process separately each cell/tissue to get active pathway ranking metric
  cell.pathway.activity.score<-list()
  pathway.activity.score<-list()
  for(i in 1:length(active.pathway.path)){
    #get each cell/tissue active pathway paths
    each.cell.active.pathway.path<-active.pathway.path[[i]]

    #take the cell/tissue name from the pathway to get that cell/tissue processed.query.data
    tmp.cell.name<-names(active.pathway.path[i])

    #take the respective cell/tissue processed data
    tmp.cell.processed.data<-processed.query.data[[tmp.cell.name]]


    ##get 1st metric - expression for the pathway name genes, i.e., expression of the receptors for the cell/tissue
    #get the names of all pathways for the cell/tissue
    each.cell.active.pathway.names<-names(each.cell.active.pathway.path)
    #get the expression of all active pathway names for the cell/tissue
    tmp.pathway.genes.exp<-tmp.cell.processed.data[each.cell.active.pathway.names]
    #assign expression for the active pathway genes
    pathway.activity.score[["pathway.name.gene.exp"]]<-tmp.pathway.genes.exp
    ##


    ##get 2nd metric - specifically expressed genes count proportion for each pathway for the cell/tissue
    pathway.path.specific.gene.count.proportion<-lapply(each.cell.active.pathway.path, function(y){
      #get each pathway specifically expressed unique set of genes that are not present in the housekeeping genes
      tmp.pathway.spec.genes<-setdiff(unique(unlist(y)), homology.hk.gene)
      #calcaulate top expressed gene count proportion for each pathway and return
      return(length(tmp.pathway.spec.genes) / length(unique(unlist(y))))
    })
    #assign pathway expressed specific gene count proportion
    pathway.activity.score[["pathway.spec.gene.cout.prop"]]<-unlist(pathway.path.specific.gene.count.proportion)
    ##


    #now assign pathways RP expression and specifically expressed gene count proportion for each cell/tissue
    cell.pathway.activity.score[[tmp.cell.name]]<-pathway.activity.score
  }
  ##


  #return pathways activity score for all cell/tissue
  return(cell.pathway.activity.score)
}
####################






################################################################################################
# plot the pathway activity score values in a 2D plane for each cell type or tissue

#' @title display_pathway_activity_score
#'
#' @description
#' This function plots the pathway activity score values for each cell type or tissue in a 2D plane. This function uses generally highly expressed receptor data from the package to separate the cell/tissue specific active pathways and generally active pathways in many cells. In the figure, the black and gray colours represent cell/tissue specifically active pathways and generally active pathways in many cells respectively.
#'
#' @rdname display_pathway_activity_score
#' @name display_pathway_activity_score
#'
#' @details
#' This function plots the pathway activity score values for each cell type or tissue in a 2D plane. This function uses generally highly expressed receptor data from the package to separate the cell/tissue specific active pathways and generally active pathways in many cells. In the figure, the black and gray colours represent cell/tissue specifically active pathways and generally active pathways in many cells respectively.
#'
#' @return NULL
#'
#' @param pathway.activity.score The ranking metric result returned by 'get_pathway_ranking_metric' function.
#' @param homology.data A list containing the homology data (i.e., homology table, homology matrix and species2 name) for the two species.
#' @param highly.expressed.rp A vector containing a list of receptor proteins that are generally highly expressed in many cells and/or tissues. Default is rp.median.exp.4.
#' @param rp.expression.range A vector containing an expression range (in the format of c(min, max)) of the receptor proteins to show at the y axis.
#'
#' @export
#'
#' @examples
#' #Pre-process the 'tooth.epi.E13.5' data
#' tooth.epi.E13.5.processed.data<-preprocess_querydata_new(cell.tissue.data = tooth.epi.E13.5, exp.cutoff.th = 5.0, species="mmusculus")
#' #Generate the homology data for the two species - 'hsapiens' and 'mmusculus'
#' mouse.homology.data<-generate_homology_data(species1 = "hsapiens", species2 = "mmusculus")
#' #Generate the mouse homology pathway path data
#' mouse.homology.pathway.path<-generate_homology_pathways(species.homology.data = mouse.homology.data, pathway.path = pathway.path.new)
#' #Identify active pathway paths of the processed query data
#' tooth.epi.E13.5.active.pathway<-identify_active_pathway_path_new(pathway.path = mouse.homology.pathway.path, processed.query.data = tooth.epi.E13.5.processed.data)
#' #Get pathway activity score (i.e., pathway name gene expression and pathway specific gene count proportion) of the processed query data
#' tooth.epi.E13.5.active.pathway.score<-get_pathway_activity_score_new(active.pathway.path = tooth.epi.E13.5.active.pathway, processed.query.data = tooth.epi.E13.5.processed.data, homology.data = mouse.homology.data)
#' #Plot the activity score result (i.e., pathway name gene expression versus pathway specific gene count proportion) in a 2D plane (black=specifically active, gray=generally active)
#' display_pathway_activity_score(pathway.activity.score = tooth.epi.E13.5.active.pathway.score, homology.data = mouse.homology.data)
#' #To separate the top ranked pathways we can do this
#' abline(v=0.2, h=10, lty=2, col="black")
#'

display_pathway_activity_score<-function(pathway.activity.score, homology.data, highly.expressed.rp = rp.median.exp.4, rp.expression.range = c(5,15)){
  ##first check for species, if species='hsapiens', keep the highly.expressed.rp as it was, otherwise get homology housekeeping genes
  if(homology.data$species=="hsapiens"){
    homology.highly.expressed.rp<-highly.expressed.rp
  }else{
    homology.highly.expressed.rp<-generate_homology_gene_set(homology.data[[2]], highly.expressed.rp)
  }


  ##in each loop, the result of one query cell/tissue type is processed and plotted in a 2D plane
  for(i in 1:length(pathway.activity.score)){
    #get the name of each cell type
    cell.tissue.names<-names(pathway.activity.score)[i]


    #for setting the title
    if(nchar(cell.tissue.names)>40){
      title<-paste("The result of:\n", cell.tissue.names, sep = " ")
    }
    else{
      title<-paste("The result of", cell.tissue.names, sep = " ")
    }


    #first take the pathway names for each cell/tissue
    tmp.pathway.names<-names(pathway.activity.score[[i]]$pathway.spec.gene.cout.prop)
    #get the generally higly expressed rp indices
    tmp.highly.expressed.rp.ind<-which(tmp.pathway.names %in% homology.highly.expressed.rp)
    #get the specifically expressed rp indices
    tmp.specifically.expressed.rp.ind<-which(!(tmp.pathway.names %in% homology.highly.expressed.rp))

    devAskNewPage(ask = FALSE) #to turn off the "Hit <Return> to see next plot" prompt

    #plot the result in a 2D plane - specific gene count proportion in x axis and RP expression in y axis
    plot(pathway.activity.score[[i]]$pathway.spec.gene.cout.prop, pathway.activity.score[[i]]$pathway.name.gene.exp,
         xlim = c(0,1), ylim = rp.expression.range, type= "n", bty="n", main = title,
         #xlim = c(0,1), ylim = c(2,10), type= "n", bty="n", main = title,
         #xlim = c(0,1), ylim = c(3,15), type= "n", bty="n", main = title,
         xlab = "Pathway expressed specific gene count proportion", ylab = "Pathway name gene expression")
    #print the specifically active pathway names with black colour
    #check for non-zero length tmp.specifically.expressed.rp.ind first and then text - when all rps are highly expressed
    if(length(tmp.specifically.expressed.rp.ind) != 0){
      text(pathway.activity.score[[i]]$pathway.spec.gene.cout.prop[tmp.specifically.expressed.rp.ind], pathway.activity.score[[i]]$pathway.name.gene.exp[tmp.specifically.expressed.rp.ind],
           tmp.pathway.names[tmp.specifically.expressed.rp.ind], cex = 0.5, col ="black")
    }
    #print the generally highly expressed active pathway names with gray colour
    #check for non-zero length tmp.highly.expressed.rp.ind first and then text - when all rps are specifically expressed
    if(length(tmp.highly.expressed.rp.ind) != 0){
      text(pathway.activity.score[[i]]$pathway.spec.gene.cout.prop[tmp.highly.expressed.rp.ind], pathway.activity.score[[i]]$pathway.name.gene.exp[tmp.highly.expressed.rp.ind],
           tmp.pathway.names[tmp.highly.expressed.rp.ind], cex = 0.5, col ="gray")
    }


    #for console message for each cell/tissue type ploting
    consol.msg<-paste(cell.tissue.names, "-- result plotting done!!", sep = " ")
    print(consol.msg)
  }
}
####################



