##########
##This file should be run in the first.
##This file contains necessary code to generate pathway molecules GO mapping information in a matrix format.
##The input data sets for this file are list of pathway molecules that are deposited in 'additional_file' directory of the SPAGI2 package repository.
##First create a folder to store the "molecule_go_map.R" file. In this folder also create a folder and give name "molecule_object". 
##In that "molecule_object" folder store the data sets of RData types - "LG.gene.list", "RP.gene.list", "KN.gene.list", "TF.gene.list".
##Then, install (if already not) load the required package "biomaRt" and run the code for the global variables.
##Next, run all the functions and the run the code in the main part section.
##From this file, save the RData files - 'genes.go.mat', 'hs.all.gene.go.mat', 'hs.ensembl.symbol.map' and 'molecule.ensembl.id' in the result folder of the current directory.
##These RData files will be used in the next files.
##########



#####
#If 'BiocManager' is not already installed, first install 'BiocManager' to insall the bioconductor packages.
#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#then install the biomaRt package
#BiocManager::install("biomaRt")
require(biomaRt)
#####



#####run at the start of the program and used as global
##When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed, 
##for example to an archived version
biomart.ID <- "ENSEMBL_MART_ENSEMBL"
#host.ID <-  "www.ensembl.org"
#host.ID <- "asia.ensembl.org"
host.ID <- "http://sep2019.archive.ensembl.org"
#####



#Sometimes you might get temporarily disconnected from the BioMart web server with an error message like:
#
#  Error in bmRequest(request = request, verbose = verbose) : 
#    Internal Server Error (HTTP 500). 
#
# If this happens just try to re assign the variable 'host.ID' with an archived version. 
# For example -
# host.ID <- "http://sep2019.archive.ensembl.org"
# 
# Alternatively the BioMart web service is temporarily down. Just wait for the BioMart web service up again.





#####to combine the molecules and as well as to create their types
#####here, duplicate molecules will retain the first type
combine_pathway_molecules<-function(ligand, receptor, kinase, transfactor){
  #create the types with numbers by repeating
  lg.type<-rep(1,length(ligand))
  rp.type<-rep(2,length(receptor))
  kn.type<-rep(3,length(kinase))
  tf.type<-rep(4,length(transfactor))
  
  #make the molecule type as vector and assign names with the gene names
  molecule.type.data<-c(lg.type, rp.type, kn.type, tf.type)
  molecule.names<-c(ligand, receptor, kinase, transfactor)
  names(molecule.type.data)<-molecule.names
  
  #getting only the first molecules by removing the duplicates
  molecule.type.data.clean<-molecule.type.data[!duplicated(names(molecule.type.data))]
  
  return(molecule.type.data.clean)
}
#####





#####get ensembl symbol table - this will be used to finally convert pathway molecules' ensembl id to gene names
get_ensembl_symbol_map <- function(species="hsapiens"){
  #assignmet of values to get information from ensembl database
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  
  #get list of genes with ensembl ids
  ensembl.gene <- getBM(mart=ensembl, attributes=c('ensembl_gene_id','external_gene_name'))
  
  #get the genes that have both ensembl id and name
  ensembl.gene.exist<-ensembl.gene[!(ensembl.gene$ensembl_gene_id=="" | ensembl.gene$external_gene_name==""),]
  
  return(ensembl.gene.exist)
}
#####





#####get ensembl gene information
get_ensembl_gene <- function(species="hsapiens"){
  #assignmet of values to get information from ensembl database
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  
  #get list of genes with ensembl ids
  ensembl.gene <- getBM(mart=ensembl, attributes=c('ensembl_gene_id','external_gene_name', 'description', 'go_id', 'name_1006'))
  
  #get the genes that have description
  ensembl.gene.description.exist<-ensembl.gene[!(ensembl.gene$description=="" |
                                                   ensembl.gene$external_gene_name=="" |
                                                   ensembl.gene$ensembl_gene_id=="" | 
                                                   ensembl.gene$go_id =="" |
                                                   ensembl.gene$name_1006==""),]
  
  #get the genes based on the unique names
  ensembl.gene.unique<-unique(ensembl.gene.description.exist$external_gene_name)
  
  return(ensembl.gene.unique)
}
#####





#####For other type of random genes and combining with the pathway molecules
combine_pathway_random_genes<-function(pathway.gene.type, ensembl.gene){
  #get only the other genes without the pathway genes (factor)
  ensembl.gene.without.factor<-ensembl.gene[!(ensembl.gene %in% names(pathway.gene.type))]
  
  #geting the random 2500 other genes
  other.gene<-sample(ensembl.gene.without.factor, 2500, replace = F)
  
  #set the type as '5' of the other 2500 genes
  other.gene.type<-rep(5, length(other.gene))
  names(other.gene.type)<-other.gene
  
  #combining the pathway genes (molecules) and the other genes
  pathway.other.gene.type<-c(pathway.gene.type, other.gene.type)
  
  return(pathway.other.gene.type)
}
#####





#############################################################################################
#####to get GO annotation information for the genes - convert the external gene names to ensembl gene ids
#####here used the ensembl gene ids finally as gene ids as RF uses data frame as model data, and 
#####data frame changes gene names, e.g., LFA-G to LFA.G 
get_go_mapping<-function(gene.list, species = "hsapiens", gene.id = "external_gene_name"){
  #assignmet of values to get information from ensembl database
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  
  
  ##check and convert for ensembl gene id
  if(gene.id == "ensembl_gene_id"){
    gene.list.ensembl<-gene.list
  }else if(gene.id == "external_gene_name"){
    ensembl.symbol.map <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
    
    #getting unique gene names 
    list.gene.names<-unique(names(gene.list))
    
    #get only ensembl mapping genes that exist in user provided gene set
    ensembl.genes.exist.in.user.cell.data <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% list.gene.names,]
    
    #!!some gene symbols have duplicate ensembl ids, e.g., DDR1 has 'ENSG00000137332', 'ENSG00000215522', 'ENSG00000230456', etc
    #to take only one ensembl id against a gene symbol - here we will take the first one match in the mapping
    duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.cell.data[[2]])
    ensembl.genes.exist.in.user.cell.data.without.duplication <- ensembl.genes.exist.in.user.cell.data[!duplicate.ensembl.gene.list,]
    #ens.gene.id<- ensembl.genes.exist.in.user.cell.data.without.duplication[[1]] #for only ensembl gene ids
    
    #take only the user provided genes that have mapping after removal of duplication of ensembl ids
    gene.list.ensembl<- gene.list[names(gene.list) %in% ensembl.genes.exist.in.user.cell.data.without.duplication[[2]]]
    
    #finally replace the names of the gene names by the corresponding ensembl gene ids
    names(gene.list.ensembl)<-ensembl.genes.exist.in.user.cell.data.without.duplication[[1]][match(names(gene.list.ensembl),ensembl.genes.exist.in.user.cell.data.without.duplication[[2]])]
  }else{
    print("Do not support other gene IDs at this moment!")
    return(NULL)
  }
  ##
  
  #getting the list of ensembl gene ids
  ens.gene.id<-names(gene.list.ensembl)
  
  #genes mapping for description and GO annotation
  genes.mapping<-getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006'), 
                       filters = 'ensembl_gene_id', 
                       values = ens.gene.id, 
                       mart = ensembl)
  
  #clean the genes description part
  genes.mapping$description<-unlist(lapply(genes.mapping$description,function(x){return(strsplit(x,' [[]')[[1]][1])}))
  
  #removing the records that have null go id or ensembl_gene_id
  #genes.mapping<-genes.mapping[!(genes.mapping$go_id==""),]    #only for go_id
  genes.mapping<-genes.mapping[!(genes.mapping$go_id=="" | genes.mapping$ensembl_gene_id==""),]
  
  
  #combining the 'gene.list.ensembl' and 'genes.mapping' using a list and return
  gene.list.information<-list()
  gene.list.information[["ensembl.gene.id.type"]]<-gene.list.ensembl
  gene.list.information[["go.mapping"]]<-genes.mapping
  
  return(gene.list.information)
}
#################################################################################################






#######################################################################
#####separation of each gene annotation information and store in a list
separate_go_annotation_information<-function(genes.mapping){
  individual.gene.inf<-list()
  mapping.genes.ensembl.id<-unique(genes.mapping$ensembl_gene_id)
  for(i in 1:length(mapping.genes.ensembl.id)){
    #getting all records of the selected gene 
    gene.inf.ind<-which(genes.mapping$ensembl_gene_id==mapping.genes.ensembl.id[i])
    selected.gene.inf<-genes.mapping[gene.inf.ind,]
    
    #assigning the information to each gene
    mapping<-list()
    mapping[["ensembl.gene.id"]]<-unique(selected.gene.inf$ensembl_gene_id)
    mapping[["description"]]<-unique(selected.gene.inf$description)
    mapping[["go.id"]]<-selected.gene.inf$go_id
    mapping[["go.term.name"]]<-selected.gene.inf$name_1006
    
    individual.gene.inf[[mapping.genes.ensembl.id[i]]]<-mapping
  }
  return(individual.gene.inf)
}
######################################################################






#########################################################################
##Make a matrix based on the information of each gene
##Here, used both the genes.mapping$go_id and individual.gene.inf
make_matrix_using_go_information<-function(genes.mapping.go.id, individual.gene.inf){
  mapping.mat<-matrix(0, nrow = length(individual.gene.inf), ncol = length(unique(genes.mapping.go.id)))
  colnames(mapping.mat)<-sort(unique(genes.mapping.go.id))
  rownames(mapping.mat)<-names(individual.gene.inf)
  
  #assign 1 for each gene's go ids, i.e., each gene's corresponding column position of go ids
  for(i in 1:nrow(mapping.mat)){
    #gene.go.ids<-individual.gene.inf[[i]]$go.id
    #mapping.mat[i,gene.go.ids]<-1
    mapping.mat[i,individual.gene.inf[[i]]$go.id]<-1
  }
  
  return(mapping.mat)
}
#########################################################################






################################################################################
##get the GO information for each gene with the type
get_genes_go_matrix<-function(gene.type.list){
  #get the go mapping information of the genes
  gene.information<-get_go_mapping(gene.list = gene.type.list, species = "hsapiens")
  
  #separate go mapping information for each gene 
  gene.information.individual<-separate_go_annotation_information(gene.information$go.mapping)
  
  #make a data frame using go mapping go ids and each gene information
  gene.mat<-make_matrix_using_go_information(gene.information$go.mapping$go_id, gene.information.individual)
  
  #combining the genes' type and go annotation information
  gene.with.type<-gene.information$ensembl.gene.id.type
  #taking only the genes from 'gene.with.type' that have mapping GO, i.e., in the 'gene.mat' object
  #also matching the type of each molecule in the same sequence according to the 'mapping.df'
  gene.go.mapping.exist<-gene.with.type[match(rownames(gene.mat), names(gene.with.type))]
  gene.mat.with.type<-cbind("type"=gene.go.mapping.exist, gene.mat)
  
  return(gene.mat.with.type)
}
#####################################################################################







###############################################################################################
####################################Main part##################################################
##load the molecules
load("molecule_object/LG.gene.list.RData")
load("molecule_object/RP.gene.list.RData")
load("molecule_object/KN.gene.list.RData")
load("molecule_object/TF.gene.list.RData")
##


#####
##get the go matrix information for the pathway molecule and other type ensembl genes
#combine the pathway genes (molecules) with their types
molecule.type<-combine_pathway_molecules(LG.gene.list, RP.gene.list, KN.gene.list, TF.gene.list)

#get ensembl unique gene lists for hsapiens
hs.ensembl.gene<-get_ensembl_gene()

#get randomly other 2500 genes and then combine these with the pathway genes
pathway.random.gene.type<-combine_pathway_random_genes(molecule.type, hs.ensembl.gene)

#get the matrix for machine learning classification model - 1st time data
#here, the gene symbols are converted to ensembl gene ids
#in classification, the matrix with rownames as ensembl gene ids will be used
genes.go.mat<-get_genes_go_matrix(pathway.random.gene.type)
dim(genes.go.mat)
#[1]  5245 12515
#[1]  5209 12677 with archived one
length(which(genes.go.mat[1,]==1))
#[1] 15
#save(genes.go.mat, file = "result/genes.go.mat.RData")
##


##Now get the go matrix information for all hsapiens ensembl genes
#set the type as '5' by default of the 19989 genes
hs.all.gene<-rep(5, length(hs.ensembl.gene))
names(hs.all.gene)<-hs.ensembl.gene

#get the matrix for machine learning classification model - 2nd time data
hs.all.gene.go.mat<-get_genes_go_matrix(hs.all.gene)
dim(hs.all.gene.go.mat)
#[1] 19939 18428
#[1] 19685 18398 with archived one
#save(hs.all.gene.go.mat, file = "result/hs.all.gene.go.mat.RData")
##
#####



#####the data for this section is used for classification model section
##get ensembl symbol map for hsapiens
##this will be used to finally convert pathway molecules' ensembl id to gene names
hs.ensembl.symbol.map<-get_ensembl_symbol_map()
#save(hs.ensembl.symbol.map, file="result/hs.ensembl.symbol.map.RData")
##


##to convert pathway molecules' names to ensembl gene ids
##get information from the functions
molecule.data.information<-get_go_mapping(gene.list = molecule.type)
molecule.ensembl.id<-unique(names(molecule.data.information$ensembl.gene.id.type))
#save(molecule.ensembl.id, file = "result/molecule.ensembl.id.RData")
##
#####
###############################################################################################


