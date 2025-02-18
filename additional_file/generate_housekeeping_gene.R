############################
##This file can be run at any position as it is not linked with other files. However, the other files should be run sequentially.
##This file is used to generate the housekeeping genes
##It utilizes both the human and mouse encode data sets with their descriptors.
##Here, a gene is treated as housekeeping if it expresses >=75% of cells/tissues
##Both for human and mouse encode data sets housekeeping genes are generated separately
##Then, an unique set of housekeeping genes is generated by taking union of the two sets of housekeeping genes. 
############################





############################################
# If 'BiocManager' is not already installed, first install 'BiocManager' to insall the bioconductor packages. Type the following in an R command window:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
#
# Finally load the library with
library(biomaRt)
#############################################
#
#
#
####################
##!!When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed
biomart.ID <- "ENSEMBL_MART_ENSEMBL"
host.ID <-  "www.ensembl.org"
##or,
#host.ID <- "asia.ensembl.org"
##or, for example to an archived version
#host.ID <- "http://sep2019.archive.ensembl.org"
####################






#############################
# Returns a data frame with two columns containing all gene symbols and their corresponding ensembl gene IDs for the given species
# Here, geneID.forMapping is used for the existing gene ID (that has external_gene_name mapping in the biomaRt database), such as "ensembl_gene_id".
make_ENSEMBL_symbol_map <- function(species, ensemblId.toMapping, geneID.forMapping){
  dataset.name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart = biomart.ID, dataset = dataset.name, host = host.ID)
  ENSEMBL.symbol.map <- getBM(attributes=c(ensemblId.toMapping, geneID.forMapping), mart = ensembl)
  return(ENSEMBL.symbol.map)
}
##############################





#########################################
##This function converts a vector of ensembl gene IDs to gene symbols 
convert_ensemblID_to_geneSymbol <- function(user.genes, species="hsapiens", ensemblId.toMapping="external_gene_name", geneID.forMapping = "ensembl_gene_id"){
  ensembl.symbol.map<-make_ENSEMBL_symbol_map(species, ensemblId.toMapping, geneID.forMapping)
  ensembl.genes.exist.in.user.genes <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% user.genes,]
  duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.genes[[2]])
  ensembl.genes.exist.in.user.genes.without.duplication <- ensembl.genes.exist.in.user.genes[!duplicate.ensembl.gene.list,]
  user.genes.symbol<- ensembl.genes.exist.in.user.genes.without.duplication[[1]]
  return(user.genes.symbol)
}
#########################################






#############################
##This funtion computes mean for a data matrix.
##First checks for vector data and if it is then simply make a data matrix by copying the values.
compute_mean <- function(data.matrix){
  if((is.vector(data.matrix)==TRUE) || (ncol(data.matrix)<2))
    data.matrix<-cbind(val1=data.matrix,vall2=data.matrix)
  return(rowMeans(data.matrix))
}
#############################





#####################
##This function formats the list data to make a matrix formatted data
format_list_data<-function(list.data, experiment.descriptor){
  names(list.data)<-experiment.descriptor
  mat.data<-do.call(cbind,lapply(list.data,compute_mean))
  return(mat.data)
}
#####################





######################################
##look for the distribution for all samples
look_for_distributin<-function(dd.matrix){
  cell.tissue.names<-colnames(dd.matrix)
  for(i in 1:ncol(dd.matrix)){
    title.text<-paste("Distribution of", cell.tissue.names[i], sep = " ")
    hist(dd.matrix[,i], breaks=100, main=title.text, xlim=c(0,10), ylim=c(0,2000), col="grey", xlab = "log2(FPKM + 1)", ylab = "Number of genes")
    abline(v=1.5,col="red")
    Sys.sleep(1)
  }
  return(NULL)
}
##
####################################





########################################
##find out the expressed genes: log2(value)>1.5
find_expressed_genes<-function(exp.data.log, th=1.5){
  all.expressed.data<-list()
  ct.names<-colnames(exp.data.log)
  for(i in 1:ncol(exp.data.log)){
    individual.data<-exp.data.log[,i]
    names(individual.data)<-rownames(exp.data.log)
    
    l<-lapply(individual.data,function(x){
      length(which(x>=th))  ### use the number as greater than abline used number i.e., v=1.5 
    })
    
    #find out expressed genes in all samples that satisfy the condition
    individual.data.expressed<-individual.data[which(l>=1)]
    all.expressed.data[[ct.names[i]]]<-individual.data.expressed
  }
  return(all.expressed.data)
}
##
############################################







#################################################################################################
############################################ For human compendium data ###########################
###############pre-processing the compendium data
##load human encode data
load("objects/all_FPKMs_human_cleanID.RData")
load("objects/experiment.descriptors.humanEncode.RData") 

##Make the human compendium
human.compendium.data<-format_list_data(all_FPKMs_clean_ID, experiment.descriptors.humanEncode)
dim(human.compendium.data)
#[1] 57820   144
##

##taking log2 of the data
human.compendium.data.log<-log2(human.compendium.data+1)
################



###########################
##look for the distribution
look_for_distributin(dd.matrix = human.compendium.data.log)

##Find out the expressed genes from the encode data set
human.expressed.data<-find_expressed_genes(exp.data.log = human.compendium.data.log, th = 1.5)

##taking all the expressed ensembl ids for all cells/tissues 
human.expressed.ensembl.id <- lapply(human.expressed.data, function(X){return(names(X))})

##get the frequency for each expressed genes
human.expressed.ensembl.id.frequency<-table(unlist(human.expressed.ensembl.id))

##get the HK genes that are expressed more than 75% samples, out of 144 is 108
human.hk.ensembl.id<-human.expressed.ensembl.id.frequency[which(human.expressed.ensembl.id.frequency >= 108)]
length(human.hk.ensembl.id)
#[1] 9439

##convert ensembl ids to gene names
human.hk.gene<-convert_ensemblID_to_geneSymbol(names(human.hk.ensembl.id))
############################
#############################################################################################








#################################################################################################
############################################ For mouse compendium data ###########################
###############pre-processing the compendium data
##load mouse encode data
load("objects/all_FPKMs_Mouse_cleanID.RData")
load("objects/experiment.descriptors.mouseEncode.RData") 

##Make the mouse compendium
mouse.compendium.data<-format_list_data(all_FPKMs_clean_ID_Mouse, experiment.descriptors.mouseEncode)
dim(mouse.compendium.data)
#[1] 43346    94
##

##taking log2 of the data
mouse.compendium.data.log<-log2(mouse.compendium.data+1)
################



###########################
##look for the distribution
look_for_distributin(dd.matrix = mouse.compendium.data.log)

##Find out the expressed genes from the encode data set
mouse.expressed.data<-find_expressed_genes(exp.data.log = mouse.compendium.data.log, th = 1.5)

##taking all the expressed ensembl ids for all cells/tissues 
mouse.expressed.ensembl.id <- lapply(mouse.expressed.data, function(X){return(names(X))})

##get the frequency for each expressed genes
mouse.expressed.ensembl.id.frequency<-table(unlist(mouse.expressed.ensembl.id))

##get the HK genes that are expressed more than 75% samples, out of 94 is 70.5
mouse.hk.ensembl.id<-mouse.expressed.ensembl.id.frequency[which(mouse.expressed.ensembl.id.frequency >= 70.5)]
length(mouse.hk.ensembl.id)
#[1] 9476

##convert ensembl ids to gene names
mouse.hk.gene<-convert_ensemblID_to_geneSymbol(names(mouse.hk.ensembl.id), species = "mmusculus")
############################
#############################################################################################







######################################################
##make a unique list of housekeeping gene to use for SPAGI
housekeeping.gene<-unique(c(toupper(mouse.hk.gene), human.hk.gene))
#save(housekeeping.gene, file = "results/housekeeping.gene.RData")
######################################################
