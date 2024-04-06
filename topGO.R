#### topGO enrichment analyses ####
#https://zhiganglu.com/post/topgo-ks-test/

#AIM: conduct an enrichment analysis to test for the overreprentation of functional categories in a list of differentially expressed genes.

# as input you need a test dataset = list of IDs corresponding to the differentially expressed genes, 
# and as reference set you need a table containing the IDs or expressed genes plus their corresponding GO IDs (functions), e.g. obtained from Interproscan

#load the following libraries
BiocManager::install("topGO")
library(topGO)

BiocManager::install("AnnotationHub")
library(AnnotationHub)

#install packages 'dbplyr' and 's4Vectors'
install.packages("dbplyr")
library("dbplyr")

BiocManager::install("S4Vectors", force = TRUE)
library("s4Vectors")

#set working directory = full path to directory containing your input files
setwd("C:/Users/kennedy/Desktop/Thesis/chapter 3")

#Load the annotation hub and check for the presence of the honey bee among the annotation data
hub <- AnnotationHub()
# previous line returns a question asking to create a directory (only for the first time use)
yes
#
query(hub, "Apis mellifera")

#Needed is: org.Apis_mellifera.eg.sqlite, retrieve with AH86655
#Amel <- hub[["AH86655"]]

#when updating R maybe use to retrieve
Amel <- hub[["AH102515"]]

#suggested retrieval using 'org.Apis_mellifera.eg.sqlite' 
#Amel <- hub[["AH97500"]]
Amel #Gives you version numbers and other info

#Retrieve the GO annotations
#The "all_geneids.txt" file contains the first column of the htseq output file (minus the stats at the bottom)
#this then lists all the possible geneids that can be in your DEG sets for the background
allkeys <- read.table("all_geneids.csv", sep=";", header=TRUE)
gene_to_go <- select(Amel, allkeys$GeneID, columns = c("GO"), "ALIAS")
#N.B.: If the select function gives an error, loaded R packages may be conflicting
#Restart R to clear memory!
#in order to fix the error message:
#I had to install packages 'dbplyr' and 's4Vectors'

#Remove genes without GO annotations
go_data <- na.omit(gene_to_go)

## Aggregate the list, which concatenates the second column:
mymerge <- function(x){
  all_in_one <- paste(unlist(x),
                      sep = ",",
                      collapse = ",")
  split_term <- unlist(strsplit(all_in_one,
                                split = ","))
  return(paste(unique(split_term),
               sep = ",",
               collapse = ","))
}

## Run function:
output <- aggregate(go_data[-1],
                    by = list(go_data$ALIAS), mymerge)

write.table(output,
            file = "go_term_database.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

#I included also these:
#"NurFor24dayMBForControlvsNurForQMP+.txt"
#"NewlyEmerged_All_Tissue_DEG_Overlap.txt"

#Load gene to go mapping
gene_to_go_mapping <- readMappings("go_term_database.txt")
geneUniverse <- names(gene_to_go_mapping)

#Set node size and output directory, output directory should be manually created before running the loop
node_size <- 5
output_directory <- "Results"

fileList <- c("ALNewForQMP-vsNewForQMP+.txt", "ALNurForQMP-vsNurForQMP+.txt", "ATNewForQMP-vsNewForQMP+.txt",
              "ATNurForQMP-vsNurForQMP+.txt", "Foragers21dayAL.txt", "Foragers21dayAT.txt", "Foragers21dayMB.txt",
              "MBNewForQMP-vsNewForQMP+.txt", "NurFor24dayMBForControlvsNurForQMP+.txt", "NewlyEmerged_All_Tissue_DEG_Overlap.txt")

for (file_name in fileList) {
  name <- sub('\\.txt$', '', file_name)
  degs <- read.table(file_name, sep = "\t", header = TRUE)
  genesOfInterest <- as.character(degs$Geneid)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  for (go_category in c("BP", "MF", "CC")) {
    ## Build the GOdata object in topGO
    my_go_data <- new("topGOdata",
                      description = paste("GOtest", go_category, sep = "_"),
                      ontology    = go_category,
                      allGenes    = geneList,
                      gene2GO     = gene_to_go_mapping,
                      annot       = annFUN.gene2GO,
                      nodeSize    = node_size) # Modify to reduce/increase stringency.
    ## Calculate ks test using 'weight01' algorithm:
    result_weight_ks <- runTest(object    = my_go_data,
                                algorithm = "weight01",
                                statistic = "ks")
    ## Calculate fisher exact test using 'weight01' algorithm:
    result_weight_fisher <- runTest(object    = my_go_data,
                                    algorithm = "weight01",
                                    statistic = "fisher")
    ## Combine results from statistical tests:
    result_weight_output <- GenTable(object       = my_go_data,
                                     weight_ks     = result_weight_ks,
                                     weight_fisher = result_weight_fisher,
                                     orderBy       = "weight_fisher",
                                     topNodes      = length(score(result_weight_fisher)))
    ## Correct ks test for multiple testing:
    result_weight_output$weight_ks <- as.numeric(result_weight_output$weight_ks)
    result_weight_output$weight_fisher <- as.numeric(result_weight_output$weight_fisher)
    result_weight_output$weight_ks_adjusted <- p.adjust(p = result_weight_output$weight_ks,
                                                        method = c("BH"))
    result_weight_output$weight_fisher_adjusted <- p.adjust(p = result_weight_output$weight_fisher,
                                                            method = c("BH"))
    ## Subset calls with significance higher than expected:
    result_weight_output_sig <- subset(result_weight_output,
                                       subset = (Significant > Expected) &
                                         (weight_fisher_adjusted < 0.05))
    ## Write to output:
    write.table(result_weight_output_sig,
                file = file.path(output_directory,
                                 paste(name,go_category,
                                       "sig.tsv",
                                       sep = "_")),
                row.names = FALSE,
                sep       = "\t",
                quote = FALSE)
  }
}

#"EOF within quoted string" warning can be ignored
