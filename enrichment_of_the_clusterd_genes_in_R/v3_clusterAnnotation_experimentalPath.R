setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/")
dir.create( recursive = TRUE,"./R/Enrichment")
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/")

#Copy all images ending with, txt.
system("cp /mnt/haus/marius/SingleCell10xKallistoJupyter/JupyterNotebooks/*.txt /mnt/haus/marius/SingleCell10xKallistoJupyter/R/")
files <- list.files(all.files = TRUE,pattern = "txt")
files


###########
#         #
#   PRE   #
#         #
###########


#Preprocess the data frame to follow the format of each column a cluster and each column will pass to be a list

#files[grepl(x = files,pattern = "^[a-z]{2}[//_]{1}experimental_leiden.txt$")]

genesdf <- read.table(files[grepl(x = files,pattern = "^[a-z]{2}[//_]{1}experimental_leiden.txt$")],sep="\t", header = TRUE)
genesdf
genesdf <-genesdf[,-1]
row.names(genesdf) <- NULL
genesdf




#Create the list with the dataframe
genesList <- list()
for (i in seq_along(genesdf)){
  #print(genesdf[[i]])
  genesList[[i]] <- genesdf[[i]] 
}


#Set the names of the clusters...
names(genesList) <- c(rep(paste0("cluster_",seq(1:length(genesList)))))
names(genesList)
genesList
#genes <- data.frame(matrix(unlist(genesList),nrow=length(genesList),byrow = FALSE))



#Load danio rerio, data base
DrOrgDb <- org.Dr.eg.db



#Function for transforming the symbols to entrezgeneids
clusters_annotation  <- function(cluster){
  
  #Load ensembl, using mart and danio rerio.
  ensembl=biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
  #biomaRt::listAttributes(ensembl),list the attributes interested in.
  return(data.frame(biomaRt::getBM(attributes = c("entrezgene_id"),mart = ensembl,filters = "external_gene_name",values=cluster)))
}


########################################
#Annotate the lists with entrezgeneids #
########################################
annotated_clusters <- lapply(genesList, function(x) clusters_annotation(x))

#Check if it worked yet..

#lapply(annotated_clusters,function(x) print(x$entrezgene_id))

#################################
#           bp                  #
#################################
annotated_cluster_bp <- lapply(annotated_clusters, function (x) {
  enrichGO(gene = x$entrezgene_id ,OrgDb  =DrOrgDb, ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)
})



#################################
#           mf                  #
#################################

annotated_cluster_mf <- lapply(annotated_clusters, function (x) {
  enrichGO(gene = x$entrezgene_id ,OrgDb  =DrOrgDb, ont = "MF", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)
})


#################################
#           cc                  #
#################################
annotated_cluster_cc <- lapply(annotated_clusters, function (x) {
  enrichGO(gene = x$entrezgene_id ,OrgDb  =DrOrgDb, ont = "CC", pAdjustMethod = "BH",pvalueCutoff  = 0.05,readable = T)
})




#Jsut chech the names and selection is done properly with [x], [[x]] =)
lapply(names(annotated_cluster_bp), function(x) print(annotated_cluster_bp[[x]]))


################################################################
#                                                              #
# Save the data, as tables for opening later with office/excel #
#                                                              #
################################################################



#Write the tables of the data bp.
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/")
dir.create(recursive = TRUE,"./experimental_leiden/bp")
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/bp")
lapply(names(annotated_cluster_bp), function(x) write.table(annotated_cluster_bp[[x]], file = paste(x,".txt"), sep = "\t"))



#Write the tables of the data mf.
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/")
dir.create("./mf")
setwd("./mf")
lapply(names(annotated_cluster_mf), function(x) write.table(annotated_cluster_mf[[x]], file = paste(x,".txt"), sep = "\t"))


#Write the tables of the data cc.
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/")
dir.create("./cc")
setwd("./cc")
lapply(names(annotated_cluster_cc), function(x) write.table(annotated_cluster_cc[[x]], file = paste(x,".txt"), sep = "\t"))


#lapply(names(annotated_cluster_cc), function(x) print(annotated_cluster_cc[[x]]))


######################
#    Barplot save    #
######################

#######
# BP  #
#######
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/bp")

for ( i in names(annotated_cluster_bp)) {
  print(i)
  pdf_names=paste0("./control_",i,"_bp_barplot.pdf")
  pdf(pdf_names,width = 15,height = 15)
  print(barplot(annotated_cluster_bp[[i]],drop=T,showCategory = 25,title = paste0("control_",i,"_bp")))
  dev.off()
}

#lapply(names(annotated_cluster_mf), function(x){
#  pdf_names=paste0("./control_",x,"_bp_barplot.pdf")
#  pdf(file = pdf_names)
#  #print(annotated_cluster_mf[[x]])
#    print(barplot(annotated_cluster_bp[[i]],drop=T,showCategory = 25,title = paste0("control_",i,"_bp")))
#  dev.off()
#})
#######
# MF  #
#######
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/mf")
for ( i in names(annotated_cluster_mf)) {
  print(i)
  pdf_names=paste0("./control_",i,"_mf_barplot.pdf")
  pdf(pdf_names,width = 15,height = 15)
  print(barplot(annotated_cluster_mf[[i]],drop=T,showCategory = 25,title = paste0("control_",i,"_mf")))
  dev.off()
}

#######
# CC  #
#######
setwd("/mnt/haus/marius/SingleCell10xKallistoJupyter/R/Enrichment/experimental_leiden/cc")

lapply(names(annotated_cluster_cc), function(x){
  pdf_names=paste0("./control_",x,"_cc_barplot.pdf")
  pdf(file = pdf_names)
  #print(annotated_cluster_mf[[x]])
  print(barplot(annotated_cluster_cc[[i]],drop=T,showCategory = 25,title = paste0("control_",i,"_cc")))
  dev.off()
})


print("Done")

