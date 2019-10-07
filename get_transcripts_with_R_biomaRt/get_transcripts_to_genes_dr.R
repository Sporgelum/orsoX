library(BiocManager)
library(biomaRt)


#List marts availables
listMarts()
#listEnsembl()
#listDatasets()


#Conect to ensembl and use the specific gene ensembl wanted.
ensembl=biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
listDatasets(ensembl)[listDatasets(ensembl)[,1]=="drerio_gene_ensembl",]


#List what we want to get 
listAttributes(ensembl)

#get the data into a df.
transcripts_to_genes=data.frame(biomaRt::getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id","external_gene_name"),mart = ensembl))
#?useMart
dim(transcripts_to_genes);head(transcripts_to_genes)
transcripts_to_genes[transcripts_to_genes$ensembl_transcript_id_version=="ENSDART00000147162.2",]
write.table(transcripts_to_genes,"/mnt/haus/marius/SingleCell10xKallistoJupyter/transcripts_to_genes.txt",row.names = F,col.names = F,sep = "\t",quote = F)
