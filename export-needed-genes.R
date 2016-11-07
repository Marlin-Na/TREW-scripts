

# Extract the data frame of the featues -----------------------

library(dplyr)

library(sqldf)
library(RMySQL)

#### Load the database.
#db <- dbConnect(MySQL(), dbname="TREW", username = "root", password = "yumenwh3920851")
db <- dbConnect(SQLite(), dbname= "TREW.db")
Genome_Location <- dbReadTable(db, "Genome_Location")
Sites_Info <- dbReadTable(db, "Sites_Info")
Source_Info <- dbReadTable(db, "Source_Info")
#Gene_C <- unique(tolower(Sites_Info$Gene_ID))
#### Generate the indexes.
idx_2to1 = table(Genome_Location$Meth_Site_ID)
idx_3to2 = table(Sites_Info$Source_ID)[Source_Info$DataSet_ID]
idx_3to1 = table(rep(Sites_Info$Source_ID,idx_2to1))[Source_Info$DataSet_ID]


Table1 = Genome_Location
Table2 = Sites_Info
Table3 = Source_Info

Tb1 <- function(Species = NULL,
                Target = "ALL",
                Modification = NULL,
                Gene_ID = ".",
                Include_Liftover = TRUE)
{
  # Initialize indexes
  idx2 = !vector(length = dim(Table2)[1])
  idx3 = !vector(length = dim(Table3)[1])
  # Select species
  if (!is.null(Species)) {
    idx3 <- Table3$Species == Species
  }else{}
  # Select regulators
  if (Target == "Regulator") {
    idx3 <- idx3 & !Table3$Target %in% c("YTHDF1","YTHDF2","YTHDC1")
  }else{
    if (Target == "Regulatee") {
      idx3 <- idx3 & Table3$Target %in% c("YTHDF1","YTHDF2","YTHDC1")
    }else{}
  }
  # Select Modifications.
  if (!is.null(Modification)) {
    idx3 <- idx3 & (Table3$Modification == Modification)
  }else{}

  # Select LiftOvers.
  if (!Include_Liftover) {
    idx3 <- idx3 & is.na(Table3$LiftOver)
  }else{}

  # Select Genes.
     # hit_idx <- grepl(Gene_ID, Table2$Gene_ID, ignore.case = TRUE)
  hit_idx <- rep(T,nrow(Table2)) # <------------------------------- Fake the gene id to select ALL!

  idx2 <-  hit_idx & rep(idx3,idx_3to2)

  # Merge the table.
  idx2[which(is.na(idx2))] <- FALSE

  idx1 <- rep(idx2,idx_2to1)

  Tb1 <- cbind(Table1[which(idx1),],
               Table2[rep(which(idx2),
                          idx_2to1[which(idx2)]),],
               Table3[rep(1:dim(Table3)[1],
                          idx_3to1)[which(idx1)],])
  Tb1
}

features_table <- Tb1()
# features_table <- Tb1() %>% tbl_df

# END         -------------------------------------------------



## Extract Gene ID of Features
gene_feas <- features_table %>%
    distinct(Gene_ID) %>%
    arrange(Gene_ID)

gene_feas <- gene_feas$Gene_ID

## Import GTF as TxDb

library(GenomicFeatures)

# hg19.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/hg19.gtf')
# mm10.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/mm10.gtf')
# dm6.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/dm6.gtf')
ori.genes.hg19 <- genes(makeTxDbFromGFF(file = 'tmp/gtf-genemodel/hg19.gtf'))
ori.genes.mm10 <- genes(makeTxDbFromGFF(file = 'tmp/gtf-genemodel/mm10.gtf'))
ori.genes.dm6 <- genes(makeTxDbFromGFF(file = 'tmp/gtf-genemodel/dm6.gtf'))

hg19.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/mod.hg19.gtf')
mm10.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/mod.mm10.gtf')
dm6.txdb <- makeTxDbFromGFF(file = 'tmp/gtf-genemodel/mod.dm6.gtf')

genes.hg19 <- genes(hg19.txdb)
genes.mm10 <- genes(mm10.txdb)
genes.dm6 <-  genes(dm6.txdb)

# By modification, some of the genes can thus be imported,
# but it seems that there are also lost ones
# See the following:
library(VennDiagram)
# hg19
draw.pairwise.venn(
  length(ori.genes.hg19$gene_id),
  length(genes.hg19$gene_id),
  length(intersect(ori.genes.hg19$gene_id,genes.hg19$gene_id)),
  category = c('Original','Modified'))
# mm10
draw.pairwise.venn(
  length(ori.genes.mm10$gene_id),
  length(genes.mm10$gene_id),
  length(intersect(ori.genes.mm10$gene_id,genes.mm10$gene_id)),
  category = c('Original','Modified'))
# dm6
draw.pairwise.venn(
  length(ori.genes.dm6$gene_id),
  length(genes.dm6$gene_id),
  length(intersect(ori.genes.dm6$gene_id,genes.dm6$gene_id)),
  category = c('Original','Modified'))



# Test for intersects between different genomes
intersect(genes.hg19$gene_id, genes.mm10$gene_id)  %>% length
intersect(genes.hg19$gene_id, genes.dm6$gene_id)   %>% length
intersect(genes.mm10$gene_id, genes.dm6$gene_id)   %>% length


df.genes.hg19 <- as.data.frame(genes.hg19) %>% tbl_df %>% mutate(genome_assembly='hg19')
df.genes.mm10 <- as.data.frame(genes.mm10) %>% tbl_df %>% mutate(genome_assembly='mm10')
df.genes.dm6 <- as.data.frame(genes.dm6)   %>% tbl_df %>% mutate(genome_assembly='dm6')

df.genes.all <- rbind(df.genes.hg19,df.genes.mm10,df.genes.dm6)

## Save all the gene info as file
saveRDS(df.genes.all, file = 'all.dataframe_genes.Rds', ascii = F)









# ## Subset Genes of the features
#
# TODO !!!!! Also have to make sure that every gene in gene_feas have coorsbonding entry in df.genes !!!
#
# It turns out that some of the genes in gtf file can not be transferred into txdb

df.genes.match <- subset(df.genes.all, gene_id %in% gene_feas)
dim(df.genes.match)

df.genes.match.Ncase <- subset(df.genes.all, tolower(gene_id) %in% tolower(gene_feas))
dim(df.genes.match.Ncase)


notin.feas <- gene_feas[!(gene_feas %in% df.genes.match$gene_id)]
length(notin.feas)

notin.feas.Ncase <- gene_feas[!(gene_feas %in% df.genes.match.Ncase$gene_id)]
length(notin.feas.Ncase)

identical(notin.feas, notin.feas.Ncase)


notin.feas


# TODO There are still 15 genes located at multiple chromosomes

