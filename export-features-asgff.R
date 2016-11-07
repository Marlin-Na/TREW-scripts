# #-- The following code copied from TREW shiny app ---------------------------------------
#
# library(sqldf)
# library(RMySQL)
#
# #### Load the database.
# #db <- dbConnect(MySQL(), dbname="TREW", username = "root", password = "yumenwh3920851")
# db <- dbConnect(SQLite(), dbname= "TREW.db")
# Genome_Location <- dbReadTable(db, "Genome_Location")
# Sites_Info <- dbReadTable(db, "Sites_Info")
# Source_Info <- dbReadTable(db, "Source_Info")
# #Gene_C <- unique(tolower(Sites_Info$Gene_ID))
# #### Generate the indexes.
# idx_2to1 = table(Genome_Location$Meth_Site_ID)
# idx_3to2 = table(Sites_Info$Source_ID)[Source_Info$DataSet_ID]
# idx_3to1 = table(rep(Sites_Info$Source_ID,idx_2to1))[Source_Info$DataSet_ID]
#
#
# Table1 = Genome_Location
# Table2 = Sites_Info
# Table3 = Source_Info
#
# Tb1 <- function(Species = NULL,
#                 Target = "ALL",
#                 Modification = NULL,
#                 Gene_ID = ".",
#                 Include_Liftover = TRUE)
# {
#   # Initialize indexes
#   idx2 = !vector(length = dim(Table2)[1])
#   idx3 = !vector(length = dim(Table3)[1])
#   # Select species
#   if (!is.null(Species)) {
#     idx3 <- Table3$Species == Species
#   }else{}
#   # Select regulators
#   if (Target == "Regulator") {
#     idx3 <- idx3 & !Table3$Target %in% c("YTHDF1","YTHDF2","YTHDC1")
#   }else{
#     if (Target == "Regulatee") {
#       idx3 <- idx3 & Table3$Target %in% c("YTHDF1","YTHDF2","YTHDC1")
#     }else{}
#   }
#   # Select Modifications.
#   if (!is.null(Modification)) {
#     idx3 <- idx3 & (Table3$Modification == Modification)
#   }else{}
#
#   # Select LiftOvers.
#   if (!Include_Liftover) {
#     idx3 <- idx3 & is.na(Table3$LiftOver)
#   }else{}
#
#   # Select Genes.
#      # hit_idx <- grepl(Gene_ID, Table2$Gene_ID, ignore.case = TRUE)
#   hit_idx <- rep(T,nrow(Table2)) # <------------------------------- Fake the gene id to select ALL!
#
#   idx2 <-  hit_idx & rep(idx3,idx_3to2)
#
#   # Merge the table.
#   idx2[which(is.na(idx2))] <- FALSE
#
#   idx1 <- rep(idx2,idx_2to1)
#
#   Tb1 <- cbind(Table1[which(idx1),],
#                Table2[rep(which(idx2),
#                           idx_2to1[which(idx2)]),],
#                Table3[rep(1:dim(Table3)[1],
#                           idx_3to1)[which(idx1)],])
#   Tb1
# }
#
# #-- Copy End  ----------------------------------------------------------------------------------


library(shiny)
library(DT)
library(readr)
Table1 <- read_tsv("Table1.txt")
Table2 <- read_tsv("Table2.txt")
Table3 <- read_tsv("Table3.txt")
idx_2to1 <- inverse.rle(read_rds("idx_2to1.rle.rds"))
idx_3to2 <- read_rds("idx_3to2.rds")
idx_3to1 <- read_rds("idx_3to1.rds")

Tb1 <- function(Species = "All",
                Target = "All",
                Modification = "All",
                Gene_ID = ".",
                Include_Liftover = "Yes")
{
  # Initialize indexes
  idx2 = !vector(length = dim(Table2)[1])
  idx3 = !vector(length = dim(Table3)[1])
  # Select species
  if (Species != "All") {
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
  if (Modification != "All") {
    idx3 <- idx3 & (Table3$Modification == Modification)
  }else{}

  # Select LiftOvers.
  if (Include_Liftover != "Yes") {
    idx3 <- idx3 & is.na(Table3$LiftOver)
  }else{}

  # Select Genes.
#### Modified here to select all genes -------------------------------------------------------------------------------------------------
# hit_idx <- grepl(Gene_ID, Table2$Gene_ID, ignore.case = TRUE)
  hit_idx <- rep(T, length(Table2$Gene_ID))


  if(sum(hit_idx) == 0){stop("Not found your gene, please input Gene Symbol....")}
  idx2 <-  hit_idx & rep(idx3,idx_3to2)

  # Merge the table.
  idx2[which(is.na(idx2))] <- FALSE

  idx1 <- rep(idx2,idx_2to1)

  Tb1 <- cbind(Table1[which(idx1),],
               Table2[rep(which(idx2),
                          idx_2to1[which(idx2)]),],
               Table3[rep(1:dim(Table3)[1],
                          idx_3to1)[which(idx1)],])
  cat("Tb1 run once\n")
  Tb1
}


count.x <- function(tb1.x,tb2.y){
  tb2.y <- tb2.y[,c(1,2)]
  tb <- table(tb1.x$Target,tb1.x$Gene_ID)
  idx <- which(tb2.y[,1]%in%rownames(tb) & tb2.y[,2]%in%colnames(tb))
  hits <- vector("numeric",nrow(tb2.y))
  hits[idx] <- tb[cbind(as.character(tb2.y[idx,1]),as.character(tb2.y[idx,2]))]
  hits
}

Tb2 <- function(Tb1,Test = FALSE){
  Tb2 <- unique(data.frame(Target = Tb1$Target,
                           Gene_ID = Tb1$Gene_ID,
                           Target_type = Tb1$Target_type,
                           Modification = Tb1$Modification))
  Tb2$Record_num <- count.x(Tb1,Tb2)
  Tb2$Consistent_num <- count.x(Tb1[which(Tb1$Consistency > 0),],Tb2)
  Tb2$Positive_num <- count.x(Tb1[which(Tb1$Diff_p_value < .05),],Tb2)
  Tb2$Positive_percent <- paste(round(100*(Tb2$Positive_num/Tb2$Record_num),2),"%",sep = "")
  if(Test){
    binom_list <- apply(Tb2[,c("Positive_num","Record_num")],1,function(x) binom.test(x[1],x[2],p = 0.2654169))
    Tb2$binom_p <- sapply(binom_list,function(x) x$p.value)
    Tb2$binom_CI <- sapply(binom_list,function(x) paste("[",round(x$conf.int[1],3),",",round(x$conf.int[2],3),"]", sep = ""))
  }else{}
  cat("Tb2 run once\n")
  Tb2
}

###Table 3 is the specific table that must be inferred from tb1 and tb2

Tb3 <- function(Tb1,Tb2,Select_Number = 1:dim(Tb2)[1],Return_All = "No")
{
  Tb2_s <- Tb2[Select_Number,]

  Tb3 <- Tb1[which(Tb1$Target %in% Tb2_s$Target &
                     Tb1$Modification %in% Tb2_s$Modification &
                     Tb1$Gene_ID %in% Tb2_s$Gene_ID),]

if (Return_All != "Yes") {
    Tb3 <- Tb3[,c(32,2,31,41,39,34,33,14,12,13)]
  } else {
    Tb3 <- Tb3[,c(32,2,31,12,1,6,5,3,4,13,41,33,34,35,36,39,40,42,9,10,15,16,11,14,17,18,19,20,21,22,23,24,25,26,27,28,38,43,37)]
  }
  cat("Tb3 run once\n")
  Tb3
}

Tb3_DT <- function(Tb3)
{
  DT::datatable(Tb3, rownames= FALSE, filter = "top",
                extensions = list("ColReorder" = NULL,
                                  "Buttons" = NULL,
                                  "FixedHeader" = NULL,
                                  "Scroller" = NULL),
                options = list(
                  scrollX = TRUE,
                  deferRender = TRUE,
                  scrollY = 400,
                  scroller = TRUE,
                  dom = 'BRrlftpi',
                  autoWidth=TRUE,
                  fixedHeader = TRUE,
                  lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                  ColReorder = TRUE,
                  buttons =
                    list(
                      'copy',
                      'print',
                      list(
                        extend = 'collection',
                        buttons = c('csv', 'excel'),
                        text = 'Download'
                      ),
                      I('colvis')
                    )
                ))
}

Tb2_DT <- function(Tb2){
  DT::datatable(Tb2,
                rownames = FALSE,
                options = list(
                  scrollX = TRUE,
                  deferRender = TRUE,
                  scrollY = 340,
                  dom = 'lrtip',
                  autoWidth = TRUE,
                  ColReorder = TRUE)
  )
}


#####---------------------------------------------------------------










library(dplyr)
library(rtracklayer)
library(GenomicRanges)

all_full_table <- Tb1()
all_full_table <- tbl_df(all_full_table)

all_full_table <- all_full_table %>%
    select(
        Chromosome,
        Range_start,
        Range_width,
        Strand,

        Meth_Range_ID,
        Meth_Site_ID,
        Diff_p_value,
        Diff_fdr,
        Diff_log2FoldChange,
        Gene_ID,
        Log2_RPKM_Wt,
        Log2_RPKM_Treated,
        Genome_assembly,
        Modification,
        Technique,
        Target,
        Target_type,
        Perturbation,
        Cell_line,
        Treatment,
        Species,
        LiftOver,

        DataSet_ID
    )

# Extract with which columns
#
# |   | Column Name             |                      Sample Value |
# |---+-------------------------+-----------------------------------|
# | T | Meth_Range_ID           |                              8333 |
# | T | Meth_Site_ID            |                              6101 |
# |   | Methylation_ID          |                              6101 |
# | T | Diff_p_value            |                 0.549540873857625 |
# | T | Diff_fdr                |                 0.621260615022094 |
# | T | Diff_log2FoldChange     |                             0.218 |
# | T | Gene_ID                 |                            HS6ST1 |
# |   | Source_ID               |         m6A_Ftokd_MEFHS_lift.hg19 |
# |   | Consistency             |                                 1 |
# | T | Log2_RPKM_Wt            |                 0.729256748836259 |
# | T | Log2_RPKM_Treated       |                  3.33121584334439 |
# | X | Overlap_UTR5            |                                 0 |
# | X | Overlap_CDS             |                                 1 |
# | X | Overlap_UTR3            |                                 0 |
# | X | Overlap_mRNA            |                                 1 |
# | X | Overlap_lncRNA          |                                 0 |
# | X | Overlap_sncRNA          |                                 0 |
# | X | Overlap_tRNA            |                                 0 |
# | X | Overlap_miRNA           |                                 0 |
# | X | Overlap_miRNATS         |                                 0 |
# | X | Distance_ConsensusMotif |                                 0 |
# | X | Distance_StartCodon     |                             50133 |
# | X | Distance_StopCodon      |                                52 |
# | ? | DataSet_ID              |         m6A_Ftokd_MEFHS_lift.hg19 |
# | T | Genome_assembly         |                              hg19 |
# | T | Modification            |                               m6A |
# | T | Technique               |                             MeRIP |
# | T | Target                  |                               Fto |
# | T | Target_type             |                            eraser |
# | T | Perturbation            |                                KD |
# |   | Date_of_process         |                          5-Jul-16 |
# |   | Paper                   | Dynamic m(6)A mRNA methylation... |
# | T | Cell_line               |                               MEF |
# | T | Treatment               |                        Heat Shock |
# | T | Species                 |                      Mus musculus |
# | T | LiftOver                |                 Lift mm10 to hg19 |
# |   | Computation_pepline     |      trim_galore,tophat,exomepeak |


as_gr <- function(table) {
    dplyr::rename(table, seqnames = Chromosome,start = Range_start, width = Range_width, strand = Strand) %>%
    mutate(end = start + width - 1) %>% as("GRanges")
}


## By Modification and export to gff

for (modi in unique(all_full_table$Modification)) {
    modi_df <- filter(all_full_table, Modification == modi)

    for (geno in unique(modi_df$Genome_assembly)) {

        modi_geno_gr <- filter(modi_df, Genome_assembly == geno) %>% as_gr

        folder <- 'TREW-scripts/gff-features/by_modification/'
        path <- paste0(folder,geno,'_',modi,'.','gff3')

        export.gff3(modi_geno_gr, path)
    }
}


# m6A_full_table <- filter(all_full_table, Modification == 'm6A')
# m1A_full_table <- filter(all_full_table, Modification == 'm1A')
# m5C_full_table <- filter(all_full_table, Modification == 'm5C')
# Psi_full_table <- filter(all_full_table, Modification == 'Psi')
#
#
#
# dm6_m6A_full_gr <- filter(m6A_full_table, Genome_assembly == 'dm6') %>% as_gr
# dm6_m1A_full_gr <- filter(m1A_full_table, Genome_assembly == 'dm6') %>% as_gr ## Empty!!!
# dm6_m5C_full_gr <- filter(m5C_full_table, Genome_assembly == 'dm6') %>% as_gr ## Empty!!!
# dm6_Psi_full_gr <- filter(Psi_full_table, Genome_assembly == 'dm6') %>% as_gr ## Empty!!!
#
#
# mm10_m6A_full_gr <- filter(m6A_full_table, Genome_assembly == 'mm10') %>% as_gr
# mm10_m1A_full_gr <- filter(m1A_full_table, Genome_assembly == 'mm10') %>% as_gr ## Empty!!!
# mm10_m5C_full_gr <- filter(m5C_full_table, Genome_assembly == 'mm10') %>% as_gr
# mm10_Psi_full_gr <- filter(Psi_full_table, Genome_assembly == 'mm10') %>% as_gr ## Empty!!!
#
#
# hg19_m6A_full_gr <- filter(m6A_full_table, Genome_assembly == 'hg19') %>% as_gr
# hg19_m1A_full_gr <- filter(m1A_full_table, Genome_assembly == 'hg19') %>% as_gr
# hg19_m5C_full_gr <- filter(m5C_full_table, Genome_assembly == 'hg19') %>% as_gr
# hg19_Psi_full_gr <- filter(Psi_full_table, Genome_assembly == 'hg19') %>% as_gr
#
#
# ## To GFF
#
# export.gff3(hg19_m6A_full_gr,'../gff_features/hg19_m6A.gff3')
# export.gff3(hg19_m1A_full_gr,'../gff_features/hg19_m1A.gff3')
# export.gff3(hg19_m5C_full_gr,'../gff_features/hg19_m5C.gff3')
# export.gff3(hg19_Psi_full_gr,'../gff_features/hg19_Psi.gff3')
#
# export.gff3(mm10_m6A_full_gr,'../gff_features/mm10_m6A.gff3')
# export.gff3(mm10_m1A_full_gr,'../gff_features/mm10_m1A.gff3') # Failed
# export.gff3(mm10_m5C_full_gr,'../gff_features/mm10_m5C.gff3')
# export.gff3(mm10_Psi_full_gr,'../gff_features/mm10_Psi.gff3') # Failed
#
# export.gff3(dm6_m6A_full_gr,'../gff_features/dm6_m6A.gff3')
# export.gff3(dm6_m1A_full_gr,'../gff_features/dm6_m1A.gff3') # Failed
# export.gff3(dm6_m5C_full_gr,'../gff_features/dm6_m5C.gff3') # Failed
# export.gff3(dm6_Psi_full_gr,'../gff_features/dm6_Psi.gff3') # Failed



####### ------------------------------------------------------------

datasets <- group_by(all_full_table, DataSet_ID)

modi2datasets <- distinct(datasets, Modification)

## Examine that one dataset will have only one type of modification
distinct(count(modi2datasets),n)




li.datasets <- split(all_full_table, all_full_table$DataSet_ID) %>%
    lapply(as_gr)

li.hg19 <-
    filter(all_full_table, Genome_assembly == 'hg19') %>%
    split(filter(all_full_table, Genome_assembly == 'hg19')$DataSet_ID) %>%
    lapply(as_gr)

li.mm10 <-
    filter(all_full_table, Genome_assembly == 'mm10') %>%
    split(filter(all_full_table, Genome_assembly == 'mm10')$DataSet_ID) %>%
    lapply(as_gr)

li.dm6 <-
    filter(all_full_table, Genome_assembly == 'dm6') %>%
    split(filter(all_full_table, Genome_assembly == 'dm6')$DataSet_ID) %>%
    lapply(as_gr)



## Export each dataset to gff

for (datasetname in names(li.hg19)) {
    gr.set <- li.hg19[[datasetname]]
    export.gff3(gr.set, paste0('TREW-scripts/gff-features/hg19/',datasetname))
}

for (datasetname in names(li.mm10)) {
    gr.set <- li.mm10[[datasetname]]
    export.gff3(gr.set, paste0('TREW-scripts/gff-features/mm10/',datasetname))
}

for (datasetname in names(li.dm6)) {
    gr.set <- li.dm6[[datasetname]]
    export.gff3(gr.set, paste0('TREW-scripts/gff-features/dm6/',datasetname))
}





