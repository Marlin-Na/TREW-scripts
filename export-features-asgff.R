
library(dplyr)
library(GenomicRanges)


all.files <- list.files("gff-features/",full.names = T,recursive = T)
for (i in all.files) file.remove(i)

#setwd("~/msync/projects/TREW/TREW-scripts/")
db <- src_sqlite(path = './TREW_0.2.0.db')
db

t1 <- tbl(db, 'Genome_Location')
t2 <- tbl(db, 'Raw_data_records')
t3 <- tbl(db, 'Sites_Info')
t4 <- tbl(db, 'Source_Info')

joined <- t1 %>%
    left_join(t3) %>%
    left_join(t4) %>%
    # Do not join `Raw_data_records`:
    #     One Experiment_ID can correspond to multiple GEO_ID
    # left_join(t2) %>%
    collect(n = Inf)
joined


colnames(joined)

needed <- joined %>%
    # Select the needed meta columns
    dplyr::select(
        Range_ID,
        Site_ID,
        Start,
        End,
        Strand,
        Chromosome,
      # Note_t1,
        Pvalue,
      # Padj,
      # P_treated,
      # P_control,
      # Counts_WT_IP,
      # Counts_WT_INPUT,
      # Counts_TR_IP,
      # Counts_TR_INPUT,
      # Log2_RR,
      # Log2_OR,
      # Q,
        Experiment_ID,
      # Note_t2,
        Gene_ID,
      # Overlap_UTR5,
      # Overlap_CDS,
      # Overlap_UTR3,
      # Overlap_mRNA,
      # Overlap_lncRNA,
      # Overlap_sncRNA,
      # Overlap_tRNA,
      # Overlap_miRNA,
      # Overlap_miRNATS,
      # Distance_ConsensusMotif,
      # Distance_TSS,
      # Distance_StopCodon,
        Genome_assembly,
        Modification,
        Technique,
        Target,
        Target_type,
        Perturbation,
      # Date_of_process,
        Paper,
        Cell_line,
        Treatment,
        Species
      # Computation_pepline,
      # Note_t3,
      
      
      # These columns are in 'Raw_data_records'
        # GEO_ID,
        # IP_Input,
        # Genotype,
        # Replicate
    )

### TEMP:  Remove three tracks    ------------------------------

needed %>% group_by(Genome_assembly, Experiment_ID) %>%
    summarize(num = n()) %>% summarize(num = n())

needed <- needed %>%
    filter(
        ((!(Experiment_ID %in% c("m6A-YTHDF2-kd-MEF", "m6A-YTHDF2-kd-MEF-HS"))) &
            Genome_assembly == "mm10"
        ) |
        (Experiment_ID != "m6A-YTHDF2-oe-Hela" & Genome_assembly == "hg19") |
        (Genome_assembly == "dm6")
    )

needed %>% group_by(Genome_assembly, Experiment_ID) %>%
    summarize(num = n()) %>% summarize(num = n())

### TEMP END                    ------------------------------


exportgff <- function(df,
                      path = file.path('gff-features',
                                       unique(df$Genome_assembly),
                                       unique(df$Experiment_ID))) {
    gr <- as(df, "GRanges")
    rtracklayer::export.gff3(gr, path)

    return('Done!')
}

library(purrr)

### Export by genome and experiment id
needed %>%
    split(.$Genome_assembly) %>%
    map(~ split(.x, .x$Experiment_ID)) %>%
    at_depth(2, exportgff)




exportbymod <- function(df) {
    path <- file.path('gff-features',
                      'by_modification',
                      paste(sep = '_',
                            df$Genome_assembly,
                            df$Modification) %>% unique)
    exportgff(df, path)
}


### Export by genome and modification
needed %>%
    split(paste(sep = '___', .$Genome_assembly, .$Modification)) %>%
    map(exportbymod)








