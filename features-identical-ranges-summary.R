

library(dplyr)
library(GenomicRanges)

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











overlaps <- function (df) {
    gr <- GRanges(df)
    hits <- findOverlaps(gr,gr)
    hits <- hits[queryHits(hits) != subjectHits(hits)]
    overlapped_gr <- gr[queryHits(hits)]
    overlapped_gr
}

grl <- needed %>%
    split(needed$Experiment_ID) %>%
    lapply(FUN = overlaps) %>%
    lapply(FUN = function(x) {
        mcols(x) <- NULL
        x
    })

grl

library(purrr)
    
grl %>%
    Filter(f = function(x) length(x) >= 1)
    
    


any_identical_rg <- function(gr) {
    gr[countOverlaps(query = gr, subject = gr, type = "equal") >= 2]
}

grl2 <- grl %>%
    Filter(f = function(x) length(x) >= 1) %>%
    Map(f = any_identical_rg) %>%
    Filter(f = function(x) length(x) >= 1)

grl2 %>% map(as.data.frame) %>%
    lapply(FUN = function(df) {
        id <- mapply(function(a, b) digest::sha1(c(a,b)), a = df$start, b = df$end)
        df$overlay_id <- id
        df
    })


library(trewjb)

ddf <- needed %>%
    group_by(Species, Experiment_ID, Chromosome, Start, End) %>%
    dplyr::filter(n() >= 2) %>%
    summarize(num_identical_fea = n(), gene_annotated = list(unique(Gene_ID))) %>%
    dplyr::arrange(desc(num_identical_fea)) %>%
    mutate(genome = ifelse(Species == "Mus musculus",
                            yes = "mm10", no = "hg19")) %>%
    dplyr::mutate(jbrowse_link = 
        Vectorize(getJbrowseLink)(
            genome = genome,
            chromosome = Chromosome,
            range_start = Start,
            range_end = End,
            hight_start = Start,
            hight_end = End,
            tracks = Experiment_ID,
            show_navagation = T,
            show_tracklist = T,
            show_overview = T
        )
    )

ddf$jbrowse_link <- {
    sapply(X = ddf$jbrowse_link, function(link) {
        htmltools::a(href = link, "Link") %>% as.character()
    })
}

DT::datatable(ddf,escape = F)
dt <- DT::datatable(ddf,escape = F)
#DT::saveWidget(dt, "~/summary_identical_feature_ranges.html",selfcontained = TRUE)

htmlTable::htmlTable(ddf)




