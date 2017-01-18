

# Extract the data frame of the featues -----------------------


library(dplyr)
library(GenomicRanges)

db <- src_sqlite(path = './TREW_0.2.0.db')
db

t1 <- tbl(db, 'Genome_Location')
t2 <- tbl(db, 'Raw_data_records')
t3 <- tbl(db, 'Sites_Info')
t4 <- tbl(db, 'Source_Info')

joined <- t1 %>%
    left_join(t3) %>%
    left_join(t4) %>%
    left_join(t2) %>%
    collect(n = Inf)
joined


colnames(joined)

# --------------------------------------------------------------


# Check the ambiguity of genes on different strands

fea_genes <- joined %>%
    distinct(Genome_assembly, Chromosome, Strand, Gene_ID)

fea_genes %>%
    group_by(Gene_ID) %>%
    summarise(ambiguity = n()) %>%
    `$` (ambiguity) %>%
    table









