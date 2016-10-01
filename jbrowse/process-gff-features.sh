#!/bin/bash

# Process the features with less columns

jbrowsedir=~/jbrowse
gffdir=~/lesscol_gff_features

cd $jbrowsedir

# Process gff files
## hg19
bin/flatfile-to-json.pl --gff $gffdir/hg19_m6A.gff3 --trackLabel 'all_m6A_less_col' --out data/hg19 --trackType CanvasFeatures --key 'm6A (less col)' --sortMem 50000000000
bin/flatfile-to-json.pl --gff $gffdir/hg19_m1A.gff3 --trackLabel 'all_m1A_less_col' --out data/hg19 --trackType CanvasFeatures --key 'm1A (less col)' --sortMem 50000000000
bin/flatfile-to-json.pl --gff $gffdir/hg19_m5C.gff3 --trackLabel 'all_m5C_less_col' --out data/hg19 --trackType CanvasFeatures --key 'm5C (less col)' --sortMem 50000000000
bin/flatfile-to-json.pl --gff $gffdir/hg19_Psi.gff3 --trackLabel 'all_Psi_less_col' --out data/hg19 --trackType CanvasFeatures --key 'Psi (less col)' --sortMem 50000000000
## mm10
bin/flatfile-to-json.pl --gff $gffdir/mm10_m6A.gff3 --trackLabel 'all_m6A_less_col' --out data/mm10 --trackType CanvasFeatures --key 'm6A (less col)' --sortMem 50000000000
bin/flatfile-to-json.pl --gff $gffdir/mm10_m5C.gff3 --trackLabel 'all_m5C_less_col' --out data/mm10 --trackType CanvasFeatures --key 'm5C (less col)' --sortMem 50000000000
## dm6
bin/flatfile-to-json.pl --gff $gffdir/dm6_m6A.gff3 --trackLabel 'all_m6A_less_col' --out data/dm6 --trackType CanvasFeatures --key 'm6A (less col)' --sortMem 50000000000


# Generate name index
bin/generate-names.pl --out data/hg19 --mem 50000000000 &&
bin/generate-names.pl --out data/mm10 --mem 50000000000 &&
bin/generate-names.pl --out data/dm6 --mem 50000000000



# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
