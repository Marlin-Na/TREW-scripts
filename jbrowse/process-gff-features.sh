#!/bin/bash

# Process the features with name index of methylation id

jbrowsedir=~/jbrowse
gffdir=~/package_gff_features

cd $jbrowsedir

# Process gff files
## hg19
bin/flatfile-to-json.pl --gff $gffdir/hg19_m6A.gff3 --trackLabel 'all_m6A_with_name' --out data/hg19 --trackType CanvasFeatures --key 'm6A (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
bin/flatfile-to-json.pl --gff $gffdir/hg19_m1A.gff3 --trackLabel 'all_m1A_with_name' --out data/hg19 --trackType CanvasFeatures --key 'm1A (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
bin/flatfile-to-json.pl --gff $gffdir/hg19_m5C.gff3 --trackLabel 'all_m5C_with_name' --out data/hg19 --trackType CanvasFeatures --key 'm5C (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
bin/flatfile-to-json.pl --gff $gffdir/hg19_Psi.gff3 --trackLabel 'all_Psi_with_name' --out data/hg19 --trackType CanvasFeatures --key 'Psi (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
## mm10
bin/flatfile-to-json.pl --gff $gffdir/mm10_m6A.gff3 --trackLabel 'all_m6A_with_name' --out data/mm10 --trackType CanvasFeatures --key 'm6A (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
bin/flatfile-to-json.pl --gff $gffdir/mm10_m5C.gff3 --trackLabel 'all_m5C_with_name' --out data/mm10 --trackType CanvasFeatures --key 'm5C (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'
## dm6
bin/flatfile-to-json.pl --gff $gffdir/dm6_m6A.gff3 --trackLabel 'all_m6A_with_name' --out data/dm6 --trackType CanvasFeatures --key 'm6A (with name)' --sortMem 50000000000 --nameAttributes 'methylation_id'


# Generate name index
bin/generate-names.pl --out data/hg19 --mem 50000000000 &&
bin/generate-names.pl --out data/mm10 --mem 50000000000 &&
bin/generate-names.pl --out data/dm6 --mem 50000000000



# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
