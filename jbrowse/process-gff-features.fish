
# Process the features

## Setup dir variable
set jbrowsedir ~/jbrowse
set gffdir ~/gff-features

cd $jbrowsedir

## Function for common arguments

function processgff
    bin/flatfile-to-json.pl --trackType CanvasFeatures --sortMem 50000000000
end

## Process the by modification tracks
# hg19
processgff --gff $gffdir/by_modification/hg19_m6A.gff3 --trackLabel 'all_m6A' --out data/hg19  --key 'All m6A sites' 
processgff --gff $gffdir/by_modification/hg19_m1A.gff3 --trackLabel 'all_m1A' --out data/hg19  --key 'All m1A sites' 
processgff --gff $gffdir/by_modification/hg19_m5C.gff3 --trackLabel 'all_m5C' --out data/hg19  --key 'All m5C sites' 
processgff --gff $gffdir/by_modification/hg19_Psi.gff3 --trackLabel 'all_Psi' --out data/hg19  --key 'All Psi sites' 
# mm10
processgff --gff $gffdir/by_modification/mm10_m6A.gff3 --trackLabel 'all_m6A' --out data/mm10  --key 'All m6A sites' 
processgff --gff $gffdir/by_modification/mm10_m5C.gff3 --trackLabel 'all_m5C' --out data/mm10  --key 'All m5C sites' 
# dm6
processgff --gff $gffdir/by_modification/dm6_m6A.gff3 --trackLabel 'all_m6A' --out data/dm6  --key 'All m6A sites' 


## Process the by dataset tracks
for geno in hg19 mm10 dm6
    set feas (ls $gffdir/$geno)

    for fea in $feas
        processgff --gff $gffdir/$geno/$fea --trackLabel "$fea" --out data/$geno --key "$fea"
    end
end


## Generate name index
bin/generate-names.pl --out data/hg19 --mem 50000000000
bin/generate-names.pl --out data/mm10 --mem 50000000000
bin/generate-names.pl --out data/dm6 --mem 50000000000



# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
