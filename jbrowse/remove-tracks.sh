#!/bin/bash

# Remove gff features with the name index

jbrowsedir=~/jbrowse

cd $jbrowsedir

## hg19
bin/remove-track.pl -D --trackLabel 'all_m6A_with_name' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m1A_with_name' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m5C_with_name' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_Psi_with_name' --dir data/hg19 
## mm10
bin/remove-track.pl -D --trackLabel 'all_m6A_with_name' --dir data/mm10 
bin/remove-track.pl -D --trackLabel 'all_m5C_with_name' --dir data/mm10 
## dm6
bin/remove-track.pl -D --trackLabel 'all_m6A_with_name' --dir data/dm6 


# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
