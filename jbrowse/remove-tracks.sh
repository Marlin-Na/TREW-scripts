#!/bin/bash

# Remove all current tracks

jbrowsedir=~/jbrowse

cd $jbrowsedir

## hg19
bin/remove-track.pl -D --trackLabel 'all_m6A' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m1A' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m5C' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_Psi' --dir data/hg19 
## mm10
bin/remove-track.pl -D --trackLabel 'all_m6A' --dir data/mm10 
bin/remove-track.pl -D --trackLabel 'all_m5C' --dir data/mm10 
## dm6
bin/remove-track.pl -D --trackLabel 'all_m6A' --dir data/dm6 


## hg19
bin/remove-track.pl -D --trackLabel 'all_m6A_less_col' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m1A_less_col' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_m5C_less_col' --dir data/hg19 
bin/remove-track.pl -D --trackLabel 'all_Psi_less_col' --dir data/hg19 
## mm10
bin/remove-track.pl -D --trackLabel 'all_m6A_less_col' --dir data/mm10 
bin/remove-track.pl -D --trackLabel 'all_m5C_less_col' --dir data/mm10 
## dm6
bin/remove-track.pl -D --trackLabel 'all_m6A_less_col' --dir data/dm6 



# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
