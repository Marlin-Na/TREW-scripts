#!/usr/bin/fish

# Remove all tracks except gene model and DNA

set jbrowsedir ~/jbrowse

cd $jbrowsedir

for genome in hg19 mm10 dm6
    set tracks (ls data/$genome/tracks -I gene_model)
    for track in $tracks
        bin/remove-track.pl -D --trackLabel $track --dir data/$genome
    end
end

# Sync to the server
rsync -r -v --delete ./ /var/www/html/jbrowse/
