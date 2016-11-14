#!/usr/bin/fish
 
# Usage: split-by-chr [Source GTF file] [Target directory]
function split-by-chr
    if [ (count $argv) != 2 ]
        echo "Usage: split-by-chr [Source GTF file] [Target directory]"
    end

    set file $argv[1]
    set dir $argv[2]

    set chrs (awk '{print $1}' $file | sort -u)

    for chr in $chrs
        grep "^$chr" $file > $dir/$chr.$file
    end
end
