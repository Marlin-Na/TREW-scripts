# TODO Note that there is no 'chrM' in the three gtf files

# egrep '(chr[0-9][0-9]?|chrX|chrY)[[:space:]]' hg19.gtf > mod.hg19.gtf
# egrep '(chr[0-9][0-9]?|chrX|chrY)[[:space:]]' mm10.gtf > mod.mm10.gtf
# egrep '(chrX|chrY|chr2[LR]|chr3[LR]|chr4)[[:space:]]' dm6.gtf > mod.dm6.gtf



function forsake-nonstandard-genomes
    egrep '((chr[0-9][0-9]?|chrX|chrY)|(chrX|chrY|chr2[LR]|chr3[LR]|chr4))[[:space:]]' $argv
end


