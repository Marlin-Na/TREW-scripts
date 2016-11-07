
function extract-multi-chr-genes 
    awk '{print $10 " " $1}' $argv| sort -u| awk '{print $1}'| awk 'a[$0]++'| sed 's/;//g'| sed 's/"//g'
end
