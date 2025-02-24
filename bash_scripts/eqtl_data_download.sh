while read QTS QTD; do
    URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/${QTS}/${QTD}/${QTD}.lbf_variable.txt.gz"
    wget "$URL"
done < list_for_download.txt


while read QTS QTD; do
    URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/${QTS}/${QTD}/${QTD}.cc.tsv.gz"
    wget "$URL"
done < list_for_download.txt