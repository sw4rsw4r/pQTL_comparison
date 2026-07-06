wget https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/README

for i in `seq 1 22`
do
 wget https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr${i}_ld.tsv.gz
 wget https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr${i}_ld.tsv.gz.tbi
done
