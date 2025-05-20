#cd MFPD/
mkdir -p db
makeblastdb -in ./cluster/clades_95.fa -dbtype nucl -out ./db/clades_blastn_db
mmseqs createdb ./cluster/clades_95.fa ./db/clades_mmseq_db
mmseqs createindex ./db/clades_mmseq_db temp --search-type 3
