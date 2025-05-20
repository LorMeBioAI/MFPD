import os
import pandas as pd
import argparse
import sys
import time
import subprocess

# Add parameters
parser = argparse.ArgumentParser(description='MFPD')
parser.add_argument('--file', required=True,help='Sample table')
parser.add_argument('--pwd', type=str, default='./', help='Output path, default is current path')
parser.add_argument('--similarity', type=float, default=0.97,
                    help='Recommended similarity: 0.99 for Full ITS AND 0.97 for ITS1/2 , default=0.97')
parser.add_argument('--region', type=str, required=True,
                    help='The sequencing region of your data, subregion or full-length')
args = parser.parse_args()

valid_regions = ['subregion', 'full-length']
if args.region not in valid_regions:
    print(f"[ERROR] Invalid region '{args.region}'. Must be one of: {', '.join(valid_regions)}")
    sys.exit(1)

def progress_bar():
    for i in range(1, 101):
        print("\r", end="")
        print(f"Running progress: {i}%: ", "▋" * (i // 2), end="")
        sys.stdout.flush()
        time.sleep(0.005)

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"--- Created new folder: {path} ---")
    else:
        print(f"--- Folder already exists: {path} ---")

# Obtain the software path, convert it to an absolute path
abs_path = os.path.split(os.path.realpath(__file__))[0]
abs_raw = os.path.realpath(args.file)
abs_db = os.path.join(abs_path, "db")
abs_bin = os.path.join(abs_path, "bin")
abs_cluser = os.path.join(abs_path, "clusters")
similarity = args.similarity

print(f"""
The PWD of workplace is: {abs_path}
The PWD of raw.fq.list is: {abs_raw}
""")

# Make output path
mkdir(args.pwd)
os.chdir(args.pwd)
mkdir("00.Data")
mkdir("01.ASV.tax")
mkdir("shell")

# Read input table
raw_fq_list = pd.read_csv(abs_raw, sep="\t", header=None)
ID_list = raw_fq_list[0]
fq_list = raw_fq_list[1]
num = ID_list.shape[0]

# Used paths
PWD = os.getcwd() + "/"
abs_shell = PWD + "shell/"
mode = PWD + "shell/S01.1.{}"

# Generate the first part script (S01.1.pathogen.part1.sh)
with open(mode.format("pathogen.part1") + '.sh', 'w', encoding='utf-8') as f1:
    content1 = f"""
cd {PWD}01.ASV.tax/
Rscript {abs_shell}S01.1.dada2.r
perl {abs_bin}/dada2normal.pl table.rechim.xls
perl {abs_bin}/asv-generate.pl table.rechim.xls asv.table asv.fa
sed -i 's/"//g' asv_table.tmp.xls
sed -i 's/"//g' asv_rep.fasta
perl {abs_bin}/sample_order.pl asv_table.tmp.xls {abs_raw} asv_table.xls
rm asv_table.tmp.xls

echo "First part of the script completed successfully."
"""
    f1.write(content1)

# Generate the second part script  (S01.1.pathogen.part2.sh)
if args.region == 'subregion':
    with open(mode.format("pathogen.part2") + '.sh', 'w', encoding='utf-8') as f2:
        content2 = f"""
    cd {PWD}01.ASV.tax/
    ####blast######
    makeblastdb -in {abs_db}/MFPD.ITS.fa -dbtype nucl -out {abs_db}/MFPD
    blastn -query asv_rep.fasta -db {abs_db}/MFPD -outfmt 6 -evalue 1e-5 -perc_identity {similarity} -num_threads 6 -out blastnM8.txt
    perl {abs_bin}/best_m8.pl blastnM8.txt > blastn.txt
    perl {abs_bin}/taxonomy_lenthre.pl blastn.txt {abs_db}/MFPD.ITS.tax asv_rep.temp.tax
    perl {abs_bin}/tax_normal.pl asv_rep.fasta asv_rep.temp.tax asv_rep.tax
    perl {abs_bin}/taxaTable_byAss.pl asv_table.xls asv_rep.tax asv_taxa_table.xls
    echo "Second part of the script completed successfully."
"""
        f2.write(content2)

elif args.region == 'full-length':
    with open(mode.format("pathogen.part2") + '.sh', 'w', encoding='utf-8') as f2:
        content2 = f"""
    cd {PWD}01.ASV.tax/
    ####hierarchical search######
    mmseqs easy-search asv_rep.fasta {abs_db}/clades_mmseq_db mmseq_search.txt temp --max-seqs 10 --min-seq-id 0.9   
    awk '{{print $2}}' mmseq_search.txt | sort | uniq > ./temp/mmseq_search_nodes.txt
    awk '{{print "../cluster/cluster_"$1".txt"}}' ./temp/mmseq_search_nodes.txt | xargs cat > ./temp/mmseq_search_subset.txt

    seqkit grep -f ./temp/mmseq_search_subset.txt {abs_db}/MFPD.ITS.fa > ./temp/mmseq_search_subset.fa
    makeblastdb -in ./temp/mmseq_search_subset.fa -dbtype nucl -out ./temp/mmseq_search_subset

    blastn -query asv_rep.fasta -db ./temp/mmseq_search_subset -out blastnM8.txt -outfmt 6 -max_target_seqs 3 -evalue 1e-10 -perc_identity {similarity} -num_threads 32
    perl {abs_bin}/blast_best.pl blastnM8.txt > blastn.txt
    perl {abs_bin}/taxonomy_lenthre.pl blastn.txt {abs_db}/MFPD.ITS.tax asv_rep.temp.tax
    perl {abs_bin}/tax_normal.pl asv_rep.fasta asv_rep.temp.tax asv_rep.tax  
    perl {abs_bin}/taxaTable_byAss.pl asv_table.xls asv_rep.tax asv_taxa_table.xls 
    echo "Second part of the script completed successfully."
"""
        f2.write(content2)

# Generate the dada2 script (would be used in firtst part) (S01.1.dada2.r)
with open(mode.format("dada2") + '.r', 'w', encoding='utf-8') as r_file:
    r_content = f"""
library(dada2)
path <- "{PWD}00.Data"
fns <- list.files(path, pattern=".fq.gz")
sample.names <- sapply(strsplit(basename(fns), "[.]"), `[`, 1)
filtpath <- file.path(path, "filtered")
filterAndTrim(file.path(path,fns), file.path(filtpath,fns),  maxEE=1, truncQ=2, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
filts <- list.files(filtpath, pattern="fq", full.names=TRUE)
names(filts) <- sample.names
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {{
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}}
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "./seqtab.rds")
write.table(t(seqtab),"table.xls",sep="\t")
seqtabrechim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
write.table(t(seqtabrechim),"table.rechim.xls",sep="\t")
"""
    r_file.write(r_content)

# Generate file link script (would be used in line 153)
with open(mode.format("symbolic_links") + '.sh', 'w', encoding='utf-8') as f_link:
    for i in range(num):
        f_link.write(f'ln -s {fq_list[i]} {PWD}00.Data/{ID_list[i]}.fq.gz\n')

# Defination run function
def run_shell_script(script):
    print(f"Running script: {script}")
    result = subprocess.run(['bash', script], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"[ERROR] Script failed: {script}")
        print("-------- STDERR --------")
        print(result.stderr)
        print("------------------------")
        sys.exit(f"[FATAL] Terminating pipeline due to error in {script}")
    else:
        print(result.stdout)
        print(f"[SUCCESS] Completed: {script}")

# Run script
print("Running prepartion of the script...")
run_shell_script(mode.format("symbolic_links") + '.sh')

print("Running first part of the script...")
run_shell_script(mode.format("pathogen.part1") + '.sh')

print("First part completed. Now running the second part...")
run_shell_script(mode.format("pathogen.part2") + '.sh')

print("All tasks completed successfully!")
progress_bar()