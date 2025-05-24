#!/usr/bin/env bash

unzip clusters.zip
unzip db.zip

mkdir -p db
makeblastdb -in ./clusters/clades_95.fa -dbtype nucl -out ./db/clades_blastn_db
mmseqs createdb ./clusters/clades_95.fa ./db/clades_mmseq_db
mmseqs createindex ./db/clades_mmseq_db temp --search-type 3

# Create launcher script
cat << 'EOF' > MFPD
#!/usr/bin/env bash
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
python "$SCRIPT_DIR/MFPD.py" "$@"
EOF

chmod +x MFPD

# Link the launcher to conda's/bin
if [ -z "$CONDA_PREFIX" ]; then
    echo "[ERROR] Please activate your conda environment before running this script."
    exit 1
fi

ln -sf "$(realpath ./MFPD)" "$CONDA_PREFIX/bin/MFPD"
