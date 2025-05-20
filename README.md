# MFPD: A multiple fungal pathogen detection pipeline for One Health practices

MFPD pipeline serves for automated processing and pathogen identification of fungal ITS sequencing data. Allow users input quality controled sequences then output pathogen ASVs abundance and taxonomy.

![image](https://github.com/LorMeBioAI/MFPD/blob/main/MFPD_graphic.png)

---


## Installation



```
git clone https://github.com/LorMeBioAI/MFPD

cd MFPD

conda env create -f environment.yml

conda activate mfpd-env
```


## Preparation for full-length version
```
sh prepare.sh
unzip clusters.zip
```


## Usage: example
```
python MFPD.py --file fq.list --pwd ./ --region subregion --identity 0.97
```

--file：sample table (tab delimiter)

--pwd：output direction (default: ./)

--region：targeted sequencing region， "subregion" or "full-length"

--identity：recommended identity: 0.99 for full-length and 0.97 for ITS1/2 , default=0.97



## Input file (tab delimiter)
A_1  {yourDir}/merge/SRR8146521.fastq.gz

A_2  {yourDir}/merge/SRR8146519.fastq.gz

A_3  {yourDir}/merge/SRR8146518.fastq.gz



## Output files


├── 00.Data/         # sample links

├── 01.ASV.tax/      # results folder

    ├── asv.table.xls：ASV tabel
    
    ├── asv_rep.fasta：ASV represent sequences
    
    
└── shell/           # Automatically generated and executed shell/r scripts


## Notes

MFPD pipeline employ DADA2 package，which can also carry out error learning and denoising even in the three generation sequencing, but it is recommended to input high-quality sequences that have undergone preliminary quality control.

In the full length mode , the mmseqs2  + blast hierarchical search strategy will be called to improve the efficiency.


## Contact

Yi Shen: shiny@stu.njau.edu.cn

Xinrun Yang: xinrunyang@stu.njau.edu.cn
