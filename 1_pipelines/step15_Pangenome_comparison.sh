##############################################################################
## 15. PAN-GENOME COMPARISON (Comamonas)
## Goal: reconstruct gene-level + genome-level comparison + circular genome plot
##############################################################################

BASE=15_pangenome_comparison

##############################################################################
## 15.1 SET INPUT WORKSPACE
## (organize genomes and create working structure)
##############################################################################

mkdir -p ${BASE}/1_Comamonas
cd ${BASE}/1_Comamonas

# copy genomes into working directory (standardize input)
cp ../../MAG_bins/2_metacoag_bin_2.fa .
cp ../../Isolates/hgt3_b16.fa .
cp ../../Isolates/hgt3_b20.fa .
cp ../../Isolates/hgt3_b2.fa .


##############################################################################
## 15.2 GENE PREDICTION (PRODIGAL)
## (extract CDS proteins for pangenome analysis)
##############################################################################

mkdir -p 1_prodigal
module load prodigal/2.6.3

prodigal -i 2_metacoag_bin_2.fa -a 1_prodigal/MAG7.faa -d 1_prodigal/MAG7.fna -p single
prodigal -i hgt3_b16.fa         -a 1_prodigal/B16.faa  -d 1_prodigal/B16.fna  -p single
prodigal -i hgt3_b20.fa         -a 1_prodigal/G20.faa  -d 1_prodigal/G20.fna  -p single
prodigal -i hgt3_b2.fa          -a 1_prodigal/B2.faa   -d 1_prodigal/B2.fna   -p single


##############################################################################
## 15.3 PROTEIN TAGGING + MERGING
## (add genome IDs so clustering can track gene origin)
##############################################################################

cd 1_prodigal

sed 's/^>/\>MAG7_/' MAG7.faa > MAG7_tagged.faa
sed 's/^>/\>B16_/'  B16.faa  > B16_tagged.faa
sed 's/^>/\>G20_/'  G20.faa  > G20_tagged.faa
sed 's/^>/\>B2_/'   B2.faa   > B2_tagged.faa

cat *_tagged.faa > all_comamonas.faa


##############################################################################
## 15.4 GENE CLUSTERING (MMSEQS2)
## (define orthologous gene families across genomes)
##############################################################################

module load mmseqs2

mmseqs createdb all_comamonas.faa comDB

mmseqs cluster comDB clusterDB tmp \
    --min-seq-id 0.85 \
    -c 0.8 \
    --cov-mode 1

mmseqs createtsv comDB comDB clusterDB clusters.tsv
mmseqs result2flat comDB comDB clusterDB rep_sequences.faa


##############################################################################
## 15.5 GENOME ANNOTATION (PROKKA)
## (standard annotation for panaroo input)
##############################################################################

cd ..

mkdir -p 2_prokka
module load prokka

for g in *.fa; do
    name=$(basename "$g" .fa)

    prokka \
        --cpus 30 \
        --force \
        "$g" \
        --outdir 2_prokka/${name} \
        --prefix ${name}
done


##############################################################################
## 15.6 PANGENOME ANALYSIS (PANAROO)
## (build core/accessory gene matrix across genomes)
##############################################################################

mkdir -p 3_panaroo

source ~/.bashrc
conda activate panaroo_env

panaroo \
    -i 2_prokka/*/*.gff \
    -o 3_panaroo \
    -t 10 \
    --clean-mode moderate \
    -c 0.90


##############################################################################
## 15.7 WHOLE GENOME ALIGNMENT + CIRCOS (MAUVE → CIRCOS)
## (visualize genome rearrangements + synteny)
##############################################################################

mkdir -p 5_circos
cd 5_circos

# copy genome FASTA files
cp ../*.fa .
cp ../2_prokka/*.fna .


##############################################################################
## STEP 1: MULTIPLE GENOME ALIGNMENT (progressiveMauve)
## (detect syntenic blocks across genomes)
##############################################################################

source ~/.bashrc
conda activate mauve_env

progressiveMauve \
  --output=4genomes.xmfa \
  *.fna


##############################################################################
## STEP 2: XMFA → CIRCOS LINKS
## (convert alignment blocks into circos-compatible links)
##############################################################################


# This script converts progressiveMauve XMFA output into Circos links format

python step15.1_xmfa_to_links.py


##############################################################################
## STEP 3: FILTER LINKS
## (remove small or weak alignments)
##############################################################################

awk '$3-$2 > 5000' links.txt > links.filtered.txt


##############################################################################
## STEP 4: KARYOTYPE (pre-defined file)
## (defines genome layout for Circos plot)
##############################################################################

# static file (do not regenerate each run)
cp karyotype.txt 5_circos/


##############################################################################
## STEP 5: CIRCOS CONFIG (pre-defined file)
## (controls layout and visualization parameters)
##############################################################################

# copy or ensure config exists in working directory
cp 5_circos/circos.conf .



##############################################################################
## STEP 6: RUN CIRCOS
## (generate final circular genome comparison figure)
##############################################################################

source ~/.bashrc
conda activate circos_env

circos -conf circos.conf

echo "[DONE] Comamonas pangenome + circos pipeline completed"




############################################################################
### Bacillus - Same process 
############################################################################

############################################################################
### Other works i did: 
############################################################################


python step15.1_gene_presence.py
# gene_presence_absence_matrix.csv


#Visualize the presence/absence heatmap
python step15.2_gene_presence_heatmap.py


#I tried skani to calculate the ANI but the results are similar to the fastANI so we didn't end up using it.
source ~/.bashrc
conda create -n skani_env python=3.9 -y
conda activate skani_env
conda install -c bioconda skani

mkdir hgt_skani
skani triangle HGT_MAGs_genomes_fa/* -o hgt_skani/hgt_ani.tsv

skani search HGT_MAGs_genomes_fa/* \
             -d HGT_MAGs_genomes_fa/ \
             -o hgt_skani/hgt_ani_detailed.tsv \

#I tried mummer2circos to make the circo figures but it is not friendly to edit. 
source ~/.bashrc
conda env list
conda activate mummer2circo_env
mummer2circos -l \
    -r B16.fna \
    -q B16.fna B2.fna G20.fna MAG7.fna \
    -f -fr

