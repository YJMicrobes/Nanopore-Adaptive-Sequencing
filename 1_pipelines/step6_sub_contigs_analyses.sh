##############################################################################
## 6.1 Contigs-based taxonomy annotation: kraken
##############################################################################


module load kraken/2.1.2


kraken_db=/isg/shared/databases/kraken/v06.2025/standard

assembly_dir=6_flye_output
output_dir=7_kraken2_contigs

mkdir -p "$output_dir"



for group in colony_no_subtraction colony_subtraction H G T

do

    echo "Running Kraken2 for ${group}"


    mkdir -p "$output_dir/${group}"


    kraken2 \
      --db "$kraken_db" \
      --output "$output_dir/${group}/kraken2_output.txt" \
      --report "$output_dir/${group}/kraken2_report.txt" \
      --use-names \
      --confidence 0.1 \
      --threads 8 \
      "$assembly_dir/${group}/assembly.fasta"


    echo "Finished ${group}"

done



#Step 2 merge the output
python step6.1.2_kraken_summary.py 





##############################################################################
## 6.2 Contigs-based taxonomy checking: NCBI
##############################################################################

##############################################################################
## 6.2 Contigs-based taxonomy checking: NCBI BLAST
##############################################################################

module load blast/2.15.0
module load edirect/20250311


assembly_dir=6_flye_output
output_dir=6.2_ncbi_blast


mkdir -p "$output_dir"



for group in colony_no_subtraction colony_subtraction H G T

do

    echo "Running NCBI BLAST for ${group}"


    mkdir -p "$output_dir/${group}"



    ##########################################################################
    # Step 1. BLAST against NCBI nt database
    ##########################################################################

    blastn \
        -query "$assembly_dir/${group}/assembly.fasta" \
        -db nt \
        -remote \
        -out "$output_dir/${group}/ncbi_blast.out" \
        -outfmt 6 \
        -max_target_seqs 10 \
        -evalue 1e-5



    ##########################################################################
    # Step 2. Extract accession numbers
    ##########################################################################

    cut -f2 \
        "$output_dir/${group}/ncbi_blast.out" | \
        sort | uniq \
        > "$output_dir/${group}/accessions.txt"



    ##########################################################################
    # Step 3. Retrieve descriptions from NCBI
    ##########################################################################

    while read acc;

    do

        echo -n "${acc}\t" \
            >> "$output_dir/${group}/taxonomy.txt"


        esearch \
            -db nucleotide \
            -query "$acc" | \
        efetch \
            -format docsum | \
        xtract \
            -pattern DocumentSummary \
            -element Title \
            >> "$output_dir/${group}/taxonomy.txt"


        sleep 0.5


    done < "$output_dir/${group}/accessions.txt"



    echo "Finished ${group}"

done




##############################################################################
## 6.3 Contigs-based mapping back to mosquito genome: minimap2 + samtools (host contamination check)
##############################################################################

module load minimap2
module load samtools/1.19.2


output_dir=6.3_mosquito_minimap
mkdir -p "$output_dir"



# -------------------------------------------------------
# Step 1. Build mosquito genome index
# -------------------------------------------------------

ref=GCA_002204515.1_AaegL5.0_genomic.fna

minimap2 \
    -d "$output_dir/mosquito.genome.mmi" \
    "$ref"



# -------------------------------------------------------
# Step 2. Map reads (ALL groups, still reads-based step)
# -------------------------------------------------------

for fastq in 5_nanofilt/*.fastq; do

    base_name=$(basename "$fastq" .fastq)

    echo "Mapping $base_name to mosquito genome"

    minimap2 \
        -ax map-ont \
        "$output_dir/mosquito.genome.mmi" \
        "$fastq" \
        > "$output_dir/${base_name}.sam"

done



# -------------------------------------------------------
# Step 3. SAM → BAM processing
# -------------------------------------------------------

for samfile in "$output_dir"/*.sam; do

    base_name=$(basename "$samfile" .sam)

    samtools view -bS "$samfile" | \
    samtools sort -o "$output_dir/${base_name}.bam"

    samtools index "$output_dir/${base_name}.bam"

    samtools flagstat "$output_dir/${base_name}.bam" \
        > "$output_dir/${base_name}.flagstat.txt"

done




##############################################################################
## 6.4 Contigs-based function annotation: prokka + eggNOG-mapper
##############################################################################

source ~/.bashrc
conda activate prokka_env


input_dir=6_flye_output
output_dir=6.4_prokka_out

mkdir -p "$output_dir"



groups=(no ada h g t)


for g in "${groups[@]}"; do

    echo "Running Prokka: $g"

    prokka \
        --cpus 30 \
        --force \
        "$input_dir/${g}/assembly.fasta" \
        --outdir "$output_dir/${g}" \
        --prefix "$g"

done


module load eggnog-mapper/2.1.9
module load parallel/20240322
export EGGNOG_DATA_DIR=/isg/shared/databases/eggnog/v5.0.2

input_dir=6.4_prokka_out
output_dir=6.5_eggnog_out

mkdir -p "$output_dir"
export OMP_NUM_THREADS=10
groups=(no ada h g t)


for g in "${groups[@]}"; do

    echo "Running eggNOG: $g"

    emapper.py \
        -i "$input_dir/$g/${g}.faa" \
        -o "$output_dir/${g}/${g}_eggnog" \
        --cpu 10 \
        -m diamond \
        --output_dir "$output_dir/${g}"

done


