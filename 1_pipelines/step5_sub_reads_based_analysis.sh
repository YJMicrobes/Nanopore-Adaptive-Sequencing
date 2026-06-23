##############################################################################
## 5.1 Reads-based taxonomy annotation: Kaiju
##############################################################################

source ~/.bashrc
conda activate kaiju

# Input and output directories
input_dir=5_nanofilt/
output_dir=5.1_kaiju_reads/

mkdir -p "$output_dir"

# Kaiju database
kaiju_nodes=/scratch/bsteven/kaijudb/nodes.dmp
kaiju_db=/scratch/bsteven/kaijudb/refseq_nr/kaiju_db_refseq_nr.fmi


# Run Kaiju on filtered reads
for fastq_file in "$input_dir"/b*.fastq; do

    base_name=$(basename "$fastq_file" .fastq)

    output_file="$output_dir/${base_name}_kaiju.out"

    kaiju \
        -t "$kaiju_nodes" \
        -f "$kaiju_db" \
        -i "$fastq_file" \
        -o "$output_file" \
        -a greedy \
        -e 3 \
        -m 11 \
        -s 65 \
        -l 7 \
        -E 0.01

    echo "Processed $fastq_file -> $output_file"

done



# Convert Kaiju results to Krona format

for kaiju_file in "$output_dir"/*_kaiju.out; do

    base_name=$(basename "$kaiju_file" _kaiju.out)

    kaiju2krona \
        -t "$kaiju_nodes" \
        -n /scratch/bsteven/kaijudb/names.dmp \
        -i "$kaiju_file" \
        -o "$output_dir/${base_name}_krona.txt"

    echo "Converted $kaiju_file"

done


# Merge Krona files

ktImportText \
    "$output_dir"/*_krona.txt \
    -o "$output_dir/krona_visualization.html"

echo "Krona output: $output_dir/krona_visualization.html"




##############################################################################
## 5.2 Reads-based taxonomy annotation: NCBI BLAST
##############################################################################

module load blast/2.15.0

output_dir=5.2_ncbi_blast/

mkdir -p "$output_dir"


# BLAST against NCBI nt database
blastn \
    -query /scratch/bsteven/Colony_non_cleanReads/flye_output/assembly.fasta \
    -db nt \
    -remote \
    -out "$output_dir/ncbi_blast.out" \
    -outfmt 6 \
    -max_target_seqs 10 \
    -evalue 1e-5


# Extract accession numbers

cut -f2 "$output_dir/ncbi_blast.out" | \
    sort | uniq > "$output_dir/accessions.txt"



# Fetch taxonomy names

module load edirect/20250311


while read acc; do

    echo -n "$acc    " >> "$output_dir/taxonomy.txt"

    esearch \
        -db nucleotide \
        -query "$acc" | \
    efetch \
        -format docsum | \
    xtract \
        -pattern DocumentSummary \
        -element Title >> "$output_dir/taxonomy.txt"

    sleep 0.5

done < "$output_dir/accessions.txt"




##############################################################################
## 5.3 Reads-based mapping back to mosquito genome: minimap2 + samtools
##############################################################################

output_dir=5.3_mosquito_minimap/

mkdir -p "$output_dir"


# Build mosquito genome index

minimap2 \
    -d "$output_dir/mosquito.genome.mmi" \
    GCA_002204515.1_AaegL5.0_genomic.fna



# Map Nanopore reads

for fastq in 5_nanofilt/*.fastq; do

    base_name=$(basename "$fastq" .fastq)

    minimap2 \
        -ax map-ont \
        "$output_dir/mosquito.genome.mmi" \
        "$fastq" \
        > "$output_dir/${base_name}.sam"

done



module load samtools/1.19.2


# Convert SAM -> BAM, sort, index, statistics

for samfile in "$output_dir"/*.sam; do

    base_name=$(basename "$samfile" .sam)

    samtools view \
        -bS "$samfile" | \
    samtools sort \
        -o "$output_dir/${base_name}.bam"

    samtools index \
        "$output_dir/${base_name}.bam"

    samtools flagstat \
        "$output_dir/${base_name}.bam" \
        > "$output_dir/${base_name}.flagstat.txt"

done