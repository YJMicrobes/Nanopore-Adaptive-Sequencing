##############################################################################
## 14. GENOME COMPARISON
##############################################################################

mkdir -p 14_genome_comparison/

##########################
## 14.1 bins preparation
##########################

mkdir -p 14_genome_comparison/1_all_genomes

# Copy MAGs + isolates with clean naming

for f in MAG_bins/*.fa; do
    cp "$f" 14_genome_comparison/1_all_genomes/MAG_$(basename "$f")
done

for f in Isolates/*.fa; do
    cp "$f" 14_genome_comparison/1_all_genomes/ISO_$(basename "$f")
done


##########################
## 14.2 fastANI
##########################

cd 14_genome_comparison/1_all_genomes
ls *.fa > genomes.txt

source ~/.bashrc
conda activate fastANI_env

mkdir -p ../2_fastANI_out

fastANI \
    -ql genomes.txt \
    --rl genomes.txt \
    -o ../2_fastANI_out/all_vs_all_fastani.tsv


##########################
## 14.3 Genome abundance
##########################

mkdir -p 14_genome_comparison/3_abundance_out

module load minimap2/2.28
module load samtools/1.19.2

READS_DIR="../../5_nanofilt"


# IMPORTANT: run inside genome directory
GENOME_DIR="14_genome_comparison/1_all_genomes"
cd "$GENOME_DIR"

for GENOME in *.fa; do

    GENOME_NAME=$(basename "$GENOME" .fa)
    OUTFILE="../3_abundance_out/${GENOME_NAME}_abundance.tsv"

    echo "Processing ${GENOME_NAME}"

    ## index only if missing
    if [[ ! -f "${GENOME%.fa}.mmi" ]]; then
        minimap2 -d "${GENOME%.fa}.mmi" "$GENOME"
    fi

    if [[ ! -f "${GENOME}.fai" ]]; then
        samtools faidx "$GENOME"
    fi


    ## header
    echo -e "sample\ttotal_reads\tmapped_reads\tmapped_percent\tavg_coverage\tbreadth\tmapped_bases\tnorm_cov_per_Mreads" > "$OUTFILE"


    ## mapping loop
    for fq in "${READS_DIR}"/*.fastq; do

        sample=$(basename "$fq" .fastq)

        bam="${sample}.bam"
        sorted="${sample}.sorted.bam"

        minimap2 -ax map-ont "${GENOME%.fa}.mmi" "$fq" | \
            samtools view -b -o "$bam" -

        samtools sort -o "$sorted" "$bam"
        samtools index "$sorted"
        rm -f "$bam"

        total_reads=$(awk 'END{print NR/4}' "$fq")

        mapped_reads=$(samtools view -F 0x4 -c "$sorted")

        mapped_percent=$(awk -v m="$mapped_reads" -v t="$total_reads" \
            'BEGIN{if(t>0) printf "%.4f", (m/t)*100; else print 0}')

        read avg_cov breadth <<< $(samtools depth -a "$sorted" | \
            awk '{sum+=$3; if($3>0) covered++} END{if(NR>0) printf "%.4f %.4f", sum/NR, covered/NR; else print "0 0"}')

        mapped_bases=$(samtools view -F 0x4 "$sorted" | \
            awk '{sum+=length($10)} END{print sum+0}')

        norm_cov=$(awk -v ac="$avg_cov" -v tr="$total_reads" \
            'BEGIN{if(tr>0) printf "%.6f", ac/(tr/1e6); else print 0}')

        echo -e "${sample}\t${total_reads}\t${mapped_reads}\t${mapped_percent}\t${avg_cov}\t${breadth}\t${mapped_bases}\t${norm_cov}" >> "$OUTFILE"

    done

done


##########################
## 14.4 Distance tree: fasttree
##########################

source ~/.bashrc
conda activate gtdbtk_env
export GTDBTK_DATA_PATH=/isg/shared/databases/gtdb/v226.0/release226

mkdir -p 14_genome_comparison/3_gtdbtk_output

gtdbtk classify_wf \
  --genome_dir 14_genome_comparison/1_all_genomes/ \
  --out_dir 14_genome_comparison/3_gtdbtk_output \
  --extension fa \
  --skip_ani_screen \
  --cpus 20 \
  --force


mkdir -p 14_genome_comparison/4_fasttree/

gzip -dc 14_genome_comparison/3_gtdbtk_output/align/gtdbtk.bac120.user_msa.fasta.gz \
> 14_genome_comparison/4_fasttree/gtdbtk.bac120.user_msa.fasta

module load fasttree/2.1.10

FastTree -nt 14_genome_comparison/4_fasttree/gtdbtk.bac120.user_msa.fasta \
> 14_genome_comparison/4_fasttree/bac120_user_tree.nwk


##########################
## 14.5 Function comparison: prokka + eggNOG-mapper
##########################

source ~/.bashrc
conda activate prokka_env

input_dir=14_genome_comparison/1_all_genomes/
output_dir=14_genome_comparison/5_prokka_out

mkdir -p "$output_dir"

for g in "$input_dir"/*.fa; do

    name=$(basename "$g" .fa)

    prokka \
        --cpus 30 \
        --force \
        --outdir "$output_dir/$name" \
        --prefix "$name" \
        "$g"

done


module load eggnog-mapper/2.1.9
module load parallel/20240322
export EGGNOG_DATA_DIR=/isg/shared/databases/eggnog/v5.0.2

input_dir=14_genome_comparison/5_prokka_out
output_dir=14_genome_comparison/6_eggnog_out

mkdir -p "$output_dir"

export OMP_NUM_THREADS=10

for faa in "$input_dir"/*/*.faa; do

    base=$(basename "$faa" .faa)

    emapper.py \
        -i "$faa" \
        -o "$output_dir/${base}_eggnog" \
        --cpu 10 \
        -m diamond

done





