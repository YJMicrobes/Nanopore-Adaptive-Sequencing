##############################################################################
## 12. GTDB-Tk classification (per group MAGs)
##############################################################################

source ~/.bashrc
conda create -n gtdbtk_env
conda activate gtdbtk_env
conda install bioconda::gtdbtk
gtdbtk --version

#https://ecogenomics.github.io/GTDBTk/
#download database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz  gtdb_db/
mv gtdbtk_data.tar.gz  gtdb_db/

tar -xvzf gtdbtk_data.tar.gz 

source ~/.bashrc
conda activate gtdbtk_env
export GTDBTK_DATA_PATH=/scratch/bsteven/database/gtdbtk_db/release226/
echo $GTDBTK_DATA_PATH

#but the database size is big to put in the individual HPC. So the admin installed the database and we can run like this:



mkdir -p 12_GTDB_out
source ~/.bashrc
conda activate gtdbtk_env
# GTDB database (shared server version preferred)
export GTDBTK_DATA_PATH=/isg/shared/databases/gtdb/v226.0/release226

groups=(colony_no_subtraction colony_subtraction H G T)

for g in "${groups[@]}"; do

    echo "Running GTDB-Tk for ${g}"


    ##########################################################################
    # Input MAGs from DAS Tool
    ##########################################################################

    GENOME_DIR=10_dastool/${g}/output/${g}_dastool_DASTool_bins
    OUT_DIR=12_GTDB_out/${g}

    mkdir -p "$OUT_DIR"



    ##########################################################################
    # Run GTDB-Tk
    ##########################################################################

    gtdbtk classify_wf \
        --genome_dir "$GENOME_DIR" \
        --out_dir "$OUT_DIR" \
        --extension fa \
        --cpus 48 \
        --skip_ani_screen \
        --force

done