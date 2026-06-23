##############################################################################
## 4. nanoplot quality check
##############################################################################


mkdir 3_nanoplot

#Pre-work 1:
#check the barcodes that used in the experiment (metadata)
#Match the barcodes.fastq files with the experiment metadata
#Rename the .fastq files and save them as the raw .fastq reads
#Delete the unknown barcodes.fastq files

mkdir 4_rawreads 


#Pre-work 2:
#Check if your system has NanoPlot
module avail NanoPlot
#If not, you can build conda env and install the tool
conda install NanoPlot


source ~/.bashrc
conda env list
conda activate nanoplot_env

NanoPlot \
  --fastq 4_rawreads/b5.fastq \
  -o 3_nanoplot/ \
  --prefix b5

# You can run one by one or merge them as one to run or write a loop
# In the end, you can download the report to check the result