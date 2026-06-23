##############################################################################
## 5. NanoFilt
##############################################################################


mkdir 5_nanofilt

#Pre-work 1:
#Check if your system has NanoFilt
module avail NanoFilt

#If not, you can build conda env and install the tool
#If you want, you can install NanoPlot and NanoFilt within the same env since we usually run them together


source ~/.bashrc
conda env list
conda activate nanoplot_env

NanoFilt -q 10 -l 1000 4_rawreads/b5.fastq  > 5_nanofilt/b5_filt.fastq 

# You can run seperately or write a loop