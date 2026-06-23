
##############################################################################
## 1. upload and download data
##############################################################################
#login in HPC
#open local terminal 
ssh username@servername

#if need the public key setup to login in, check here:
https://cbc-uconn.github.io/hpc-docs/connect.html#ssh-keys
#follow the recipe and then open the terminal to login in


#upload
#scp pod5/* username@servername:/pod5/
scp *.pod5 username@servername:/pod5


#download
scp username@servername:/filename ~/folder

#rsync -av --ignore-existing release220/ /scratch/bsteven/gtdb_database/

#or use fillzella
#servername
#usename
#code
#port: 22