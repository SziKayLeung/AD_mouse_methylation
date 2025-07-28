# Isabel Castanho # I.Castanho@exeter.ac.uk

# Transfer RRBD files

nohup wget -r LINK & # example: nohup wget -r ftp://Project_2605:fuln29mocZ8dh@ftp.sequencing.exeter.ac.uk/ &

ps #optional

more nohup.out #optional




# transfer cuffdiff2 (2 samples removed)
# Sample P22 removed due to contamination and low alignment rate. See krona plot.
# Sample L19 removed due to low yield

/mnt/data1/isabel/RRBS

nohup wget -r ftp://Project_2605:fuln29mocZ8dh@ftp.sequencing.exeter.ac.uk/ &

ps #optional

more nohup.out #optional

mv /mnt/data1/isabel/RRBS/ftp.sequencing.exeter.ac.uk/01_miseq_nano/ /mnt/data1/isabel/RRBS/01_miseq_nano/