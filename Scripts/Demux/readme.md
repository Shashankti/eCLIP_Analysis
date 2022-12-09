## Required modules for demultiplexing

* Since the Undetermined files from the bcl2fastq conversion had the names of the barcodes present in the header
and not along with the sequence, the script `demuxbyname.sh` from BBMapSuite tools was used to generate fastq fils
which were split according to the names from barcode samples.

* The resulting fastq files were merged according to the lanes before proceeding with the analysis.

*Link to the tools used*

*Script used to run on hpc*
