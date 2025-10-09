# ChloroScan_code
Codes used in the ChloroScan benchmarking and data preparations. 

 - Bash commands used in generating synthetic metagenomes and benchmarking ChloroScan: ``codes/bash_codes.md``. Codes such as downloading NCBI dataset genomes, running binning benchmarking and running binning tools (MetaBAT2/ChloroScan) are here.

 ``Note``: ChloroScan v0.1.5 is used here for generating these results, the command to download is: ``pip install chloroscan==0.1.5``.

 - Codes used in generating marker gene database: ``codes/marker_gene_database_commands.md``.
 
 - OrthoFinder outputs and the scripts and files involved in generating the marker gene database are dumped here: https://doi.org/10.26188/28722788. A copy of result files and intermediary files for generating marker gene database is in ``marker_gene_database_complete.tar.gz``. Before running codes in ``./codes/marker_gene_database_commands.md``, please download the files from figshare. The command to download is:

```sh
    wget --referer=https://figshare.unimelb.edu.au --user-agent="Mozilla/5.0" -O "marker_gene_database_complete.tar.gz" https://figshare.unimelb.edu.au/ndownloader/files/58584934

```
