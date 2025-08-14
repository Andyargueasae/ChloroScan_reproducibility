#!/bin/bash

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=I_chloroscan

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G

# set your minimum acceptable walltime=days-hours:minutes:seconds
#SBATCH -t 00:40:00

#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --account=punim1293

# Specify your email address to be notified of progress.
#SBATCH --mail-user=yuhao.tong@student.unimelb.edu.au
#SBATCH --mail-type=ALL

module purge
module load GCCcore/11.3.0
module load Python/3.10.4

# # add anaconda3, see if conda can work.
# module load Anaconda3/2024.02-1
# eval "$(conda shell.bash hook)"

source /data/scratch/projects/punim1293/andy/chloroscan/.venv/bin/activate

# use the pip installed version of chloroscan
# pip install chloroscan
#export PATH="~/.cache/pypoetry/virtualenvs/chloroscan-C4-ifWe--py3.10/bin/chloroscan":$PATH

cd /home/yuhtong/scratch/andy/chloroscan

# Next assembly: /home/yuhtong/scratch/andy/punim1293/corgi_metag_protists_SAMEA2732613_sizefrac_20_180/assembly/SAMEA2732613-assembled.fa
ASSEMBLY="/home/yuhtong/scratch/andy/packed_262genome_8Gbp/contigs/anonymous_gsa.fasta"
ALIGNMENT="/home/yuhtong/scratch/andy/packed_262genome_8Gbp/bam"
# DEPTH=""
OUTPUT="/home/yuhtong/scratch/andy/packed_262genome_8Gbp/chloroscan_output_annot_2db_sc97_prefer_pfam"
BATCH_NAME="chloroscan_262"
BINNY_OUTPUTDIR="${BATCH_NAME}_binny_output_annot_pfam_sc97_prefer_pfam"
BINNY_DATABASE="/home/yuhtong/scratch/andy/marker_gene_database_90_95_97_annotation/sc97/2db_universal_0.6"

# prevent the use of other python packages except conda env pacakges.  
export PYTHONNOUSERSITE=1
export VIRTUAL_ENV=/home/yuhtong/scratch/andy/chloroscan/.venv
export PATH="$VIRTUAL_ENV/bin:$PATH"

echo $(which python)

chloroscan run \
    --Inputs-assembly=$ASSEMBLY \
    --Inputs-alignment=$ALIGNMENT \
    --Inputs-batch-name=$BATCH_NAME \
    --outputdir=$OUTPUT \
    --use-conda \
    --conda-prefix="conda_env" \
    --conda-frontend="mamba" \
    --cores=12 \
    --verbose \
    --corgi-pthreshold=0.5 \
    --binning-database=$BINNY_DATABASE \
    --binning-clustering-hdbscan-min-sample-range=1,5,10 \
    --binning-bin-quality-purity=90 \
    --binning-outputdir=$BINNY_OUTPUTDIR \
    --corgi-min-length=1000 \
    --force \
    --binning-bin-quality-min-completeness=50 $OUTPUT/working/binny

# db_path from binny (specified by chloroscan as binny-database) provides the taxon_marker_set_lineage_sorted.tsv. 
# While mantis.cfg gives the hmm file to annotate predicted genes.
