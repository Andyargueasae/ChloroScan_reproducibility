# Bash commands used in ChloroScan

The following commands are used to generate synthetic metagenomes, including assembling, mapping reads.
We also provide commands of running amber, metabat2, ChloroScan, FragGeneScanRs and orthofisher.


## 1. Generate Simulated Metagenomes:

 - Download prokaryote genomes from NCBI dataset 

```sh    
datasets download accession_id --include genome 
```

 - CAMISIM run on data.

We used ``CAMISIMv2`` which is a nextflow pipeline, link: https://github.com/CAMI-challenge/CAMISIM/tree/dev.
We firstly prepare the metagenomes based on the config files.
When running the main workflow, use the following command:
```sh
nextflow run main.nf 
```

Two synthetic metagenomes used in the study can be found at: https://doi.org/10.26188/28748540.

## 2. Mapping reads to assemblies for synthetic metagenomes and real metagenomes. 
Turn sam to sorted.bam files.

```sh
minimap2 -ax sr –1 {read1.fq.gz} –2 {read2.fq.gz} > sample_x.sam 
samtools view -bS sample_x.sam > sample_x.bam
samtools sort sample_x.bam -o sample_x.sorted.bam
```

## 3. Benchmarking ChloroScan: 

Run metabat2: 
1. generate metabat2 depth profiles from bam files provided
```sh
Jgi_summarize_bam_contig_depths --outputDepth /path/to/coverage.txt input1.sorted.bam input2.sorted.bam
```

2. Run MetaBAT2: 
```sh
metabat2 –i contigs.fasta -a /path/to/coverage.txt -o metabat2 –m 1500 –s 50000
```

ChloroScan commands for benchmarking datasets:
```sh
ASSEMBLY=/path/to/contigs.fasta
ALIGNMENT=/path/to/bam_dir
BATCH_NAME="your_batch_name"
OUTPUT=/path/to/output_dir
BINNY_OUTPUTDIR=/path/to/binny_output_directory

chloroscan run --Inputs-assembly=$ASSEMBLY \
    --Inputs-alignment=$ALIGNMENT \
    --Inputs-batch-name=$BATCH_NAME \
    --outputdir=$OUTPUT \
    --use-conda \
    --conda-prefix="conda_env" \
    --conda-frontend="mamba" \
    --cores=12 \
    --verbose \
    --corgi-pthreshold=0.5 \
    --binning-clustering-hdbscan-min-sample-range=1,5,10 \
    --binning-bin-quality-purity=90 \
    --binning-outputdir=$BINNY_OUTPUTDIR \
    --corgi-min-length=1000 \
    --force \
    --binning-bin-quality-min-completeness=50
```

For real metagenome runs, we altered the --binning-bin-quality-min-completeness from 50% to 70%. The command is:

```sh
ASSEMBLY=/path/to/contigs.fasta
ALIGNMENT=/path/to/bam_dir
BATCH_NAME="your_batch_name"
OUTPUT=/path/to/output_dir
BINNY_OUTPUTDIR=/path/to/binny_output_directory

chloroscan run --Inputs-assembly=$ASSEMBLY \
    --Inputs-alignment=$ALIGNMENT \
    --Inputs-batch-name=$BATCH_NAME \
    --outputdir=$OUTPUT \
    --use-conda \
    --conda-prefix="conda_env" \
    --conda-frontend="mamba" \
    --cores=12 \
    --verbose \
    --corgi-pthreshold=0.5 \
    --binning-clustering-hdbscan-min-sample-range=1,5,10 \
    --binning-bin-quality-purity=90 \
    --binning-outputdir=$BINNY_OUTPUTDIR \
    --corgi-min-length=1000 \
    --force \
    --binning-bin-quality-min-completeness=70
```

## 4. Run amber: 
AMBER is a binning benchmarking tool to compare the performance between 2 or more metagenomic binners with golden standard mapping information.  
```sh
amber.py -g gsa_mapping.tsv -l “binny,metabat2” binny_mapping.tsv metabat2_mapping.tsv -t 12  
```

# 5. Run FragGeneScanRs on resulting bins. 
After getting new bins from the sample SAMEA2732360 and SAMEA2189670, we planned to blast their rbcL genes. So we firstly predicted their genes via FragGeneScanRs: 

```sh
FragGeneScanRs -s bins.fa -t illumina_5 -p 4 -a bins_genes.faa -n bins_genes.fna
```

# 6. Run orthofisher to annotate rbcL genes.
This step has two file inputs:
 - bins.txt: a one/two column table storing paths to bins. 
 - hmms.txt: a one-column table storing path of rbcL.hmm.
Other parameters such as e-value, score and bitscore remained as default values.

```sh
orthofisher -f bins.txt -m hmms.txt -o output_directory
```

Results of identified rbcL sequences are saved in output_directory/all_sequences, and the long_summary.txt and short_summary.txt will tell the information of sequence-bin pairing. 
