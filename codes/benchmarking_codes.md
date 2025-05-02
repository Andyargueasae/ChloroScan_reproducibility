# Bash commands used in ChloroScan

The following commands are used to generate synthetic metagenomes, including assembling, mapping reads.
And the commands of amber, metaquast and metabat2.


## 1. Generate Simulated Metagenomes: 
Download genomes from NCBI dataset 
```sh    
datasets download accession_id --include genome 
```

## 2. CAMISIM run on data.

We used ``CAMISIMv2`` which is a nextflow pipeline, link: https://github.com/CAMI-challenge/CAMISIM/tree/dev.
We firstly prepare the metagenomes based on the config files.
When running the main workflow, use the following command:
```sh
nextflow run main.nf 
```

## 3. Mapping reads to assemblies. 
Turn sam to sorted.bam files.

```sh
minimap2 -ax sr –1 {read1.fq.gz} –2 {read2.fq.gz} > sample_x.sam 
samtools view -bS sample_x.sam > sample_x.bam
samtools sort sample_x.bam -o sample_x.sorted.bam
```

## 4. Coassemble reads from the multisample synthetic metagenome: 
```sh
megahit -1 sample_0_01.fq.gz,sample_1_01.fq.gz,sample_2_01.fq.gz,sample_3_01.fq.gz \
-2 sample_0_02.fq.gz, sample_1_02.fq.gz, sample_2_02.fq.gz, sample_3_02.fq.gz -t 20 --out-dir \ megahit_coassembly_synethic --presets meta-large 
```

## 5. Metaquast alignment between filtered contigs and source genomes: 
```sh
metaquast path/to/query.fasta -r path/to/source_genomes/ -t 12 –o output_dir 
```

## 6. Run amber: 
```sh
amber.py -g gsa_mapping.tsv -l “binny,metabat2” binny_mapping.tsv metabat2_mapping.tsv -t 12  
```

## 7. Benchmarking ChloroScan: 
Run ChloroScan for synthetic metagenomes:

```sh
chloroscan run --Inputs-assembly=$ASSEMBLY --Inputs-alignment=$ALIGNMENT --Inputs-batch-name=$BATCH_NAME --outputdir=$OUTPUT --use-conda --cores=12 --verbose --corgi-pthreshold=0.5 --binning-clustering-hdbscan-min-sample-range=1,5,10 --binning-bin-quality-purity=90 --binning-outputdir=$BINNY_OUTPUTDIR --corgi-min-length=1500 --force --binning-bin-quality-min-completeness=50 $OUTPUT/working/binny 
```
note: When running with real metagenomes, the corgi-pthreshold is set to 0.8 and the contig length cutoff is set to 1000.  

## 8. metabat2: 
1. generate metabat2 depth profiles from bam files provided
```sh
Jgi_summarize_bam_contig_depths --outputDepth /path/to/coverage.txt input1.sorted.bam input2.sorted.bam
```

2. Run metabat2: 
```sh
metabat2 –i contigs.fasta -a /path/to/coverage.txt -o metabat2 –m 1500 –s 50000
```
