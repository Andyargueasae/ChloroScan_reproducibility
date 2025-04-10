#Generate Simulated Metagenomes: 

# CAMISIM run on data: 
nextflow run main.nf 

# Mapping reads to assemblies:  
minimap2 -ax sr –1 {read1.fq.gz} –2 {read2.fq.gz} > sample_x.sam 

#Coassemble reads from the multisample synthetic metagenome: 
megahit -1 sample_0_01.fq.gz, sample_1_01.fq.gz, sample_2_01.fq.gz, sample_3_01.fq.gz -2 sample_0_02.fq.gz, sample_1_02.fq.gz, sample_2_02.fq.gz, sample_3_02.fq.gz -t 20 --out-dir megahit_coassembly_synethic --presets meta-large 

#Metaquast alignment between filtered contigs and source genomes: 
metaquast query_assembly.fasta -r source_genomes/ -t 12 –o output_dir 

#Run amber: 
amber.py -g gsa_mapping.tsv -l “binny,metabat2” binny_mapping.tsv metabat2_mapping.tsv -t 12  

#Benchmarking ChloroScan: 
# Run ChloroScan 
chloroscan run --Inputs-assembly=$ASSEMBLY --Inputs-alignment=$ALIGNMENT --Inputs-batch-name=$BATCH_NAME --outputdir=$OUTPUT --use-conda --cores=12 --verbose --corgi-pthreshold=0.5 --binning-clustering-hdbscan-min-sample-range=1,5,10 --binning-bin-quality-purity=90 --binning-outputdir=$BINNY_OUTPUTDIR --corgi-min-length=1500 --force --binning-bin-quality-min-completeness=50 $OUTPUT/working/binny 

#When running with real metagenomes, the corgi-pthreshold is set to 0.8 and the contig length cutoff is set to 1000.  

#Generate metabat2 depth profiles: 

Jgi_summarize_bam_contig_depths --outputDepth coverage.txt input1.sorted.bam input2.sorted.bam [...] 

#Run metabat2: 

metabat2 –i contigs.fasta -a coverage.txt -o metabat2 –m 1500 –s 50000  

#Download genomes from NCBI dataset 

datasets download accession_id --include genome 