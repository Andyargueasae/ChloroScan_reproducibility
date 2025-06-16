# Marker gene database construction process. 

## 1. Run OrthoFinder:
``Input data type``: 458 fasta files storing proteins for each gene stored in a directory. 
``output directory``: orthofinder's results containing all orthogroups and basic information summary, such as the sinle-copy orthogroups.

```sh
orthofinder -f /path/to/458_genomes_selected \
-t 12  -n “job name” -og 
```
## 2. Identify orthogroups containing desired genes, which are saved in ortho_A2K folder, under the directory /ortho_A2K/Orthogroups_selected.

To simplify the downstream process, I renamed each sequence's header to "species_name|protein_id|gene_name".

```sh
for i in ortho_A2K/Orthogroups_selected/*.fasta; do
    python rename_orthogroup.py $i backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/selected_orthogroups/$(basename $i .fasta).fasta
```

Then, some orthogroups might have inparalogs which are identical to each other, try to deduplicate orthogroups using dedup.sh.

```sh
bash dedup.sh --fpath backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/selected_orthogroups --opath backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/intermediary_files_for_preprocessing_alignments/deduplicated_genes --scriptpath backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/scripts_used_for_preprocessing_orthogroup_sequences
```

``Note``: some orthogroups have genes annotated with two different sequences, which are not inparalogs, we deduplicated them further by retaining only the longer seqeunce and labelled them as "dedup2". This is a manual process, and the results are provided. 

## 3. Process input sequences into alignment using mafft. (need a for loop)

```sh
alignment_dir="backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/intermediary_files_for_preprocessing_alignments/deduplicated_aligned_genes"

for i in $(ls backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/intermediary_files_for_preprocessing_alignments/deduplicated_genes); do
    mafft --auto --threads 8 $i > backup_for_manuscript1_IMPORTANT/section_1_marker_gene_database/intermediary_files_for_preprocessing_alignments/deduplicated_aligned_genes/$(basename $i .fasta).aln   
done
```

remove the header for each sequence that only retains taxon name.

## 4. trim the each gene's alignment using clipkit. 

```sh
for i in $(ls deduplicated_aligned_taxon_only/*.fa); do    
    clipkit $i -m smart-gap -o trimmed_deduplicated_aligned_taxon_only/trimmed_$(basename $i .fa)fa
done
```

## 5. Construct supermatrix (concatenate 22 genes’ alignments into one) 

collect links of all alignments into a file called alignment_list.tsv.
```sh
phykit create_concatenation_matrix --alignment alignment_list.tsv --prefix A2K.phylo 
```
One of the output: A2K.phylo.fa will be the concatenated alignment of 22 genes and will be the input of iqtree2.  

## 6. Infer a maximum Likelihood Tree from the constructed supermatrix.
```sh
iqtree -s /path/to/A2K.phylo.fa -bb 1000 -m TEST -nt 8 -redo -pre A2K
```

After running the codes above, a phylogenetic tree for 458 plastid genomes is generated. 
Next, following checkm's rationale (https://pmc.ncbi.nlm.nih.gov/articles/PMC4484387/), the tree is decorated with marker genes for each internal node.   

## 7. The lineage-specific marker set calculation for each internal node of the tree from iqtree2 outputs: 

```sh
python construct_marker_set.py --input-tree "/path/to/A2K.tax_mod.rerooted.reannotated.treefile"  --node-wise True --output-tree marker_sets_verify.tree --species-genome-dict /path/to/species_genome_effective_dict.pkl 
```

## 8. Then, the marker genes will be collocated into marker sets, following settings in original checkm.  

```sh
python taxon_annotation_ms.py --taxon-list "taxon_list.txt" --input-tree "marker_sets_verify.tree" --endosymbiosis-map "./Endosymbiosis_dict.pkl" --output "./taxon_marker_sets_verify.tsv"
```

The marker sets we have now should contain a universal marker gene set to be completely functional. So we calculate the marker set containing single-copy marker genes found in > 95% of genomes and colocalize them into marker sets.  

It is already computed as when doing node-wise computations the node that contains all genomes as leaves is included, we only need to colocalize them into marker sets and store them in a tsv file.
```sh
python taxon_annotation_ms.py --universal-marker-set True --input-tree /path/to/A2K.tax_mod.rerooted.reannotated.treefile --universal-marker-set-out "/path/to/universal_marker_set.tsv"
```
The marker set is then concatenated at the top of the other marker sets.

We also included a list of 34 genes proposed by Janouškovec et al. (2010), in comparsion with our marker sets in estimating plastid qualities.

Now the taxon_marker_sets.tsv file is finished, the left step is to create profile hmms. 

## 9. Align included genes.

After identifying marker genes that will be used in the database, we align the marker genes identified in the orthogroups (A manual selection).

```sh
for selected_OGs in $(ls selected_OGs); do
    mafft --auto $selected_OGs > /path/to/marker_gene_alignments/gene_name.aln
done
```

## 10. Build profile HMM.
Then, we create the hmm file for these genes using hmmer3.
```sh
bash scripts_used_for_preprocessing_orthogroup_sequences/hmmbuild.sh --aln_dir path/to/aln_dir --hmm_dir path/to/hmm_dir
```

Finally, create a directory called database with such structure:
binny_Chloroscan/A2K_database

```sh
├── hmms
│   └── checkm_pf
│       ├── checkm_filtered_pf.hmm
│       ├── checkm_filtered_pf.hmm.h3f
│       ├── checkm_filtered_pf.hmm.h3i
│       ├── checkm_filtered_pf.hmm.h3m
│       ├── checkm_filtered_pf.hmm.h3p
│       └── chunks
│           ├── checkm_filtered_pf_chunk_0.hmm
│           ├── checkm_filtered_pf_chunk_0.hmm.h3f
│           ├── checkm_filtered_pf_chunk_0.hmm.h3i
│           ├── checkm_filtered_pf_chunk_0.hmm.h3m
│           └── checkm_filtered_pf_chunk_0.hmm.h3p
├── pfam
│   └── tigrfam2pfam.tsv
└── taxon_marker_sets_lineage_sorted.tsv
```

Here, ``taxon_marker_sets_lineage_sorted.tsv`` is the tsv file that provides all marker sets available to use, and the ``hmms`` directory stores combined hmms in ``checkm_filtered_pf_chunk_0.hmm``. We retained the names of those files from the original database used by the original version of binny that targets prokaryote MAGs. The folder ``pfam`` contains the mapping information of PFAM families to TIGRFAM families from the original CheckM database, it is no longer used.   

Move the directory into binny's directory (the version used in ChloroScan). Finally, change the binny_mantis.cfg file to enable binny's workflow to use this database. 