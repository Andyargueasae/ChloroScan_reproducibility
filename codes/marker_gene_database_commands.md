# Marker gene database construction process. 

``Note``: All essential files used in this section can be found in ``marker_gene_database_files.tar.gz``, in the same github repository.

First, clone the repository and unzip the tar.gz file. 

```sh
git clone https://github.com/Andyargueasae/ChloroScan_reproducibility.git
cd ChloroScan_reproducibility/

wget --referer=https://figshare.unimelb.edu.au --user-agent="Mozilla/5.0" \
    -O "marker_gene_database_complete.tar.gz" https://figshare.unimelb.edu.au/ndownloader/files/57355894

tar zxvf marker_gene_database_complete.tar.gz

SCRIPT_DIR="./marker_gene_database/scripts_used_for_preprocessing_orthogroup_sequences"
PYTHON_SCRIPT_DIR="./marker_gene_database/python_scripts_for_generating_marker_sets"

# Make a directory to save all intermediary files in here. 
mkdir -p reproducibility
```

Previous exploration using the 913 plastid genomes of unique species found 22 genes: atpA, atpB, atpH, petB, psaB, psaC, psbA, psbE, psbH, rpl14, rpl16, rpl2, rpl20, rpoA, rpoC1, rps11, rps12, rps14, rps19, rps3, rps7, rps8. 

 - These genes are present in >= 97% of genomes, so we chose them to conduct the downstream analyses. We used these genes to filter the genomes, and finally chose 458 genomes (with one genome per genus) to conduct downstream analyses. 

References can be found in ``marker_gene_database/data_exploration_plastids/A2K.ipynb``.

A table mapping selected genbank id to the species name is in ``marker_gene_database/list_of_genbanks_used_in_tree.txt``.

## 1. Run OrthoFinder:

``Input data type``: 
 - 458 fasta files storing proteins for each gene stored in a directory.
 
``output data type``: 
 - orthofinder's results containing all orthogroups and basic information summary, such as the sinle-copy orthogroups.

```sh
orthofinder -f ChloroScan_reproducibility/marker_gene_database/458_genomes_selected \
-t 20  -n "plastid_marker_gene_db" -og 
```
The runtime is about 3 hours. A backup result using the same input genomes is dumped in figshare: https://figshare.unimelb.edu.au/articles/dataset/ChloroScan_Supplementary_files/28722788, with the name ortho_A2K.tar.gz. 

## 2. Identify orthogroups containing desired genes, preprocess them by renaming headers and deduplications.

``Input data type``: 
 - each selected orthogroup's fasta file.

``output data type``: 
 - selected orthogroup's fasta file with header in "species_name|protein_id|gene_name".

```sh

mkdir -p reproducibility/selected_orthogroups

for i in ./marker_gene_database/Orthogroups_phylo/*.fasta; do
    python $SCRIPT_DIR/rename_orthogroup.py $i ./reproducibility/selected_orthogroups/$(basename $i .fasta).fasta
done
```

Then, some orthogroups might have inparalogs which are identical to each other, try to deduplicate orthogroups using dedup.sh.

``Input data type``: 
 - each renamed and selected orthogroup's fasta file.

``output data type``: 
 - deduplicated genes, with inparalogs and genes of same annotation dereplicated to retain one sequence per gene.

```sh

mkdir -p reproducibility/deduplicated_genes

bash $SCRIPT_DIR/dedup.sh --fpath reproducibility/selected_orthogroups --opath reproducibility/deduplicated_genes --scriptpath $SCRIPT_DIR
```

``Note``: some orthogroups have genes annotated with two different sequences, which are not inparalogs, we deduplicated them further by retaining only the longer seqeunce and labelled them as "dedup2".

```sh 
NUM_GENOMES=458
mkdir -p reproducibility/deduplicated_genes2
for file in $(ls reproducibility/deduplicated_genes); do
    if [[ $(cat reproducibility/deduplicated_genes/$file | grep -c ">") -gt $NUM_GENOMES ]]; then
        # do something.
        echo "$file has extra copies, potentially the same gene annotated with two sequences."
        echo "Should remove the shorter one".
        python $SCRIPT_DIR/remove_extra_annotations.py -i reproducibility/deduplicated_genes/$file -o reproducibility/deduplicated_genes2/$file
    else
        echo "$file has matching number of entries, skip."
        cp reproducibility/deduplicated_genes/$file reproducibility/deduplicated_genes2/$file
    fi
done
```

## 3. Process input sequences into alignment using mafft, and rename the sequence headers again.

``Input data type``: 
 - renamed, deduplicated genes' fasta file.

``Output data type``: 
 - alignment file produced by mafft.

``Package``: 
 - Mafft v7.526 (2024/Apr/26).


```sh
mkdir -p reproducibility/deduplicated_aligned_genes
for i in $(ls reproducibility/deduplicated_genes2); do
    mafft --auto --thread 8 reproducibility/deduplicated_genes2/$i > reproducibility/deduplicated_aligned_genes/$(basename $i .fasta).aln
done
```

rename the header for each sequence that only retains taxon name.

```sh
mkdir -p reproducibility/renamed_alignments 

for file in $(ls reproducibility/deduplicated_aligned_genes/); do
    python $SCRIPT_DIR/rename_alignment.py reproducibility/deduplicated_aligned_genes/$file reproducibility/renamed_alignments/$file
done
```
Doing steps above is to make sure that no bugs will arise in step 4, as duplicated genes and uneven alignment length causes phykit create_concatenation_matrix to fail.

## 4. Construct supermatrix (concatenate 22 genes’ alignments into one) 

``Input data``: 
 - a tsv file storing 22 genes' alignment file link. 

``Output data``: 
 - supermatrix-the concatenated alignment from 22 genes' aln file. 

``Package``:
 - phykit v1.19.8.

```sh
# prepare alignment_list.tsv.
mkdir -p reproducibility/supermatrix

for file in reproducibility/renamed_alignments/*;do
    realpath $file;
done | sed 's/\x1b\[[0-9;]*[a-zA-Z]//g' > reproducibility/alignment_list.tsv

phykit create_concatenation_matrix --alignment reproducibility/alignment_list.tsv --prefix reproducibility/supermatrix/A2K.phylo 
```

One of the output: A2K.phylo.fa will be the concatenated alignment of 22 genes and will be the input of iqtree2.  

## 5. Infer a maximum Likelihood Tree from the constructed supermatrix.

``Input data``: 
 - concatenated alignment of 22 genes from 458 input genomes. 

``Output data``: 
 - a maximum-likelihood tree infering phylogenetic relationships of 458 input genomes.  

``Package``:
 - IQ-TREE v2.3.6.

```sh
mkdir -p reproducibility/tree
iqtree -s reproducibility/supermatrix/A2K.phylo.fa -bb 1000 -m TEST -nt 8 -redo -pre reproducibility/tree/A2K
```
The whole process will take ca. 5 hours to finish.

After running the codes above, a phylogenetic tree for 458 plastid genomes is generated. 
Next, following CheckM's rationale (https://pmc.ncbi.nlm.nih.gov/articles/PMC4484387/), the tree is decorated with marker genes for each internal node.

Before going to single copy marker gene calculation, the tree has to be modified by: renaming branches to add taxonomic lineages similar to greengene format: kingdom_phylum_class_order_family_genus_species, this simplifies the identification of node's taxon. The modified tree is provided with the link: ChloroScan_reproducibility/marker_gene_database/taxon_annotated_treefiles_and_marker_set_decorated_tree_files/A2K.tax_mod.rerooted.reannotated.treefile. 

Here we provide the pdf file visualizing the topology for this file: ChloroScan_reproducibility/A2K_tree_structure.pdf. 

## 6. The lineage-specific marker set calculation for each internal node of the tree from iqtree2 outputs: 

``Input data``:
 - ``treefile``: the resulting maximum likelihood tree.
 - ``species_genome_effective_dict.pkl``: a python dictionary pairing genomes with gene contents and genbank id for each source genomes.  

``Output data``:
 - ``output tree``: the original tree file decorated with single-copy marker genes for each internal node.   

``note``: the tree should be firstly preprocessed and rerooted by setting Paulinella as the root. 

```sh
mkdir -p reproducibility/marker_gene_database

python $PYTHON_SCRIPT_DIR/construct_marker_set.py --input-tree "marker_gene_database/taxon_annotated_treefiles_and_marker_set_decorated_tree_files/A2K.tax_mod.rerooted.reannotated.treefile"  --node-wise True --output-tree reproducibility/marker_gene_database/marker_sets_verify.tree --species-genome-dict marker_gene_database/required_files_for_generating_marker_gene_database/species_genome_effective_dict.pkl 
```

## 7. Then, the marker genes will be collocated into marker sets, following settings in original checkm. 

``Input data``:
 - ``taxon list``: a table showing the number of target lineages for calculating marker gene set list.
 - ``input tree``: the output from step 6. 
 - ``endosymbiosis dict``: the dictionary pairing the lineages with endosymbiosis event. 
 - ``discarded``: the file containing species with errornous tree placement, may be mislabelled samples, thus they are ignored from the analyses.

``Output data``:
 - ``taxon_marker_set.tsv``: a tsv file working as the reference table for binny to estimate the MAG quality based on single-copy marker genes. 

``note``: we chose the taxa names by checking their monophyly and node bootstrapping support, taxa with their node satisfying both and contains big enough number of leaves shall be included. Due to a relatively smaller number of genomes included compared with CheckM's tree, we only calculated marker gene sets for several high-rank groups.    

```sh
python $PYTHON_SCRIPT_DIR/taxon_annotation_ms.py --taxon-list "marker_gene_database/required_files_for_generating_marker_gene_database/taxon_list.txt" --input-tree "reproducibility/marker_gene_database/marker_sets_verify.tree" --endosymbiosis-map "marker_gene_database/required_files_for_generating_marker_gene_database/Endosymbiosis_dict.pkl" --output "reproducibility/marker_gene_database/taxon_marker_sets.tsv" --discarded marker_gene_database/required_files_for_generating_marker_gene_database/20240829113905_discarded_1.txt --genbank-dir marker_gene_database/genbanks_dir/ 
```

The marker sets we have now should contain a universal marker gene set to be completely functional. So we calculate the marker set containing single-copy marker genes found in > 97% of genomes and colocalize them into marker sets.

```sh
python marker_gene_database/python_scripts_for_generating_marker_sets/taxon_annotation_ms.py --universal-marker-set True --input-tree reproducibility/marker_gene_database/marker_sets_verify.tree --universal-marker-set-out "reproducibility/marker_gene_database/taxon_marker_sets_universal.tsv" --genbank-dir marker_gene_database/genbanks_dir/ --species-genome-dict marker_gene_database/required_files_for_generating_marker_gene_database/species_genome_effective_dict.pkl
```

The marker set is then concatenated at the top of the other marker sets, via a script "postprocess_ms.py"

```sh
python marker_gene_database/python_scripts_for_generating_marker_sets/postprocess_ms.py reproducibility/marker_gene_database/taxon_marker_sets_universal.tsv reproducibility/marker_gene_database/taxon_marker_sets.tsv reproducibility/taxon_marker_sets_lineage_sorted.tsv
```

We also included a list of 34 genes with each of them regarded as single marker set, proposed by Janouškovec et al. (2010) into our resulting taxon marker sets file. When we run binny, this marker set can be compared with our data in estimating completeness and purity. 

```sh
phylum	plastidcore	algae;plastidcore	1	34	34	[set(['acsF']), set(['atpA']), set(['atpB']), set(['atpH']), set(['atpI']), set(['clpC']), set(['petB']), set(['petD']), set(['petG']), set(['petN']), set(['psaA']), set(['psaB']), set(['psaC']), set(['psaD']), set(['psbA']), set(['psbB']), set(['psbC']), set(['psbD']), set(['psbE']), set(['psbF']), set(['psbH']), set(['psbI']), set(['psbJ']), set(['psbK']), set(['psbN']), set(['psbT']), set(['rpl11']), set(['rpl14']), set(['rpl16']), set(['rps12']), set(['rps19']), set(['rps5']), set(['sufB']), set(['tufA'])]
```
Simply copy and paste this line to the taxon_marker_sets_lineage_sorted.tsv, just below the line of algae (domain ms) marker set list.
Now the taxon_marker_sets.tsv file is finished, the other component is to create profile hmms. 

## 8. Align included genes.

After identifying marker genes that will be used in the database, we align the marker genes identified in the orthogroups.
``Note``: The orthogroups output from orthofinder is not named but numbered, so before this process we have to find the name of the genes. 
Unlike Pfams or TigrFams which have been summarised quite long ago, plastid marker genes hasn't been treated with equivalent efforts. Thus, we used the Orthogroups.tsv file generated from orthofinder, to identify the representative name of each orthogroup by finding the annotated gene name taking the most of the orthogroup with standard nomenclature, thus naming the group. To simplify, we offer the alignments generated by MAFFT of the identified genes directly in ChloroScan_reproducibility/marker_gene_database/marker_gene_alignments, these genes are present in the taxon_marker_sets.tsv present above. 

## 9. Build profile HMM.
Then, we create the hmm file for these genes using hmmer3.

``Input data``: 
 - alignment directory containing all genes present in marker sets.

``Output data``: 
 - directory of profile hidden markov models for each gene.

``Package``:
 - hmmer v3.4.

```sh
mkdir -p reproducibility/hmms
bash $SCRIPT_DIR/hmmbuild.sh --aln_dir marker_gene_database/marker_gene_alignments --hmm_dir reproducibility/hmms
```

## 10. Press HMMS and create database.

Binny's database requires the HMMs to be concatenated into a large HMM file, pressed and indexed for annotations. Thus, we firstly concatenate the hmms and press them using hmmpress command:

``Package``:
 - hmmer v3.4.

```sh
cat ChloroScan_reproducibility/marker_gene_database/new_hmms/* >> checkm_filtered_pf.hmm
hmmpress checkm_filtered_pf.hmm
```

## Results

Finally, create a directory called database with such structure:

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

we offer an example in ``ChloroScan_reproducibility/marker_gene_database/marker_gene_database``.

Here, ``taxon_marker_sets_lineage_sorted.tsv`` is the tsv file we just created that provides all marker sets available to use, and the ``hmms`` directory stores combined hmms in ``checkm_filtered_pf_chunk_0.hmm``. Hmms are then separated into chunks required by binny. Because the total number of HMMs is not impacting the performance of annotation, we did not split more chunks from it. We retained the names of those files from the original database used by the original version of binny. The folder ``pfam`` contains the mapping information of PFAM families to TIGRFAM families from the original CheckM database, it is no longer used since we only have plastid proteins in hmm files. But to maintain the original CheckM directory structure we did not remove this file. 

Then we moved this to binny and implemented it to ChloroScan. It is automatically loaded while running ChloroScan conda environments building.