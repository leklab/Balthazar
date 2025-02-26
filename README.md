**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies. The manuscript is currently available on [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3). This workflow involves a 2-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated with short-read Next-generation sequencing. 

We have organized scripts used for SMuRF across three repositories:
* **Balthazar (this repository)** is used to create the oligo library needed for the saturation mutagenesis as well as perform QC for the resulting construct.    
* [**Gagarmel repository**](https://github.com/leklab/Gargamel) is used to process the raw data from NGS to quantify the enrichment of the variants.    
* [**Azrael repository**](https://github.com/leklab/Azrael) is used to generate and analyze the functional scores, and plot the results.   

## Requirements   
* Python (tested on v2.7.16 or v3.6)
* BioPython (tested on v1.78)   
* R (tested on v4.1.2)    
* BWA (tested on v0.7.15)    
* SAMtools    
* FASTX-Toolkit (tested on v0.0.14)    
* picard 2.9.0 (tested on v0.7.17)

### R libraries required    
```
ggplot2    
```

## 1. Generate all possible SNVs of the GOIs

### Input files

1. Fasta file containing coding sequence of gene of interest
2. Clinvar data file subsetted to gene of interest

### FKRP example
Files associated with this example are in `all_possible_SNVs/fkrp`
```
# Generate a file containing all possible SNVs
python mutation_classfication.py fkrp_1to1485.fasta all_possible_SNVs

# Add ClinVar annotation to this file
python add_ClinVar_to_allpossibleSNVs_FKRP.py all_possible_SNVs clinvar_result_FKRP_04202023.txt fkrp_1to1485.fasta 
all_possible_SNVs_with_clinvar_FKRP_04202023

# Create a barplot in R
Rscript --vanilla all_possible_SNVs_plot.r all_possible_SNVs
```
### LARGE1 example
Files associated with this example are in `all_possible_SNVs/large1`
```
# Generate a file containing all possible SNVs
python mutation_classfication.py large1_1to2268.fasta all_possible_SNVs

# Add ClinVar annotation to this file
python add_ClinVar_to_allpossibleSNVs_LARGE1.py all_possible_SNVs clinvar_result_04202023_LARGE1.txt large1_1to2268.fasta all_possible_SNVs_with_clinvar_LARGE1_04202023

# Create a barplot in R
Rscript --vanilla all_possible_SNVs_plot.r all_possible_SNVs
```

### A generalized script is also provided, which is applicable to other genes of interest.

Please find the **scripts** and the **readme.md** in the cloning alert directory

```

# all_possible_SNVs_test.txt were generated with SNV_mutation_classification.py.
# all_possible_NNKmut_test.txt were generated with NNK_mutation_classification.py.

python SNV_mutation_classification.py test_input.fa all_possible_SNVs_test.txt 7311 8796
python NNK_mutation_classification.py test_input.fa all_possible_NNKmut_test.txt 7311 8796

```


## 2. Call the variants that may be less represented in the pool

This package identifies two groups of variants that might be underrepresented in the pool due to the limitation of the cloning method.

*group1*: will introduce a cut before the variants which will remove the variant from the top strand  
*group2*: will introduce a cut right after the variants which might affect the elongation of the top strand  

**readme.md** in the cloning alert directory

## 3. Generate all oligos to synthesize

**readme.md** in the oligo generator directory

## 4. Plasmid pool QC

The jobs were run on Yale HPC Ruddle (RIP)
```
# FKRP example
Sbatch QC_FKRP.slurm

# LARGE1 example
Sbatch QC_LARGE1.slurm
```
