# Balthazar

## 1. Generate all possible SNVs of the GOIs


*clinvar_result_FKRP_04202023.txt* and *clinvar_result_04202023_LARGE1.txt* was manually downloaded from ClinVar webpage.

### for FKRP

python 22.1.25\ mutation\ classfication\ v3.py fkrp_1to1485.fasta all_possible_SNVs
python add_ClinVar_to_allpossibleSNVs_FKRP.py all_possible_SNVs clinvar_result_FKRP_04202023.txt fkrp_1to1485.fasta all_possible_SNVs_with_clinvar_FKRP_04202023
Rscript --vanilla all_possible_SNVs_plot.r all_possible_SNVs

### for LARGE1

python 22.1.25\ mutation\ classfication\ v3.py large1_1to2268.fasta all_possible_SNVs
python add_ClinVar_to_allpossibleSNVs_LARGE1.py all_possible_SNVs clinvar_result_04202023_LARGE1.txt large1_1to2268.fasta all_possible_SNVs_with_clinvar_LARGE1_04202023
Rscript --vanilla all_possible_SNVs_plot.r all_possible_SNVs


## 2. Cloning alert

This package identifies two groups of variants that might be underrepresented in the pool due to the limitation of the cloning method.

*group1* will introduce a cut before the variants which will remove the variant from the top strand
*group2* will introduce a cut right after the variants which might affect the elongation of the top strand

python cloning_alert_FKRP.py fkrp_1to1485.fasta all_possible_SNVs_FKRP out_FKRP
python cloning_alert_LARGE1.py large1_1to2268.fasta all_possible_SNVs_LARGE1 out_LARGE1

## 3. generate all oligos to synthesize

**important notes** are in the python scripts

python PALS_C_oligos_FKRP.py FKRP_oligo_to_order_out.fa
python PALS_C_oligos_LARGE1.py LARGE1_oligo_to_order_out.fa

## 4. plasmid pool QC

The jobs were run on Yale HPC Ruddle (RIP)

Sbatch QC_FKRP.slurm
Sbatch QC_LARGE1.slurm
