#!/bin/bash

# base locale
homeLOC="/home/aubonnier"
# base GO - fonctionne
homeGO="/home/genouest/cnrs_umr6553/aubonnier"
# scripts GO - fonctionne
scriptab="$homeGO/scriptab"
# dossier principal des exports EXP - fonctionne
exp="$homeGO/EXP"
# sous-dossiers - fonctionne
output="$exp/output"
log="$homeGO/genouest_log"

#arbre d'especes - fonctionne
#SPtree="/projects/synomics/DeCoSTAR/results/speciesTree.newick" ce n'est pas un dossier
# Nom fichier => SpeciesTree.phyloxml
SPtree2="/projects/synomics/DeCoSTAR/results"

#Variable projet 1) reconciled_gene_trees - fonctionne
RECNE="/projects/synomics/data_sharing_Audrey/reconciled_gene_trees/newick_trees"
RECNHX="/projects/synomics/data_sharing_Audrey/reconciled_gene_trees/nhx_trees"
RECXML="/projects/synomics/data_sharing_Audrey/reconciled_gene_trees/xml_trees"
#Variable projet 2) simple_gene_tree - fonctionne
SIMG="/projects/synomics/data_sharing_Audrey/simple_gene_trees"
#Stat synomics - fonctionne
SYN="/projects/synomics"


### General MXE ###
#plantlp="/projects/plantlp"
#snp_calling="$plantlp/01_SNP_CALLING"
#vcf_processing="$plantlp/02_VCF_PROCESSING"
#ne="$plantlp/03_NE_INFERENCE"
#softwares="$plantlp/softwares"
#scratch="/scratch/mbrault"

### Ne analysis ###
#ne_scripts="$ne/scripts"
#ne_data="$ne/data"
#ne_results="$ne/results"  

### SNP calling ###
#snp_scripts="$snp_calling/scripts"
#snp_data="$snp_calling/data"
#snp_results="$snp_calling/results"

### VCF processing ###
#vcf_scripts="$vcf_processing/scripts"
#vcf_data="$vcf_processing/data"
#vcf_results="$vcf_processing/results"
#vcf_log="$vcf_processing/genouest_log"
