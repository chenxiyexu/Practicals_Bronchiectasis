#!/usr/bin/env bash
set -euo pipefail
#
#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 
#
#Special thanks to Pneumobiome from IRYCIS for their support and for allowing me to 
# research such an interesting topic, and above all special thanks to Iván Linares Ambohades for 
# letting me use his scripts as a base and for his help in the development of this pipeline.
#
#The full analysis requires 12 steps, which comprise of 7 files/scripts. The rest of the analysis 
# continues in R. 
#
#This script/file covers the first 3 steps in the microbiome analysis
# using Qiime2 version 2026.7.0:
#    1.   Import data -> demux.qza
#    2.   Import metadata (.tsv) -> metadata.qzv
#    3.   Trim primers/adapters using Cutadapt -> trimmed_demux.qza
#    4.   Base quality assessment -> trimmed_demux.qzv (manual review in Qiime2 View)
#    5.   Quality filtering/denoising using DADA2 -> ASV_table.qza AND Denoising_stats.tsv
#    6.   Train Naive Bayes classifier for taxonomic assignation -> SILVA or GREENGENES classifier.qza (or Pre-trained)
#    7.   ASV taxonomy classification using classifier -> Taxonomy.qza
#    8.   Rarefaction analysis.
#    8a.  Finding maximum depth.
#    8b.  ALPHA RAREFACTION (observed_features + shannon).
#    8c.  Build the MASTER per-sample table (reads at every stage + metadata).
#    8d.  Static plots + candidate depths (using R).
#    9.   Choosing sampling depth (manually).
#    9a.  Rarefy the feature table.
#    9b.  Collapse rarefied table at every taxonomic level and export.
#    10.  Export ASV and genus level rarefied and raw tables.
#    10a. Export non-rarefied ASV table.
#    10b. Collapse and export non-rarefied genus table
#    10c. Export rarefied ASV table
#    10d. Collapse and export rarefied genus table
#    10e. Export ASV taxonomy mapping
#   *11.  Building phylogenetic tree from ASV representative sequences.
#    12. Core metrics (only at ASV-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#    13. Core metrics (at genus-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#
#Resulting .qzv files must be viewed using https://view.qiime2.org/, in this
# case for manual quality assessment to decide which bases/positions to 
# remove from the reads. 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Global_functions.sh"
trap "error_check ${LINENO}" ERR

########################################################
## Computing power assigned to the script
########################################################
#Change this depending on which PC you use 
# (nproc --all #this shows total amount of cores/CPUs)
THREADS=8

########################################################
## Variables needed!!!
########################################################
QZA_DIR="${SCRIPT_DIR}/Output_qza"
QZV_DIR="${SCRIPT_DIR}/Output_qzv"
METADATA="${SCRIPT_DIR}/metadata_Qiime_adapted.tsv"
FINAL="${SCRIPT_DIR}/Final_results"



#########################################################
### STEP 11: Building phylogenetic tree from ASV representative sequences.
#########################################################
#This is REQUIRED for two things later in the R analysis:
#  - Faith's Phylogenetic Diversity (alpha diversity)
#  - Weighted UniFrac (beta diversity)
#
#Neither of those metrics can be computed after collapsing to genus/family/
# etc: collapsing merges many ASVs (each its own tip on the tree) into one
# row, so there's no longer a 1-to-1 match between abundance-table rows and
# tree tips. Because of this, Faith's PD and Weighted UniFrac are ALWAYS
# computed at the ASV level in 06_analysis_16S.R, regardless of whether
# you have TAXA_LEVEL set to "genus" or "asv" for everything else.
#
#Run this any time after 02_denoise_and_classify.sh (it only needs
# ASV_seqs.qza, so it does not need to be re-run if you change rarefaction
# depth later -- the tree only depends on which ASVs exist, not on rarefied
# counts).
#########################################################
log "STEP 11: Building phylogenetic tree from ASV representative sequences."
check_file "${QZA_DIR}/ASV_seqs.qza" "ASV representative sequences (from script 2_Denoise-Classify.sh)"

#Create a final folder for the tree .tsv outputs. 
TREE_TSV_DIR="${FINAL}/Tree"
mkdir -p "${TREE_TSV_DIR}"

#Create a folder for tree related .qza outputs.
TREE_QZA_DIR="${QZA_DIR}/Tree"
mkdir -p "${TREE_QZA_DIR}"

log "Aligning ASVs (mafft) -> masking -> building tree (fasttree) -> rooting..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QZA_DIR}/ASV_seqs.qza" \
  --p-n-threads "${THREADS}" \
  --o-alignment "${TREE_QZA_DIR}/Aligned_ASV_seqs.qza" \
  --o-masked-alignment "${TREE_QZA_DIR}/Masked_aligned_ASV_seqs.qza" \
  --o-tree "${TREE_QZA_DIR}/Unrooted_tree.qza" \
  --o-rooted-tree "${TREE_QZA_DIR}/Rooted_tree.qza" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${TREE_QZA_DIR}/Rooted_tree.qza" "Rooted phylogenetic tree"

qiime empress tree-plot \
  --i-tree "${TREE_QZA_DIR}/Rooted_tree.qza" \
  --o-visualization "${QZV_DIR}/tree.qzv" \
  2>&1 | tee -a "${LOGFILE}" || \
  log "NOTE: 'qiime empress' plugin not installed -- skipping optional tree.qzv view (harmless, the Newick export below is what matters)."

log "Exporting tree as plain Newick."

#Temporary folder (the qiime2 function exports to a folder, not a single file, so we need to rename it).
export_dir="$(mktemp -d)"

qiime tools export \
  --input-path "${TREE_QZA_DIR}/Rooted_tree.qza" \
  --output-path "${export_dir}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${export_dir}/tree.nwk" "Exported Newick tree"

#We change the name here because it is actually a rooted tree. 
#On a personal note, do NOT change tree.nwk to Tree.nwk, the qiime2 command creates a folder 
# with a FIXED NAME!!!! You like names starting with CAPITAL LETTERS but be careful!!!!!
cp "${export_dir}/tree.nwk" "${TREE_TSV_DIR}/Rooted_tree.nwk"

rm -rf "${export_dir}"

check_file "${TREE_TSV_DIR}/Rooted_tree.nwk" "Newick tree copied to tablas/tree/"

log "FINISHED STEP 11: Building phylogenetic tree from ASV representative sequences."



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 6_Build-tree.sh finished running!                                                                   "
log "                                                                                                           "
log "TREE BUILDING COMPLETE. Key outputs:                                                                       "
log "  Artifacts (.qza):                                                                                        "
log "    - ${TREE_QZA_DIR}/Aligned_ASV_seqs.qza                (Aligned ASV sequences using Mafft)                   "
log "    - ${TREE_QZA_DIR}/Masked_aligned_ASV_seqs.qza         (Masked aligned ASV sequences)                        "
log "    - ${TREE_QZA_DIR}/Unrooted_tree.qza                   (Unrooted phylogenetic tree)                          "
log "    - ${TREE_QZA_DIR}/Rooted_tree.qza                     (Rooted phylogenetic tree)                            "
log "  Visualizations (.qzv - view at view.qiime2.org):                                                         "
log "    - ${QZV_DIR}/tree.qzv                            (Phylogenetic tree visualization)                     "  
log "  Exported Summaries:                                                                                      "
log "    - ${TREE_TSV_DIR}/Rooted_tree.nwk                (Rooted tree used in R with ape::read.tree() for Faith's PD and Weighted UniFrac)"
log "==========================================================================================================="
log "NEXT STEPS:                                                                                                "
log "  1. Use ${TREE_TSV_DIR}/Rooted_tree.nwk in R for downstream phylogenetic diversity calculations.          "
log "  2. Proceed to downstream R integration and final diversity analyses.                                     "
log "==========================================================================================================="