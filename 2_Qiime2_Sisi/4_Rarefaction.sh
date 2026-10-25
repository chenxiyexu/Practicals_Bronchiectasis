#!/usr/bin/env bash
set -euo pipefail
#
#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 
#
#Special thanks to TODO 
#
#The full analysis requires TODO steps, which comprise of TODO files. The rest of the analysis 
# continues in R. 
#
# This script/file covers the first 3 steps in the microbiome analysis
# using Qiime2 version 2026.7.0:
#    1.   Import data -> demux.qza
#    2.   Import metadata (.tsv) -> metadata.qzv
#    3.   Trim primers/adapters using cutadapt -> trimmed_demux.qza
#    4.   Base quality assessment -> trimmed_demux.qzv (manual review in Qiime2 View)
#    5.   Quality filtering/denoising using DADA2 -> AND Denoising_stats.tsv
#    6.   Train Naive Bayes classifier for taxonomic assignation -> silva_144_v3v5_classifier.qza
#    7.   ASV taxonomy classification using DADA2 -> Taxonomy.qza
#    8.   Rarefaction analysis.
#    8a.  Finding maximum depth.
#    8b.  ALPHA RAREFACTION (observed_features + shannon)
#    8c.  Build the MASTER per-sample table (reads at every stage + metadata)
#    8d.  Static plots + candidate depths (using R)
#   *9.   Choosing sampling depth (manually)
#   *9a.  Rarefy the feature table
#   *9b.  Collapse rarefied table at every taxonomic level and export
#    10.  Export ASV and genus level rarefied and raw tables.
#    10a. Export non-rarefied ASV table
#    10b. Collapse and export non-rarefied genus table
#    10c. Export rarefied ASV table
#    10d. Collapse and export rarefied genus table
#    10e. Export ASV taxonomy mapping
#    11.  Building phylogenetic tree from ASV representative sequences.
#    12. Core metrics (only at ASV-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#
#
#Resulting .qzv files must be viewed using https://view.qiime2.org/, in this
# case for manual quality assessment to decide which bases/positions to 
# remove from the reads. 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Global_functions.sh"
trap 'error_check "$?" "$LINENO" "$BASH_COMMAND"' ERR

########################################################
## Computing power assigned to the script
########################################################
# Change this depending on which PC you use 
# (nproc --all #this shows total amount of cores/CPUs)
THREADS=8

########################################################
## Variables needed!!!
########################################################
QZA_DIR="${SCRIPT_DIR}/Output_qza"
QZV_DIR="${SCRIPT_DIR}/Output_qzv"
METADATA="${SCRIPT_DIR}/metadata_Qiime_adapted.tsv"
FINAL="${SCRIPT_DIR}/Final_results"
CLASSIFIER="${SCRIPT_DIR}/Classifier/silva_138-2_v3v4_classifier.qza"

RAREDIR="${FINAL}/Rarefaction"
MAX_DEPTH=$(tail -n +2 "${RAREDIR}/sample_depths.tsv" | cut -f2 | sort -n | tail -1 | cut -d"." -f1)





#########################################################
### STEP 9: Choosing sampling depth (manually). 
#########################################################
log "STEP 9: Choosing sampling depth (manually). "

#Ask the user to introduce a valid SAMPLING_DEPTH (positive integer!!).
echo "Please revise rarefaction_depth_candidates.tsv and the rarefaction figures to choose the sampling depth!"
echo "Bear in mind it must be a positive interger (for example: 7000):"
read -r SAMPLING_DEPTH

#Check if SAMPLING_DEPTH is a positive integer. 
if ! [[ "${SAMPLING_DEPTH}" =~ ^[0-9]+$ ]]; then
  echo "SAMPLING_DEPTH must be a positive integer."
  exit 1111
fi

#Check if SAMPLING_DEPTH is smaller than MAX_DEPTH.
if [ "${SAMPLING_DEPTH}" -gt "${MAX_DEPTH}" ]; then
  log "FAIL - SAMPLING_DEPTH (${SAMPLING_DEPTH}) cannot exceed MAX_DEPTH (${MAX_DEPTH})."
  exit 1111
fi

log "The chosen sampling depth was: ${SAMPLING_DEPTH}."
log "FINISHED STEP 9: Choosing sampling depth (manually)."

check_file "${QZA_DIR}/ASV_table.qza" "ASV feature table (from step 2)"
check_file "${QZA_DIR}/Taxonomy.qza" "Taxonomy assignments (from step 2)"



#########################################################
### STEP 9a: Rarefy the feature table
#########################################################
log "STEP 9a: Rarefying feature table at sampling depth = ${SAMPLING_DEPTH}"

qiime feature-table rarefy \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --o-rarefied-table "${QZA_DIR}/ASV_table_rarefied.qza" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${QZA_DIR}/ASV_table_rarefied.qza" "Rarefied ASV feature table"

qiime feature-table summarize \
  --i-table "${QZA_DIR}/ASV_table_rarefied.qza" \
  --m-metadata-file "${METADATA}" \
  --o-feature-frequencies "${QZA_DIR}/Feature_frequencies_rarefied.qza" \
  --o-sample-frequencies "${QZA_DIR}/Sample_frequencies_rarefied.qza" \
  --o-summary "${QZV_DIR}/ASV_table_rarefied.qzv" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${QZV_DIR}/ASV_table_rarefied.qzv" "Rarefied feature table summary"

log "FINISHED STEP 9a: Open ${QZV_DIR}/.qzv to confirm which"
log "   samples were dropped (any sample with fewer than "
log "   ${SAMPLING_DEPTH} reads disappears from every table below)."



#########################################################
### STEP 9b: Collapse rarefied table at every taxonomic level and export
#########################################################
log "STEP 9b: Collapse rarefied table at every taxonomic level and export."

TAXA_QZA_DIR="${QZA_DIR}/taxa_tables_rarefied"
TAXA_TSV_DIR="${FINAL}/taxa_tables_rarefied"

mkdir -p "${TAXA_QZA_DIR}" "${TAXA_TSV_DIR}"

declare -A LEVEL_NAMES=(
  [1]="domain" [2]="kingdom" [3]="phylum" [4]="class"
  [5]="order" [6]="family"  [7]="genus"
)

log "Collapsing the rarefied table at every taxonomic level."
for level in 1 2 3 4 5 6 7; do
  rank_name="${LEVEL_NAMES[$level]}"
  log "  Level ${level} (${rank_name})..."

  collapsed_qza="${TAXA_QZA_DIR}/table-level${level}-${rank_name}-rarefied.qza"
  qiime taxa collapse \
    --i-table "${QZA_DIR}/ASV_table_rarefied.qza" \
    --i-taxonomy "${QZA_DIR}/Taxonomy.qza" \
    --p-level "${level}" \
    --o-collapsed-table "${collapsed_qza}" \
    2>&1 | tee -a "${LOGFILE}"

  check_file "${collapsed_qza}" "Rarefied collapsed table (level ${level}, ${rank_name})"

  export_dir="$(mktemp -d)" #Temporary folder
  qiime tools export \
    --input-path "${collapsed_qza}" \
    --output-path "${export_dir}" \
    2>&1 | tee -a "${LOGFILE}"

  #File name cannot be changed since Qiime2 here creates a folder with the file inside. 
  check_file "${export_dir}/feature-table.biom" "Exported rarefied BIOM (level ${level}, ${rank_name})"

  tsv_out="${TAXA_TSV_DIR}/table-level${level}-${rank_name}-rarefied.tsv"
  biom convert -i "${export_dir}/feature-table.biom" -o "${tsv_out}" --to-tsv \
    2>&1 | tee -a "${LOGFILE}"

  rm -rf "${export_dir}"

  check_file "${tsv_out}" "Rarefied TSV table (level ${level}, ${rank_name})"
done

log "FINISHED STEP 9b: Collapse rarefied table at every taxonomic level and export"



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 4_Rarefaction.sh finished running!                                                                  "
log "                                                                                                           "
log "RAREFACTION AND TAXONOMIC COLLAPSING COMPLETE. Key outputs:                                                "
log "  Artifacts & Visualizations (.qza / .qzv):                                                                "
log "    - ${QZA_DIR}/ASV_table_rarefied.qza             (Rarefied ASV feature count table)                     "
log "    - ${QZV_DIR}/ASV_table_rarefied.qzv             (Rarefied feature table summary)                       "
log "  Exported Summaries & Tables:                                                                             "
log "    - ${FINAL}/taxa_tables_rarefied/                (Rarefied tables collapsed across levels 1-7)          "
log "==========================================================================================================="
log "NEXT STEPS:                                                                                                "
log "  1. Proceed to Script 5 for final table organization, raw/rarefied ASV exports, and taxonomy mapping.     "
log "==========================================================================================================="