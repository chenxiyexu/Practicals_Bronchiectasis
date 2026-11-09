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
#   *8.   Rarefaction analysis.
#   *8a.  Finding maximum depth.
#   *8b.  ALPHA RAREFACTION (observed_features + shannon).
#   *8c.  Build the MASTER per-sample table (reads at every stage + metadata).
#   *8d.  Static plots + candidate depths (using R).
#    9.   Choosing sampling depth (manually).
#    9a.  Rarefy the feature table.
#    9b.  Collapse rarefied table at every taxonomic level and export.
#    10.  Export ASV and genus level rarefied and raw tables.
#    10a. Export non-rarefied ASV table.
#    10b. Collapse and export non-rarefied genus table
#    10c. Export rarefied ASV table
#    10d. Collapse and export rarefied genus table
#    10e. Export ASV taxonomy mapping
#    11.  Building phylogenetic tree from ASV representative sequences.
#    12. Core metrics (only at ASV-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#    13. Core metrics (at genus-level, for exploratory purposes, for a more complete analysis go to the R scripts).
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
CLASSIFIER="${SCRIPT_DIR}/Classifier/silva_138-2_v3v4_classifier.qza"


#########################################################
### STEP 8: Rarefaction analysis
#########################################################
log "STEP 8: Rarefaction analysis."
#Manually inspect and choose which depth of sequencing to keep for analysis. 

check_file "${QZA_DIR}/ASV_table.qza" "ASV feature table"
check_file "${QZA_DIR}/dada2_stats.qza" "DADA2 denoising stats"
check_file "${METADATA}" "Sample metadata"

RAREDIR="${FINAL}/Rarefaction"
mkdir -p "${RAREDIR}"



#########################################################
### STEP 8a: Finding maximum depth
#########################################################
#Find the sample with the most amount of sequencing depth to set 
# the rarefaction curve´s max depth. 
log "STEP 8a: Finding maximum depth of the samples."

ASV_EXPORT_DIR="$(mktemp -d)"
qiime tools export \
  --input-path "${QZA_DIR}/ASV_table.qza" \
  --output-path "${ASV_EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${ASV_EXPORT_DIR}/feature-table.biom" "Exported ASV-level BIOM"

biom summarize-table \
  -i "${ASV_EXPORT_DIR}/feature-table.biom" \
  -o "${RAREDIR}/biom-summary.txt"
rm -rf "${ASV_EXPORT_DIR}"

check_file "${RAREDIR}/biom-summary.txt" "BIOM per-sample summary"


#awk must work with single quotes '' because if double quotes are 
# used, Bash will try to interpret $ as a variable and crash!!  
log "Writing a clean per-sample depth table (sample_depths.tsv)."
{
  echo -e "sample_id\tdepth"
  grep -A100000 "Counts/sample detail" "${RAREDIR}/biom-summary.txt" \
    | tail -n +2 \
    | awk -F: 'NF==2 {
        gsub(/^[ \t\r]+|[ \t\r]+$/, "", $1);
        val=$2;
        gsub(/^[ \t\r]+|[ \t\r]+$/, "", val);

        #Handle thousands and decimal separators
        if (val ~ /,[0-9]{3}\./) {
          # English format: 12,710.000.
          gsub(/,/, "", val);
        } else if (val ~ /,[0-9]+$/) {
          #European format: 12.345,000.
          gsub(/\./, "", val);
          gsub(/,/, ".", val);
        }

        #Remove any remaining thousands separators only.
        if (val ~ /^[0-9]+(\.[0-9]+)?$/)
          print $1 "\t" val
      }'
} > "${RAREDIR}/sample_depths.tsv"

check_file "${RAREDIR}/sample_depths.tsv" "Clean per-sample depth table"

N_DEPTHS=$(($(wc -l < "${RAREDIR}/sample_depths.tsv") - 1))
if [ "${N_DEPTHS}" -le 0 ]; then
  log "FAIL - sample_depths.tsv has no data rows. The "Counts/sample detail""
  log "FAIL - section in ${RAREDIR}/biom-summary.txt has a number format this"
  log "FAIL - parser doesn´t recognize. Check it manually with:"
  log "FAIL - grep -A20 \"Counts/sample detail\" ${RAREDIR}/biom-summary.txt"
  exit 1111
fi

MAX_DEPTH=$(tail -n +2 "${RAREDIR}/sample_depths.tsv" | cut -f2 | sort -n | tail -1 | cut -d"." -f1)

if [ -z "${MAX_DEPTH}" ] || [ "${MAX_DEPTH}" -le 0 ]; then
  log "FAIL - Could not determine a valid max depth from sample_depths.tsv."
  exit 1
fi

log "FINISHED STEP 8a: Deepest sample has ${MAX_DEPTH} reads (${N_DEPTHS} samples total); using this as max-depth."



#########################################################
### STEP 8b: ALPHA RAREFACTION (observed_features + shannon)
#########################################################
log "STEP 8b: Running qiime diversity alpha-rarefaction (max depth = ${MAX_DEPTH})."

qiime diversity alpha-rarefaction \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --p-max-depth "${MAX_DEPTH}" \
  --p-steps 20 \
  --p-metrics observed_features \
  --p-metrics shannon \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${QZV_DIR}/ALPHA_rarefaction.qzv" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${QZV_DIR}/ALPHA_rarefaction.qzv" "Alpha rarefaction visualization"

qiime tools export \
  --input-path "${QZV_DIR}/ALPHA_rarefaction.qzv" \
  --output-path "${RAREDIR}/alpha_rarefaction_export" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${RAREDIR}/alpha_rarefaction_export/observed_features.csv" "Raw observed_features rarefaction data"
check_file "${RAREDIR}/alpha_rarefaction_export/shannon.csv" "Raw shannon rarefaction data"

log "FINISHED STEP 8b: Alpha rarefaction curves can be found at ${FINAL}/Rarefaction"



#########################################################
### STEP 8c: Build the MASTER per-sample table (reads at every stage + metadata)
#########################################################
log "STEP 8c: Building MASTER per-sample table (reads + taxonomy + metadata)."

BS_ARGS=(--metadata "${METADATA}" --output "${RAREDIR}/MASTER_sample_summary.tsv")
[ -s "${FINAL}/reads_00_raw.tsv" ]        && BS_ARGS+=(--raw-counts "${FINAL}/reads_00_raw.tsv")
[ -s "${FINAL}/reads_01_trimmed.tsv" ]    && BS_ARGS+=(--trimmed-counts "${FINAL}/reads_01_trimmed.tsv")
[ -s "${FINAL}/Denoising_stats.tsv" ]         && BS_ARGS+=(--dada2-stats "${FINAL}/Denoising_stats.tsv")
[ -s "${FINAL}/Domain_assignment.tsv" ]   && BS_ARGS+=(--domain-table "${FINAL}/Domain_assignment.tsv")

python3 "./build_sample_summary.py" "${BS_ARGS[@]}" \
  2>&1 | tee -a "${LOGFILE}"
  
check_file "${RAREDIR}/MASTER_sample_summary.tsv" "Master per-sample summary table"

log "FINISHED STEP 8c: Built MASTER per-sample table (reads + taxonomy + metadata) at ${RAREDIR}/MASTER_sample_summary.tsv"



#########################################################
### STEP 8d: Static plots + candidate depths (using R)
#########################################################
log "STEP 8d: Generating static rarefaction plots and candidate depths (R)."

mkdir -p "${RAREDIR}/figures"
FIGURAS_DIR="${RAREDIR}/figures"

#To generate the graphs we use an R script. 
Rscript "${SCRIPT_DIR}/rarefaction_plots.R" "${SCRIPT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${FIGURAS_DIR}/rarefaction_observed_features.png" "Observed-features rarefaction plot"
check_file "${FIGURAS_DIR}/rarefaction_shannon.png" "Shannon rarefaction plot"
check_file "${FIGURAS_DIR}/samples_retained_vs_depth.png" "Samples-retained-vs-depth plot"
check_file "${RAREDIR}/rarefaction_depth_candidates.tsv" "Candidate rarefaction depths table"

log "FINISHED STEP 8d: Generating static rarefaction plots and candidate depths (R)."



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 3_Rarefaction-Analysis.sh finished running!                                                         "
log "                                                                                                           "
log "RAREFACTION ANALYSIS COMPLETE. Key outputs:                                                                "
log "  Visualizations (.qzv - view at view.qiime2.org):                                                         "
log "    - ${QZV_DIR}/ALPHA_rarefaction.qzv                   (Interactive alpha rarefaction plot)              "
log "  Exported Summaries & Plots:                                                                              "
log "    - ${RAREDIR}/MASTER_sample_summary.tsv               (Reads at every stage + all metadata)             "
log "    - ${RAREDIR}/rarefaction_depth_candidates.tsv        (Candidate depths: samples kept, diversity)       "
log "    - ${FIGURAS_DIR}/rarefaction_observed_features.png                                                     "
log "    - ${FIGURAS_DIR}/rarefaction_shannon.png                                                               "
log "    - ${FIGURAS_DIR}/samples_retained_vs_depth.png                                                         "
log "==========================================================================================================="
log "NEXT STEPS:                                                                                                "
log "  1. Review the rarefaction plots and candidate depths table.                                              "
log "  2. Choose a valid sampling depth to use for rarefaction.                                                 "
log "  3. Proceed to Script 4 for feature table rarefaction and taxonomic collapses.                            "
log "==========================================================================================================="