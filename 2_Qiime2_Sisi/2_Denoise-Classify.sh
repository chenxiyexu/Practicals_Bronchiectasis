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
#This script/file covers the first 3 steps in the microbiome analysis
# using Qiime2 version 2026.7.0:
#    1.   Import data -> demux.qza
#    2.   Import metadata (.tsv) -> metadata.qzv
#    3.   Trim primers/adapters using cutadapt -> trimmed_demux.qza
#   *4.   Base quality assessment -> trimmed_demux.qzv (manual review in Qiime2 View)
#   *5.   Quality filtering/denoising using DADA2 -> AND Denoising_stats.tsv
#   *6.   Train Naive Bayes classifier for taxonomic assignation -> silva_144_v3v5_classifier.qza
#   *7.   ASV taxonomy classification using DADA2 -> Taxonomy.qza
#    8.   Rarefaction analysis.
#    8a.  Finding maximum depth.
#    8b.  ALPHA RAREFACTION (observed_features + shannon)
#    8c.  Build the MASTER per-sample table (reads at every stage + metadata)
#    8d.  Static plots + candidate depths (using R)
#    9.   Choosing sampling depth (manually)
#    9a.  Rarefy the feature table
#    9b.  Collapse rarefied table at every taxonomic level and export
#    9c.  Export raw/uncollapsed ASV level rarefied table + per-ASV taxonomy needed to analysis in R "asv" mode
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
# CLASSIFIER="${SCRIPT_DIR}/Classifier/silva_138-2_v3v4_classifier.qza"





#########################################################
### STEP 4: Base quality assessment in Qiime2 View
#########################################################
log "STEP 4: Base quality assessment in Qiime2 View and choosing which bases to remove."
#Manually inspect and choose which positions to remove
#TRUNC_BEFORE = before x position, remove all bases.
#TRUNC_AFTER = after x position, remove all bases.

echo "Please visualise the results of demux_trimmed.qzv at https://view.qiime2.org"
echo ""

echo "For the FORWARD reads."
echo "Please choose which positions to remove from the FORWARD reads."
echo "Before which position number should we discard the reads? \
(Enter 0 if starting quality is good or position N to discard bases 1 to N!)"
read -r TRUNC_BEFORE_F

echo "After which position number should we discard the FORWARD reads? \
(position chosen NOT included! For example if N = 240 we keep \
bases 1 to 240 and discard 241 to the end of the reads!)"
read -r TRUNC_AFTER_F

echo "For the REVERSE reads."
echo "Please choose which positions to remove from the REVERSE reads."
echo "Before which position number should we discard the reads? \
(Enter 0 if starting quality is good or position N to discard bases 1 to N!)"
read -r TRUNC_BEFORE_R

echo "After which position number should we discard the REVERSE reads? \
(position chosen NOT included! For example if N = 240 we keep \
bases 1 to 240 and discard 241 to the end of the reads!)"
read -r TRUNC_AFTER_R


##For the FORWARD reads.
#Checking if both variables are positive integers (0>=)
if ! [[ "${TRUNC_BEFORE_F}" =~ ^[0-9]+$ ]] || ! [[ "${TRUNC_AFTER_F}" =~ ^[0-9]+$ ]]; then
  echo "ERROR! TRUNC_BEFORE: ${TRUNC_BEFORE_F} and TRUNC_AFTER: ${TRUNC_AFTER_F} must be positive integers!"
  echo "Exiting the program due to ERROR!"
  exit 1111
fi


#Checking if TRUNC_BEFORE is smaller than TRUNC_AFTER
if [ "${TRUNC_BEFORE_F}" -ge "${TRUNC_AFTER_F}" ]; then
    echo "ERROR! TRUNC_BEFORE: ${TRUNC_BEFORE_F} must be a smaller number than TRUNC_AFTER: ${TRUNC_AFTER_F}"
    echo "Exiting the program due to ERROR!"
    exit 1111
fi

##For the REVERSE reads.
#Checking if both variables are positive integers (0>=)
if ! [[ "${TRUNC_BEFORE_R}" =~ ^[0-9]+$ ]] || ! [[ "${TRUNC_AFTER_R}" =~ ^[0-9]+$ ]]; then
  echo "ERROR! TRUNC_BEFORE: ${TRUNC_BEFORE_R} and TRUNC_AFTER: ${TRUNC_AFTER_R} must be positive integers!"
  echo "Exiting the program due to ERROR!"
  exit 1111
fi

#Checking if TRUNC_BEFORE is smaller than TRUNC_AFTER
if [ "${TRUNC_BEFORE_R}" -ge "${TRUNC_AFTER_R}" ]; then
    echo "ERROR! TRUNC_BEFORE: ${TRUNC_BEFORE_R} must be a smaller number than TRUNC_AFTER: ${TRUNC_AFTER_R}"
    echo "Exiting the program due to ERROR!"
    exit 1111
fi


log "FINISHED STEP 4: Base quality assessment in Qiime2 View."
log "For the FORWARD reads ====================================================="
log "Chosen length to discard BEFORE AND INCLUDING THIS POSITION: ${TRUNC_BEFORE_F}"
log "Chosen length to discard AFTER NOT INCLUDING THIS POSITION: ${TRUNC_AFTER_F}"

log ""

log "For the REVERSE reads ====================================================="
log "Chosen length to discard BEFORE AND INCLUDING THIS POSITION: ${TRUNC_BEFORE_R}"
log "Chosen length to discard AFTER NOT INCLUDING THIS POSITION: ${TRUNC_AFTER_R}"


#########################################################
### STEP 5: Quality filtering/denoising using DADA2
#########################################################
log "STEP 5: Quality filtering/denoising using DADA2 (paired-end)."

#--p-trim-left "${TRUNC_BEFORE}" ##Remove bases from the start of the read up to TRUNC_BEFORE
#--p-trunc-len "${TRUNC_AFTER}"  ##Remove bases after TRUCN_AFTER position until the end of the read

#TODO CHECK PSEUDO POOLING
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${QZA_DIR}/demux_trimmed.qza" \
  --p-trim-left-f "${TRUNC_BEFORE_F}" \
  --p-trim-left-r "${TRUNC_BEFORE_R}" \
  --p-trunc-len-f "${TRUNC_AFTER_F}" \
  --p-trunc-len-r "${TRUNC_AFTER_R}" \
  --p-n-threads "${THREADS}" \
  --p-pooling-method "pseudo" \
  --o-representative-sequences "${QZA_DIR}/ASV_seqs.qza" \
  --o-table "${QZA_DIR}/ASV_table.qza" \
  --o-denoising-stats "${QZA_DIR}/dada2_stats.qza" \
  --o-base-transition-stats "${QZA_DIR}/dada2_base_transition_stats.qza" \
  --verbose 2>&1 | tee -a "${LOGFILE}"

check_file "${QZA_DIR}/ASV_seqs.qza" "DADA2 representative sequences (ASVs)"
check_file "${QZA_DIR}/ASV_table.qza" "DADA2 feature table"
check_file "${QZA_DIR}/dada2_stats.qza" "DADA2 denoising stats"
check_file "${QZA_DIR}/dada2_base_transition_stats.qza" "DADA2 base transition stats"

#Transform into a .qzv to view at https://view.qiime2.org
qiime metadata tabulate \
  --m-input-file "${QZA_DIR}/dada2_stats.qza" \
  --o-visualization "${QZV_DIR}/dada2_stats.qzv"

qiime feature-table summarize \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --o-feature-frequencies "${QZA_DIR}/Feature_frequencies.qza" \
  --o-sample-frequencies "${QZA_DIR}/Sample_frequencies.qza" \
  --o-summary "${QZV_DIR}/ASV_table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "${QZA_DIR}/ASV_seqs.qza" \
  --o-visualization "${QZV_DIR}/ASV_seqs.qzv"

#Exporting the denoising stats (input/filtered/denoised/non-chimeric) as a table.
log "Exporting filtering steps stats to Denoising_stats.tsv."

#We create the folder which will contain the final results. 
mkdir -p "Final_results"
FINAL="${SCRIPT_DIR}/Final_results"

#This is used because Qiime2 automatically creates a folder with a file 
# named stats.tsv inside, for more control we use their folder and then
# delete it to be able to change the name of the final file for clarity. 
DADA2_EXPORT_DIR="$(mktemp -d)" 

qiime tools export \
  --input-path "${QZA_DIR}/dada2_stats.qza" \
  --output-path "${DADA2_EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

cp "${DADA2_EXPORT_DIR}/stats.tsv" "${FINAL}/Denoising_stats.tsv"
rm -rf "${DADA2_EXPORT_DIR}"

check_file "${FINAL}/Denoising_stats.tsv" "Exported DADA2 Denoising_stats.tsv"

log "FINISHED STEP 5: Quality filtering/denoising using DADA2."
log "Check in ${FINAL}/Denoising_stats.tsv for reads lost per sample/step."



#########################################################
### STEP 6: Train Naive Bayes classifier for taxonomic 
###         assignation 
### The classifier only has to be trained once!!!!
#########################################################
log "STEP 6: Training Naive Bayes taxonomic classifier."

echo "Please enter the path to the directory with the trained classifier (should be a .qza)!"
echo "Or directly press ENTER if it does not exist!"
read -r INPUT_TAX_DB_DIR

#If the classifier alreaady exists we just copy it to a custom directory for 
# clarity and to avoid overwriting it.
if [[ -n "${INPUT_TAX_DB_DIR}" ]]; then
  #Checking if the classifier exists
  if [[ ! -f "${INPUT_TAX_DB_DIR}" ]]; then
    echo "ERROR: ${INPUT_TAX_DB_DIR} classifier file not found!"
    exit 1111
  fi

  mkdir -p "./Classifier"
  cp "${INPUT_TAX_DB_DIR}" "./Classifier/Pre-trained_classifier.qza"

  CLASSIFIER="${PWD}/Classifier/Pre-trained_classifier.qza"

  log "Using existing trained classifier: ${INPUT_TAX_DB_DIR}, copied to ./Classifier/Pre-trained_classifier.qza"

else
  #If the classifier does NOT exist we train it using the SILVA database 138-2 release for SSU!
  mkdir -p "./Classifier"
  cd "./Classifier" || exit 1

  log "Training Naive Bayes classifier using SILVA database 138-2 release for SSU via RESCRIPt."

  ## 1) Downloading the SILVA database
  log "Downloading SILVA SSU 138.2 raw sequences and taxonomy."
  qiime rescript get-silva-data \
    --p-version "138.2" \
    --p-target "SSURef_NR99" \
    --p-include-species-labels \
    --o-silva-sequences silva_ssu_seqs_138-2.qza \
    --o-silva-taxonomy silva_ssu_tax_138-2.qza \
    2>&1 | tee -a "${LOGFILE}"

  check_file "silva_ssu_seqs_138-2.qza" "SILVA 138-2 SSU database raw sequences."
  check_file "silva_ssu_tax_138-2.qza" "SILVA 138-2 SSU database raw taxonomy."

  ## 2) Convert raw RNA (U) sequences to DNA (T)
  log "Converting SILVA 138-2 RNA (U) sequences to DNA (T)."
  qiime rescript reverse-transcribe \
    --i-rna-sequences silva_ssu_seqs_138-2.qza \
    --o-dna-sequences silva_ssu_seqs_138-2_DNA.qza \
    2>&1 | tee -a "${LOGFILE}"
    
  check_file "silva_ssu_seqs_138-2_DNA.qza" "SILVA 138-2 SSU database DNA raw sequences."

  ## 3) Remove all low quality sequences
  log "Removing low quality sequences!"
  qiime rescript cull-seqs \
    --i-sequences silva_ssu_seqs_138-2_DNA.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 8 \
    --p-n-jobs "${THREADS}" \
    --o-clean-sequences silva_138-2_seqs_cleaned.qza \
    2>&1 | tee -a "${LOGFILE}"
  
  check_file "silva_138-2_seqs_cleaned.qza" "Culled SILVA sequences of low quality."

  ## 4) Filtering sequences that are abnormally short for their domain
  #CAREFUL! This step is very heavy for the RAM!
  log "Filtering sequences by taxon-specific length."
  qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva_138-2_seqs_cleaned.qza \
    --i-taxonomy silva_ssu_tax_138-2.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva_138-2_cleaned_filtered.qza \
    --o-discarded-seqs silva_138-2_discarded.qza \
    2>&1 | tee -a "${LOGFILE}"
  
  check_file "silva_138-2_cleaned_filtered.qza" "Length-filtered SILVA sequences."


  ## 5) Extracting V3-V4 regions with in-silico PCR.

  #The primers used here must match the ones used in your own experiments!!!!!
  log "Extracting V3-V4 region with in-silico PCR."

  #Minimum maybe should be 370 or something around 390 (real). 
  #Here 2>&1 | tee -a "${LOGFILE}" for some reason causes an error
  # so its best to leave it out. 

  #NOTE: I ran it once without taking off the overhangs and somehow got Phylum Planctomycetota as a result
  # be careful in this step to only leave THE PRIMER SEQUENCES!!!!
  # --p-f-primer "FILL WITH RELEVANT PRIMERS BUT NO OVERHANGS!!!!" \       OVERHANG NOT INCLUDED: 
  # --p-r-primer "FILL WITH RELEVANT PRIMERS BUT NO OVERHANGS!!!!" \       OVERHANG NOT INCLUDED: 
  
  qiime feature-classifier extract-reads \
    --i-sequences silva_138-2_cleaned_filtered.qza \
    --p-f-primer "FILL WITH RELEVANT PRIMERS BUT NO OVERHANGS!!!!" \
    --p-r-primer "FILL WITH RELEVANT PRIMERS BUT NO OVERHANGS!!!!" \
    --p-min-length 350 \
    --p-max-length 600 \
    --p-n-jobs "${THREADS}" \
    --o-reads silva_138-2_v3v4_seqs.qza \
    --o-read-extraction-stats "${QZA_DIR}/Read_v3v4_extraction_stats.qzv"
  
  check_file "silva_138-2_v3v4_seqs.qza" "Extracted V3-v4 reference sequences."

  ## 6) Dereplicating sequences and taxonomy
  log "Dereplicating V3-V4 sequences and taxonomy."
  qiime rescript dereplicate \
    --i-sequences silva_138-2_v3v4_seqs.qza \
    --i-taxa silva_ssu_tax_138-2.qza \
    --p-mode "lca" \
    --p-threads "${THREADS}" \
    --o-dereplicated-sequences silva_138-2_v3v4_dereplicated.qza \
    --o-dereplicated-taxa silva_138-2_v3v4_tax_derep.qza \
    2>&1 | tee -a "${LOGFILE}"

  check_file "silva_138-2_v3v4_dereplicated.qza" "Dereplicated V3-V4 sequences."
  check_file "silva_138-2_v3v4_tax_derep.qza" "Dereplicated V3-V4 taxonomy."

  ## 7) Training the Naive Bayes classifier!
  #CAREFUL! This step is very heavy for the RAM!
  echo "We will proceed to train the classifier!"
  echo "Bear in mind this step requires a lot of computing power so it might take from a few minutes up to hours!"

  #The scikit-learn version must be 1.7.1!!
  log "Training the Naive Bayes classifier."
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva_138-2_v3v4_dereplicated.qza \
    --i-reference-taxonomy silva_138-2_v3v4_tax_derep.qza \
    --o-classifier "silva_138-2_v3v4_classifier.qza" \
    --verbose \
    2>&1 | tee -a "${LOGFILE}"

  check_file "silva_138-2_v3v4_classifier.qza" "Trained V3-v4 Naive Bayes classifier!"

  cd .. || exit 1

  CLASSIFIER="${SCRIPT_DIR}/Classifier/silva_138-2_v3v4_classifier.qza"

  log "FINISHED STEP 6: Newly trained classifier saved permanently at ${CLASSIFIER}."

fi


#########################################################
### STEP 7: ASV taxonomy classification using DADA2  
#########################################################
log "STEP 7: ASV taxonomy classification using DADA2."


qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads "${QZA_DIR}/ASV_seqs.qza" \
  --p-n-jobs "${THREADS}" \
  --o-classification "${QZA_DIR}/Taxonomy.qza" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${QZA_DIR}/Taxonomy.qza" "Taxonomy assignments."

qiime metadata tabulate \
  --m-input-file "${QZA_DIR}/Taxonomy.qza" \
  --o-visualization "${QZV_DIR}/Taxonomy.qzv"

check_file "${QZV_DIR}/Taxonomy.qzv" "Taxonomy table visualization."

#Create an interactive barplot of the taxonomic assignments per sample.
qiime taxa barplot \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --i-taxonomy "${QZA_DIR}/Taxonomy.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${QZV_DIR}/Taxa_barplot.qzv"

check_file "${QZV_DIR}/Taxa_barplot.qzv" "Taxonomy barplot."

#Summarizing the number of reads assigned to each domain 
# (Bacteria, Archaea, Eukaryota, Unassigned) per sample.
log "Building per-sample domain-level assignment summary \
(how many reads per sample landed in Bacteria/Archaea/Eukaryota vs. \
Unassigned)."

DOMAIN_QZA="${QZA_DIR}/Table-level1-domain.qza"
qiime taxa collapse \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --i-taxonomy "${QZA_DIR}/Taxonomy.qza" \
  --p-level 1 \
  --o-collapsed-table "${DOMAIN_QZA}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${DOMAIN_QZA}" "Domain-level (level 1) collapsed table"

#We create a temporary directory to export the domain-level 
# assignment table as a .tsv file.
DOMAIN_EXPORT_DIR="$(mktemp -d)"
qiime tools export \
  --input-path "${DOMAIN_QZA}" \
  --output-path "${DOMAIN_EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

#biom is a package that handles BIOM format files! Qiime2 uses it to 
# represent biological samples by observation contingency tables.

#Domain_assignment shows a row name called #OTU ID but it is 
# actually our sample´s id. 
biom convert \
  -i "${DOMAIN_EXPORT_DIR}/feature-table.biom" \
  -o "${FINAL}/Domain_assignment.tsv" --to-tsv \
  2>&1 | tee -a "${LOGFILE}"

rm -rf "${DOMAIN_EXPORT_DIR}"

check_file "${FINAL}/Domain_assignment.tsv" "Per-sample domain-level assignment table"

log "FINISHED STEP 7: Finished ASV taxonomy classification using DADA2. \
Check ${FINAL}/Domain_assignment.tsv for per-sample domain-level \
assignment table."



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 2_Denoise-Classify.sh finished running!                                                             "
log "                                                                                                           "
log "DENOISING AND TAXONOMIC CLASSIFICATION COMPLETE. Key outputs:                                              "
log "  Artifacts (.qza):                                                                                        "
log "    - ${QZA_DIR}/ASV_table.qza                       (ASV feature count table)                             "
log "    - ${QZA_DIR}/ASV_seqs.qza                        (ASV representative sequences)                        "
log "    - ${QZA_DIR}/Taxonomy.qza                        (Taxonomic assignments)                               "
log "  Visualizations (.qzv - view at https://view.qiime2.org):                                                 "
log "    - ${QZV_DIR}/dada2_stats.qzv                     (DADA2 filtering stats per sample)                    "
log "    - ${QZV_DIR}/ASV_table.qzv                       (ASV count distribution per sample)                   "
log "    - ${QZV_DIR}/Taxa_barplot.qzv                    (Interactive taxonomic composition barplot)           "
log "  Exported Summaries:                                                                                      "
log "    - ${FINAL}/Denoising_stats.tsv                   (Reads retained after filtering and chimera removal)  "
log "    - ${FINAL}/Domain_assignment.tsv                 (Read counts grouped by Bacteria/Archaea/Eukaryota)   "
log "==========================================================================================================="
log "NEXT STEPS:                                                                                                "
log "  1. Inspect ${QZV_DIR}/Taxa_barplot.qzv at https://view.qiime2.org to evaluate taxonomy distributions.    "
log "  2. Run Python merge script to unite metadata, denoising stats, and domain tables into one master summary."
log "  3. Proceed to Rarefaction Analysis (alpha/beta diversity steps).                                         "
log "==========================================================================================================="