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
#   *4.   Base quality assessment -> trimmed_demux.qzv (manual review in Qiime2 View)
#   *5.   Quality filtering/denoising using DADA2 -> ASV_table.qza AND Denoising_stats.tsv
#   *6.   Train Naive Bayes classifier for taxonomic assignation -> SILVA or GREENGENES classifier.qza (or Pre-trained)
#   *7.   ASV taxonomy classification using classifier -> Taxonomy.qza
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

#Using pseudopooling for better rare ASV detection. 
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

mkdir -p "${SCRIPT_DIR}/Classifier"

##########################################################
### OPTION A) USING AN EXISTING PRE-TRAINED CLASSIFIER
##########################################################
#This avoids downloading reference databases and training a new model.
if [[ -n "${INPUT_TAX_DB_DIR}" ]]; then

  #Check that the supplied path points to an existing file.
  if [[ ! -f "${INPUT_TAX_DB_DIR}" ]]; then
    echo "ERROR: ${INPUT_TAX_DB_DIR} classifier file not found!"
    exit 1
  fi

  #Check the file extension.
  if [[ "${INPUT_TAX_DB_DIR}" != *.qza ]]; then
    echo "ERROR: The classifier must be a .qza file!"
    exit 1
  fi

  #Copy the selected classifier to a standard location.
  #The original file is not modified.
  cp "${INPUT_TAX_DB_DIR}" \
     "${SCRIPT_DIR}/Classifier/Pre-trained_classifier.qza"

  #Store the absolute path for use in downstream classification steps.
  CLASSIFIER="${SCRIPT_DIR}/Classifier/Pre-trained_classifier.qza"

  #Record which classifier was selected in the log.
  log "Using existing trained classifier: ${INPUT_TAX_DB_DIR}"

else

  ###########################################################
  ### OPTION B) TRAINING A NEW CLASSIFIER
  ###########################################################
  #Ask which reference database should be used for training.
  #SILVA and Greengenes2 require different reference preparation steps.
  echo "Which database should be used to train the classifier?"
  echo "1) SILVA 138.2"
  echo "2) GREENGENES2 2024.09"
  read -r -p "Enter SILVA or GREENGENES: " TAX_DB

  #Convert the input to uppercase so that lowercase input is also accepted.
  TAX_DB=$(echo "${TAX_DB}" | tr '[:lower:]' '[:upper:]')

  #Validate the input and convert numeric choices to database names.
  case "${TAX_DB}" in
    SILVA|1)
      TAX_DB="SILVA"
      ;;
    GREENGENES|GREENGREENES|GREENGENES2|2)
      TAX_DB="GREENGENES"
      ;;
    *)
      #Stop instead of accidentally training with an unintended database.
      echo "ERROR: Please enter SILVA or GREENGENES!"
      exit 1111
      ;;
  esac

  #Save the choice to automate the later parts of the code. 
  echo "${TAX_DB}" > TEMP_FILE_DB_CHOICE.txt

  #Create and enter the classifier output directory.
  #All intermediate files will be stored here to keep everything organized.
  mkdir -p "${SCRIPT_DIR}/Classifier"
  cd "${SCRIPT_DIR}/Classifier" || exit 1111

  #Define the forward and reverse primers used in the experiment.
  #Do NOT include sequencing adapters or library-preparation overhangs.
  #The primers must correspond to the region amplified in your samples!!
  F_PRIMER="ACTCCTACGGGAGGCAGCA"
  R_PRIMER="GGACTACHVGGGTWTCTAAT"

  #Prevent training from starting with invalid placeholder primer sequences.
  if [[ "${F_PRIMER}" == "FILL_WITH_FORWARD_PRIMER" ||
        "${R_PRIMER}" == "FILL_WITH_REVERSE_PRIMER" ]]; then
    echo "ERROR! Set F_PRIMER and R_PRIMER before training!!!!"
    exit 1111
  fi


  #==========================================================
  # DATABASE A) SILVA 138.2
  #==========================================================
  if [[ "${TAX_DB}" == "SILVA" ]]; then


    #--------------------------------------------------------
    # 1) Downloading the SILVA reference sequences and taxonomy
    #--------------------------------------------------------
    #Download the SILVA 138.2 SSURef_NR99 reference database.
    #SSURef_NR99 contains small-subunit ribosomal RNA sequences,
    # including bacterial, archaeal and eukaryotic reference sequences.
    log "Downloading SILVA 138.2 SSU reference sequences and taxonomy."

    qiime rescript get-silva-data \
      --p-version "138.2" \
      --p-target "SSURef_NR99" \
      --p-include-species-labels \
      --o-silva-sequences silva_ssu_seqs_138-2.qza \
      --o-silva-taxonomy silva_ssu_tax_138-2.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_ssu_seqs_138-2.qza" "SILVA sequences"
    check_file "silva_ssu_tax_138-2.qza" "SILVA taxonomy"


    #--------------------------------------------------------
    # 2) Converting RNA sequences to DNA sequences
    #--------------------------------------------------------
    #SILVA reference sequences may contain RNA uracil (U).
    #Convert U to T because the downstream DNA classifier workflow
    # and the amplicon reads use DNA nucleotide sequences.
    log "Converting SILVA RNA sequences to DNA sequences."

    qiime rescript reverse-transcribe \
      --i-rna-sequences silva_ssu_seqs_138-2.qza \
      --o-dna-sequences silva_ssu_seqs_138-2_DNA.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_ssu_seqs_138-2_DNA.qza" "SILVA DNA sequences"


    #--------------------------------------------------------
    # 3) Removing low-quality reference sequences
    #--------------------------------------------------------
    #Remove sequences with excessive ambiguous bases or long homopolymers.
    #This helps reduce problematic reference sequences before training.
    #The thresholds are quality-control choices and may remove some valid
    # biological sequences, so adjust them!!
    #CAREFUL: Processing a large reference database can use substantial RAM.
    log "Removing low-quality SILVA sequences."

    qiime rescript cull-seqs \
      --i-sequences silva_ssu_seqs_138-2_DNA.qza \
      --p-num-degenerates 5 \
      --p-homopolymer-length 8 \
      --p-n-jobs "${THREADS}" \
      --o-clean-sequences silva_138-2_seqs_cleaned.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_138-2_seqs_cleaned.qza" "Cleaned SILVA sequences"

    #--------------------------------------------------------
    # 4) Filtering sequences by taxon-specific length
    #--------------------------------------------------------
    #Remove sequences that are unusually short for their assigned domain.
    #The minimum lengths correspond to the order of labels:
    #Archaea = 900 bp, Bacteria = 1200 bp, Eukaryota = 1400 bp.
    #These thresholds are intended for near-full-length SSU references,
    #not for the short V3-V4 amplicons extracted later.
    #CAREFUL: This operation can be very RAM-intensive.
    log "Filtering SILVA sequences by taxon-specific length."

    qiime rescript filter-seqs-length-by-taxon \
      --i-sequences silva_138-2_seqs_cleaned.qza \
      --i-taxonomy silva_ssu_tax_138-2.qza \
      --p-labels Archaea Bacteria Eukaryota \
      --p-min-lens 900 1200 1400 \
      --o-filtered-seqs silva_138-2_cleaned_filtered.qza \
      --o-discarded-seqs silva_138-2_discarded.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_138-2_cleaned_filtered.qza" "Length-filtered SILVA sequences"


    #--------------------------------------------------------
    # 5) Extracting the V3-V4 amplicon region
    #--------------------------------------------------------
    #Extract only the region between the forward and reverse primers.
    #This makes the reference sequences comparable to your sample ASVs.
    #IMPORTANT: The primer sequences must match your actual experiment.
    #The length limits below are examples for a V3-V4 assay.
    #Adjust them if your expected amplicon length differs.
    log "Extracting the V3-V4 region from SILVA sequences."

    qiime feature-classifier extract-reads \
      --i-sequences silva_138-2_cleaned_filtered.qza \
      --p-f-primer "${F_PRIMER}" \
      --p-r-primer "${R_PRIMER}" \
      --p-min-length 350 \
      --p-max-length 600 \
      --p-n-jobs "${THREADS}" \
      --o-reads silva_138-2_v3v4_seqs.qza \
      --o-read-extraction-stats \
        "${QZA_DIR}/Read_v3v4_extraction_stats.qzv" \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_138-2_v3v4_seqs.qza" "Extracted SILVA V3-V4 sequences"


    #--------------------------------------------------------
    # 6) Dereplicating the extracted sequences and taxonomy
    #-------------------------------------------------------
    #Collapse duplicate sequences and consolidate their taxonomy.
    #LCA (lowest common ancestor) handling helps resolve cases where
    #identical sequences have different taxonomic assignments.
    #The output contains unique reference sequences and corresponding taxa.
    log "Dereplicating SILVA V3-V4 sequences and taxonomy."

    qiime rescript dereplicate \
      --i-sequences silva_138-2_v3v4_seqs.qza \
      --i-taxa silva_ssu_tax_138-2.qza \
      --p-mode "lca" \
      --p-threads "${THREADS}" \
      --o-dereplicated-sequences silva_138-2_v3v4_dereplicated.qza \
      --o-dereplicated-taxa silva_138-2_v3v4_tax_derep.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "silva_138-2_v3v4_dereplicated.qza" "Dereplicated SILVA sequences"
    check_file "silva_138-2_v3v4_tax_derep.qza" "Dereplicated SILVA taxonomy"


    #--------------------------------------------------------
    # 7) Training the Naive Bayes classifier
    #--------------------------------------------------------
    #Train a scikit-learn Naive Bayes classifier using the extracted
    #reference sequences and their matching taxonomic annotations.
    #The resulting classifier will assign taxonomy to your sample ASVs.
    #CAREFUL: This step can be very RAM-intensive and may take minutes
    #to hours, depending on the number of sequences and available hardware.
    log "Training the SILVA Naive Bayes classifier."

    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads silva_138-2_v3v4_dereplicated.qza \
      --i-reference-taxonomy silva_138-2_v3v4_tax_derep.qza \
      --o-classifier silva_138-2_v3v4_classifier.qza \
      --verbose \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    #Verify that the trained classifier was created successfully.
    check_file "silva_138-2_v3v4_classifier.qza" "Trained SILVA classifier"

    #Store the final classifier path for downstream pipeline steps.
    CLASSIFIER="${SCRIPT_DIR}/Classifier/silva_138-2_v3v4_classifier.qza"



  #==========================================================
  # DATABASE B) GREENGENES2 2024.09
  #==========================================================
  else


    #--------------------------------------------------------
    # 1) Downloading the Greengenes2 reference data
    #--------------------------------------------------------
    #Using the Greengenes2 2024.09 release.
    #Greengenes2 distributes multiple reference files, including
    # backbone and larger sequence collections.
    #For classifier training, the sequence and taxonomy files must
    # refer to the same records and use compatible identifiers.
    #IMPORTANT: you must verify the exact filenames in the release directory
    # before downloading; filenames can differ between file variants!
    log "Preparing the Greengenes2 2024.09 reference database."

    GG_URL="https://ftp.microbio.me/greengenes_release/2024.09"

    #The filenames must be verified in case the official release changes
    # before using these download commands.
    wget -c "${GG_URL}/2024.09.backbone.full-length.fna.qza" \
      -O gg2_backbone_full_length.qza || exit 1111

    #Taxonomy file
    wget -c "${GG_URL}/2024.09.backbone.tax.qza" \
      -O gg2_backbone_taxonomy.qza || exit 1111

    check_file "gg2_backbone_full_length.qza" "Greengenes2 backbone sequences"
    check_file "gg2_backbone_taxonomy.qza" "Greengenes2 backbone taxonomy"


    #--------------------------------------------------------
    # 2) Extracting the V3-V4 amplicon region
    #--------------------------------------------------------
    #Extract the same primer-defined region used for SILVA.
    #Using the same primers makes the two classifiers comparable
    #for the same sequencing assay.
    #The reference sequence IDs must be preserved so they can be
    #matched correctly to their taxonomy annotations.
    log "Extracting the V3-V4 region from Greengenes2 sequences."

    qiime feature-classifier extract-reads \
      --i-sequences gg2_backbone_full_length.qza \
      --p-f-primer "${F_PRIMER}" \
      --p-r-primer "${R_PRIMER}" \
      --p-min-length 350 \
      --p-max-length 600 \
      --p-n-jobs "${THREADS}" \
      --o-reads gg2_2024-09_v3v4_seqs.qza \
      --o-read-extraction-stats \
        "${QZA_DIR}/GG2_v3v4_extraction_stats" \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "gg2_2024-09_v3v4_seqs.qza" "Extracted Greengenes2 V3-V4 sequences"


    #--------------------------------------------------------
    # 3) Matching taxonomy and dereplicating sequences
    #--------------------------------------------------------
    #IMPORTANT: Before proceeding, confirm that the taxonomy artifact
    # contains the identifiers expected by the reference sequences.
    #
    #If the taxonomy file uses a different identifier scheme, map or
    # convert it to the appropriate QIIME 2 taxonomy format first.
    #We cannot assume that all Greengenes2 taxonomy files are interchangeable!
    #The sequence and taxonomy IDs must correspond to the same records.
    #
    #After confirming compatibility, dereplicate the extracted sequences
    # and resolve conflicting taxonomy assignments using LCA.
    log "Dereplicating Greengenes2 V3-V4 sequences and taxonomy."

    qiime rescript dereplicate \
      --i-sequences gg2_2024-09_v3v4_seqs.qza \
      --i-taxa gg2_backbone_taxonomy.qza \
      --p-mode "lca" \
      --p-threads "${THREADS}" \
      --o-dereplicated-sequences gg2_2024-09_v3v4_dereplicated.qza \
      --o-dereplicated-taxa gg2_2024-09_v3v4_tax_derep.qza \
      2>&1 | tee -a "${LOGFILE}" || exit 1

    check_file "gg2_2024-09_v3v4_dereplicated.qza" "Dereplicated Greengenes2 sequences"
    check_file "gg2_2024-09_v3v4_tax_derep.qza" "Dereplicated Greengenes2 taxonomy"


    #--------------------------------------------------------
    # 4) Training the Greengenes2 Naive Bayes classifier
    #--------------------------------------------------------
    #Train the classifier using the extracted Greengenes2 sequences
    #and their matching taxonomy annotations.
    #CAREFUL: Model training can require substantial RAM and CPU time.
    #Training will fail if the reference sequence and taxonomy artifacts
    #are incompatible or if the required software versions conflict.
    log "Training the Greengenes2 2024.09 Naive Bayes classifier."

    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads gg2_2024-09_v3v4_dereplicated.qza \
      --i-reference-taxonomy gg2_2024-09_v3v4_tax_derep.qza \
      --o-classifier gg2_2024-09_v3v4_classifier.qza \
      --verbose \
      2>&1 | tee -a "${LOGFILE}" || exit 1111

    #Verify that the trained Greengenes2 classifier exists.
    check_file "gg2_2024-09_v3v4_classifier.qza" "Trained Greengenes2 classifier"

    #Store the final classifier path for use by downstream steps.
    CLASSIFIER="${SCRIPT_DIR}/Classifier/gg2_2024-09_v3v4_classifier.qza"

  fi

  #Return to the main project directory.
  cd "${SCRIPT_DIR}" || exit 1111

  log "Classifier training finished: ${CLASSIFIER}"

fi


#########################################################
### STEP 7: ASV taxonomy classification using classifier  
#########################################################
log "STEP 7: ASV taxonomy classification using classifier."


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