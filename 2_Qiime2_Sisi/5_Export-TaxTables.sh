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
#    4.   Base quality assessment -> trimmed_demux.qzv (manual review in Qiime2 View)
#    5.   Quality filtering/denoising using DADA2 -> AND Denoising_stats.tsv
#    6.   Train Naive Bayes classifier for taxonomic assignation -> silva_144_v3v5_classifier.qza
#    7.   ASV taxonomy classification using DADA2 -> Taxonomy.qza
#    8.   Rarefaction analysis.
#    8a.  Finding maximum depth.
#    8b.  ALPHA RAREFACTION (observed_features + shannon)
#    8c.  Build the MASTER per-sample table (reads at every stage + metadata)
#    8d.  Static plots + candidate depths (using R)
#    9.   Choosing sampling depth (manually)
#    9a.  Rarefy the feature table
#    9b.  Collapse rarefied table at every taxonomic level and export
#    9c.  Export raw/uncollapsed ASV level rarefied table + per-ASV taxonomy needed to analysis in R "asv" mode
#   *10.  Export ASV and genus level rarefied and raw tables.
#   *10a. Export non-rarefied ASV table
#   *10b. Collapse and export non-rarefied genus table
#   *10c. Export rarefied ASV table
#   *10d. Collapse and export rarefied genus table
#   *10e. Export ASV taxonomy mapping
#    11.  Building phylogenetic tree from ASV representative sequences.
#    12. Core metrics (only at ASV-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#
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
FINAL="${SCRIPT_DIR}/Final_results"

#New folder! The new structure is basically: 
#Final_results
#   FINAL_tables
#       RAW
#       RAREFIED 

FINAL_TABLES_DIR="${FINAL}/FINAL_tables"
RAW_DIR="${FINAL_TABLES_DIR}/RAW"
RAREFIED_DIR="${FINAL_TABLES_DIR}/RAREFIED"

TAXONOMY_QZA="${QZA_DIR}/Taxonomy.qza"

#Existing files to use!
RAW_ASV_QZA="${QZA_DIR}/ASV_table.qza"
RAREFIED_ASV_QZA="${QZA_DIR}/ASV_table_rarefied.qza"

#Final output .tsv frequency tables!
RAW_ASV_TSV="${RAW_DIR}/ASV_table_RAW.tsv"
RAW_GENUS_TSV="${RAW_DIR}/GENUS_table_RAW.tsv"

RAREFIED_ASV_TSV="${RAREFIED_DIR}/ASV_table_rarefied.tsv"
RAREFIED_GENUS_TSV="${RAREFIED_DIR}/GENUS_table_rarefied.tsv"



########################################################
## 10. Export ASV and genus level rarefied and raw tables. 
########################################################
log "STEP 10: Exporting raw and rarefied ASV/genus tables."

#In truth more than exporting we are just transforming the exported tables from script 4_Rarefaction.sh 
# and putting them somewhere easily accessible. 

#Also!!! As a reminder, ASV does not get taxonomy so you can export it but do NOT 
# try to use taxonomy with ASV tables. 

########################################################
## Create output directories
########################################################
mkdir -p "${RAW_DIR}" "${RAREFIED_DIR}"
mkdir -p "${QZA_DIR}/Genus_taxa_tables"

########################################################
## Validate required inputs
########################################################
check_file "${RAW_ASV_QZA}" "Non-rarefied ASV feature table"
check_file "${RAREFIED_ASV_QZA}" "Rarefied ASV feature table (from the 4_Rarefaction.sh script)"
check_file "${TAXONOMY_QZA}" "ASV taxonomy assignments"



########################################################
## STEP 10a: Export non-rarefied ASV table
########################################################
log "STEP 10a: Exporting non-rarefied ASV table."

EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${RAW_ASV_QZA}" \
  --output-path "${EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${EXPORT_DIR}/feature-table.biom" \
  "Exported non-rarefied ASV BIOM table"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${RAW_ASV_TSV}" \
  --to-tsv \
  2>&1 | tee -a "${LOGFILE}"

rm -rf "${EXPORT_DIR}"

check_file "${RAW_ASV_TSV}" "Non-rarefied ASV TSV table"

log "FINISHED STEP 10a: Exporting non-rarefied ASV table."



########################################################
## STEP 10b: Collapse and export non-rarefied genus table
########################################################
log "STEP 10b: Collapsing and exporting non-rarefied genus table."

RAW_GENUS_QZA="${QZA_DIR}/Genus_taxa_tables/table-level7-genus-RAW.qza"

qiime taxa collapse \
  --i-table "${RAW_ASV_QZA}" \
  --i-taxonomy "${TAXONOMY_QZA}" \
  --p-level 7 \
  --o-collapsed-table "${RAW_GENUS_QZA}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${RAW_GENUS_QZA}" "Non-rarefied genus-level QZA table"

EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${RAW_GENUS_QZA}" \
  --output-path "${EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${EXPORT_DIR}/feature-table.biom" \
  "Exported non-rarefied genus BIOM table"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${RAW_GENUS_TSV}" \
  --to-tsv \
  2>&1 | tee -a "${LOGFILE}"

rm -rf "${EXPORT_DIR}"

check_file "${RAW_GENUS_TSV}" "Non-rarefied genus TSV table"

log "FINISHED STEP 10b: Collapsing and exporting non-rarefied genus table."



########################################################
## STEP 10c: Export rarefied ASV table
########################################################
log "STEP 10c: Exporting rarefied ASV table."

EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${RAREFIED_ASV_QZA}" \
  --output-path "${EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${EXPORT_DIR}/feature-table.biom" \
  "Exported rarefied ASV BIOM table"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${RAREFIED_ASV_TSV}" \
  --to-tsv \
  2>&1 | tee -a "${LOGFILE}"

rm -rf "${EXPORT_DIR}"

check_file "${RAREFIED_ASV_TSV}" "Rarefied ASV TSV table"

log "FINISHED STEP 10c: Exporting rarefied ASV table."



########################################################
## STEP 10d: Collapse and export rarefied genus table
########################################################
log "STEP 10d: Collapsing and exporting rarefied genus table."

RAREFIED_GENUS_QZA="${QZA_DIR}/Genus_taxa_tables/table-level7-genus-rarefied.qza"

qiime taxa collapse \
  --i-table "${RAREFIED_ASV_QZA}" \
  --i-taxonomy "${TAXONOMY_QZA}" \
  --p-level 7 \
  --o-collapsed-table "${RAREFIED_GENUS_QZA}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${RAREFIED_GENUS_QZA}" "Rarefied genus-level QZA table"

EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${RAREFIED_GENUS_QZA}" \
  --output-path "${EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${EXPORT_DIR}/feature-table.biom" \
  "Exported rarefied genus BIOM table"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${RAREFIED_GENUS_TSV}" \
  --to-tsv \
  2>&1 | tee -a "${LOGFILE}"

rm -rf "${EXPORT_DIR}"

check_file "${RAREFIED_GENUS_TSV}" "Rarefied genus TSV table"

log "FINISHED STEP 10d: Collapsing and exporting rarefied genus table."




########################################################
## STEP 10e: Export GENUS and ASV (modified) taxonomy mapping
########################################################

log "STEP 10e: Exporting GENUS and ASV taxonomy assignments."


########################################################
## PART 1: Export the original QIIME 2 taxonomy table
########################################################

#The original Taxonomy.qza is produced by the QIIME 2
# classify-sklearn function.
#
#In this dataset, the classifier only assigns taxonomy up to
# genus for many ASVs. This is not an error in the export.
#
#We preserve the original taxonomy table exactly as QIIME 2
# generated it. It contains:
#
#   Feature ID  = ASV identifier
#   Taxon       = taxonomy assigned to the ASV
#   Confidence  = classifier confidence score
#
#We call these files GENUS_taxonomy.tsv to make clear that
# the taxonomy assignments may stop at genus.
#
#The same taxonomy mapping is copied to RAW and RAREFIED
# because rarefaction changes abundance counts, not taxonomy.

TAXONOMY_EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${TAXONOMY_QZA}" \
  --output-path "${TAXONOMY_EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${TAXONOMY_EXPORT_DIR}/taxonomy.tsv" \
  "Exported taxonomy TSV (assignments up to genus where available)"

cp "${TAXONOMY_EXPORT_DIR}/taxonomy.tsv" \
  "${RAW_DIR}/GENUS_taxonomy.tsv"

cp "${TAXONOMY_EXPORT_DIR}/taxonomy.tsv" \
  "${RAREFIED_DIR}/GENUS_taxonomy.tsv"

check_file "${RAW_DIR}/GENUS_taxonomy.tsv" \
  "Raw GENUS taxonomy mapping"

check_file "${RAREFIED_DIR}/GENUS_taxonomy.tsv" \
  "Rarefied GENUS taxonomy mapping"



########################################################
## PART 2: Export the DNA sequences of each ASV
########################################################

#In order to create a TSV containing both taxonomy and the
# DNA sequence of each ASV, we also need the representative
# sequences stored in ASV_seqs.qza.
#
#QIIME 2 exports these sequences in FASTA format.
#
#Each FASTA header contains the ASV Feature ID, which allows
# us to match the DNA sequence to the corresponding taxonomy.
#
#This makes it possible to identify significant ASVs later
# and investigate their sequences using BLAST.
#
#IMPORTANT:
# Adding the sequence does not automatically assign taxonomy
# below genus. It simply gives us the DNA sequence for further
# identification!!

SEQUENCES_EXPORT_DIR="$(mktemp -d)"

qiime tools export \
  --input-path "${QZA_DIR}/ASV_seqs.qza" \
  --output-path "${SEQUENCES_EXPORT_DIR}" \
  2>&1 | tee -a "${LOGFILE}"

check_file "${SEQUENCES_EXPORT_DIR}/dna-sequences.fasta" \
  "Exported representative ASV sequences"



########################################################
## PART 3: Match each ASV with its taxonomy and DNA sequence
########################################################

#We create a modified taxonomy table with four columns:
#
#   Feature ID | ASV_sequence | Taxon | Confidence
#
#The merge is performed using the ASV Feature ID.
#
#The AWK script works as follows:
#
#   1. Read taxonomy.tsv and store each ASV's taxonomy and
#      confidence using its Feature ID as the key.
#
#   2. Read dna-sequences.fasta.
#
#   3. Extract the Feature ID from each FASTA header.
#
#   4. Join the sequence with the corresponding taxonomy.
#
#   5. Write one row per ASV to the final TSV.
#
#FASTA sequences may span several lines, so the script
# concatenates all sequence lines belonging to each ASV.
#
#We save the output initially in RAW_DIR. After validation,
# we copy the same mapping into RAREFIED_DIR because taxonomy
# and sequences are identical for both abundance tables.

ASV_TAXONOMY_SEQUENCE_TSV="${RAW_DIR}/ASV_taxonomy_with_sequences.tsv"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}

########################################################
# Read taxonomy.tsv first.
# FNR == NR identifies the first input file.
########################################################

FNR == NR {

    # Read the taxonomy table column names.
    if (FNR == 1) {
        for (i = 1; i <= NF; i++) {
            # Remove a possible Windows carriage return.
            gsub(/\r$/, "", $i)

            if ($i == "Feature ID") id_col = i
            if ($i == "Taxon") taxon_col = i
            if ($i == "Confidence") confidence_col = i
        }

        # Stop if the expected columns are missing.
        if (!id_col || !taxon_col || !confidence_col) {
            print "ERROR: Expected taxonomy columns not found." \
                > "/dev/stderr"
            failed = 1
            exit 1
        }

        next
    }

    # Ignore empty lines.
    if (NF == 0) next

    # Store taxonomy and confidence indexed by Feature ID.
    feature_id = $id_col
    gsub(/\r$/, "", feature_id)

    taxonomy[feature_id] = $taxon_col
    confidence[feature_id] = $confidence_col

    next
}

########################################################
# Read the FASTA file.
# Each header begins with > and contains the ASV Feature ID.
########################################################

/^>/ {

    # Before reading a new sequence, print the previous one.
    if (feature_id != "") {
        print_sequence()
    }

    # Extract the identifier from the FASTA header.
    feature_id = substr($0, 2)

    # Remove any optional description after the identifier.
    sub(/[[:space:]].*/, "", feature_id)

    # Start a new sequence.
    sequence = ""

    next
}

########################################################
# Append sequence lines.
# This also supports multiline FASTA sequences.
########################################################

{
    gsub(/[[:space:]\r]/, "")
    sequence = sequence $0
}

########################################################
# Finish processing the final FASTA record.
########################################################

END {

    if (feature_id != "") {
        print_sequence()
    }

    # Fail if any FASTA record could not be matched to taxonomy.
    if (unmatched > 0) {
        printf "ERROR: %d ASV sequences had no taxonomy match.\n", \
            unmatched > "/dev/stderr"
        failed = 1
    }

    # Fail if no ASV sequences were processed.
    if (sequence_count == 0) {
        print "ERROR: No ASV sequences were successfully matched." \
            > "/dev/stderr"
        failed = 1
    }

    if (failed) {
        exit 1
    }

    printf "Successfully matched %d ASV sequences to taxonomy.\n", \
        sequence_count > "/dev/stderr"
}

########################################################
# Function to write one ASV record.
########################################################

function print_sequence() {

    # Confirm that the Feature ID exists in taxonomy.tsv.
    if (!(feature_id in taxonomy)) {
        print "ERROR: No taxonomy match for " feature_id \
            > "/dev/stderr"
        unmatched++
        return
    }

    # Write the column names once.
    if (!header_printed) {
        print "Feature ID", "ASV_sequence", "Taxon", "Confidence"
        header_printed = 1
    }

    # Write one complete row per ASV.
    print feature_id, sequence, taxonomy[feature_id], confidence[feature_id]

    sequence_count++
}
' "${TAXONOMY_EXPORT_DIR}/taxonomy.tsv" \
  "${SEQUENCES_EXPORT_DIR}/dna-sequences.fasta" \
  > "${ASV_TAXONOMY_SEQUENCE_TSV}"

check_file "${ASV_TAXONOMY_SEQUENCE_TSV}" \
  "ASV taxonomy mapping with DNA sequences"


########################################################
## PART 4: Copy the modified taxonomy table to RAREFIED and RAW folders
########################################################

#Originally the ASV_taxonomy_with_sequences.tsv is in the RAW directory so we only
# need to copy it to the RAREFIED folder. 

cp "${ASV_TAXONOMY_SEQUENCE_TSV}" \
  "${RAREFIED_DIR}/ASV_taxonomy_with_sequences.tsv"

check_file "${RAW_DIR}/ASV_taxonomy_with_sequences.tsv" \
  "Raw ASV taxonomy mapping with DNA sequences"

check_file "${RAREFIED_DIR}/ASV_taxonomy_with_sequences.tsv" \
  "Rarefied ASV taxonomy mapping with DNA sequences"



########################################################
## PART 5: Clean up temporary directories
########################################################

#Remove temporary export folders after the output files
# have been created and validated.

rm -rf "${TAXONOMY_EXPORT_DIR}" "${SEQUENCES_EXPORT_DIR}"

log "FINISHED STEP 10e: Exporting GENUS and ASV taxonomy mappings."
log "FINISHED STEP 10: Exporting raw and rarefied ASV/genus tables."



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 5_Export-TaxTables.sh finished running!                                                             "
log "                                                                                                           "
log "ASV AND GENUS TABLE EXPORT COMPLETE. Key outputs:                                                          "
log "                                                                                                           "
log "  Exported abundance tables (TSV):                                                                         "
log "    - ${RAW_DIR}/ASV_table_RAW.tsv                    (Raw, non-rarefied ASV counts)                       "
log "    - ${RAW_DIR}/GENUS_table_RAW.tsv                  (Raw, non-rarefied genus counts)                     "
log "    - ${RAREFIED_DIR}/ASV_table_rarefied.tsv          (Rarefied ASV counts)                                "
log "    - ${RAREFIED_DIR}/GENUS_table_rarefied.tsv        (Rarefied genus counts)                              "
log "                                                                                                           "
log "  Original taxonomy assignments:                                                                           "
log "    - ${RAW_DIR}/GENUS_taxonomy.tsv                   (QIIME 2 taxonomy assignments)                       "
log "    - ${RAREFIED_DIR}/GENUS_taxonomy.tsv               (Same taxonomy assignments for rarefied data)       "
log "                                                                                                           "
log "  ASV taxonomy with representative DNA sequences:                                                          "
log "    - ${RAW_DIR}/ASV_taxonomy_with_sequences.tsv      (Feature ID, DNA sequence, taxonomy, confidence)     "
log "    - ${RAREFIED_DIR}/ASV_taxonomy_with_sequences.tsv (Same mapping for rarefied data)                      "
log "                                                                                                           "
log "  Output directories:                                                                                      "
log "    - ${RAW_DIR}/                                     (Non-rarefied tables and taxonomy mappings)          "
log "    - ${RAREFIED_DIR}/                                (Rarefied tables and taxonomy mappings)              "
log "                                                                                                           "
log "==========================================================================================================="
log "NEXT STEPS:"
log "  1. Check that the exported tables contain the expected samples and counts."
log "  2. Use RAW tables for analyses requiring non-rarefied counts, such as ANCOM-BC2."
log "  3. Use RAREFIED tables for analyses requiring standardized sequencing depth."
log "  4. Use ASV_taxonomy_with_sequences.tsv to retrieve the DNA sequence of a significant ASV."
log "  5. Investigate selected ASV sequences using BLAST or another sequence-identification method."
log "  6. Proceed to downstream R analysis or phylogenetic diversity calculations."
log "==========================================================================================================="
