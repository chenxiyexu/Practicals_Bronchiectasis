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
#   *1.   Import data -> demux.qza
#   *2.   Import metadata (.tsv) -> metadata.qzv
#   *3.   Trim primers/adapters using Cutadapt -> trimmed_demux.qza
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
#Change this depending on which PC you use (must be a multiple of 2)
# (nproc --all #this shows total amount of cores/CPUs in a linux terminal).
THREADS=8

########################################################
## Variables needed!!!
########################################################
MANIFEST_DIR="${SCRIPT_DIR}/../1_SRA_Download/Manifest_Qiime2.tsv"




#########################################################
### STEP 1: Importing data and adapting it to .qza
#########################################################
log "STEP 1: Importing data and adapting it to .qza."

mkdir -p "Output_qza" #We group all .qza here for legibility. 
QZA_DIR="${SCRIPT_DIR}/Output_qza"

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-path "${MANIFEST_DIR}" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path "${QZA_DIR}/demux.qza"
  
check_file "${QZA_DIR}/demux.qza" "Finished importing demultiplexed sequences!"
  
  
mkdir -p "Output_qzv" #We group here all .qzv for legibility (Qiime2 View).
QZV_DIR="${SCRIPT_DIR}/Output_qzv"

qiime demux summarize \
  --i-data "${QZA_DIR}/demux.qza" \
  --o-visualization "${QZV_DIR}/demux.qzv"  #Check the raw data in Qiime2 View

check_file "${QZV_DIR}/demux.qzv" "Please use Qiime2 view to revise the quality of the reads!" 
  
log "FINISHED STEP 1: Importing and adapting data from .fastq to .qza."



#########################################################
### STEP 2: Visualizing metadata as a .qzv
#########################################################
log "STEP 2: Importing metadata and visualizing as a .qzv"

echo "Please enter the name of your metadata file (must be in the same directory as this file): "
read -r METADATA

#Check if it exists
if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: ${METADATA} metadata file not found!"
    exit 1111
fi

log "Using metadata file: ${METADATA}"

qiime metadata tabulate \
  --m-input-file "${METADATA}" \
  --o-visualization "${QZV_DIR}/metadata.qzv" #Check the raw data in Qiime2 View

log "FINISHED STEP 2: Importing and adapting metadata."



#########################################################
### STEP 3: Trim off adapters/primers using cutadapt
#########################################################
log "STEP 3: Trim off adapters/primers using cutadapt (for Paired-end read files)."

#Includes PRIMER and OVERHANG (from Ilumina)
#--p-front-f "FILL WITH RELEVANT PRIMERS AND OVERHANGS IF NEEDED!!!!" \   Forward primer in 5'-3' direction.
#--p-front-r "FILL WITH RELEVANT PRIMERS AND OVERHANGS IF NEEDED!!!!" \   Reverse primer in 5'-3' direction.
log "NOTE: The primers used in this step are the ones used in the original experiments. If you used different primers, please change them accordingly in the script."
log "This step includes removing PRIMERS and OVERHANGS (from Illumina in this case)."


qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "${QZA_DIR}/demux.qza" \
  --p-front-f "ACTCCTACGGGAGGCAGCA" \
  --p-front-r "GGACTACHVGGGTWTCTAAT" \
  --p-match-adapter-wildcards \
  --p-discard-untrimmed \
  --p-no-indels \
  --p-cores "${THREADS:-1}" \
  --o-trimmed-sequences "${QZA_DIR}/demux_trimmed.qza" \
  --o-stats "${QZA_DIR}/trim_stats.qza" \
  --verbose 2>&1 | tee -a "${LOGFILE}"

check_file "${QZA_DIR}/demux_trimmed.qza" "Sequences with the primer trimmed!"

qiime demux summarize \
  --i-data "${QZA_DIR}/demux_trimmed.qza" \
  --o-visualization "${QZV_DIR}/demux_trimmed.qzv"
  
check_file "${QZV_DIR}/demux_trimmed.qzv" "Trimmed demux files summary visualization"

log "FINISHED STEP 3: Trimmed off primers."



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 1_Import-QualityCheck.sh finished running!                                                          "
log "                                                                                                           "
log "DATA IMPORT, METADATA PREPARATION, AND PRIMER TRIMMING COMPLETE. Key outputs:                              "
log "  Artifacts (.qza):                                                                                        "
log "    - ${QZA_DIR}/demux.qza                          (Raw imported paired-end sequence data)                "
log "    - ${QZA_DIR}/demux_trimmed.qza                  (Sequence data with primers trimmed via Cutadapt)      "
log "  Visualizations (.qzv - view at view.qiime2.org):                                                         "
log "    - ${QZV_DIR}/demux.qzv                          (Quality distribution of raw reads)                    "
log "    - ${QZV_DIR}/metadata.qzv                       (Tabulated view of sample metadata)                    "
log "    - ${QZV_DIR}/demux_trimmed.qzv                  (Quality distribution of trimmed reads)                "
log "==========================================================================================================="
log "NEXT STEPS:                                                                                                "
log "  1. Inspect ${QZV_DIR}/demux_trimmed.qzv at https://view.qiime2.org to evaluate base quality.             "
log "  2. Determine truncation positions for quality filtering.                                                 "
log "  3. Proceed to Script 2 for DADA2 denoising and taxonomic classification.                                 "
log "==========================================================================================================="