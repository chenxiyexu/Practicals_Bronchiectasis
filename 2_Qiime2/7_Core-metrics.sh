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
#   *12. Core metrics (at ASV-level, for exploratory purposes, for a more complete analysis go to the R scripts).
#   *13. Core metrics (at genus-level, for exploratory purposes, for a more complete analysis go to the R scripts).
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
TREE_QZA_DIR="${QZA_DIR}/Tree"




########################################################
## Create separate output directories
########################################################

#Main output directory for all core metrics.
CORE_METRICS_DIR="${FINAL}/Core_metrics_Qiime2"

#Separate directories for ASV-level and genus-level analyses.
ASV_CORE_DIR="${CORE_METRICS_DIR}/ASV"
GENUS_CORE_DIR="${CORE_METRICS_DIR}/GENUS"

#Delete previous results automatically.
#WARNING: This deletes all previous ASV and genus results.
rm -rf "${ASV_CORE_DIR}" "${GENUS_CORE_DIR}"

#Create only the main parent directory.
#Do NOT create ASV_CORE_DIR or GENUS_CORE_DIR
#yet, because QIIME 2 will create these directories itself.
mkdir -p "${CORE_METRICS_DIR}"

#Define separate export directories without creating them.
ASV_CORE_EXP="${ASV_CORE_DIR}/Core_metrics_Qiime2_EXPORT"
GENUS_CORE_EXP="${GENUS_CORE_DIR}/Core_metrics_Qiime2_EXPORT"

#Define separate visualization directories without creating them.
ASV_CORE_QZV="${ASV_CORE_DIR}/Core_metrics_QZV"
GENUS_CORE_QZV="${GENUS_CORE_DIR}/Core_metrics_QZV"

#Ask the user to choose the sampling depth.
echo "Enter the sampling depth for core metrics for genus and ASV-level comparisons (example: 8000): "
read -r SAMPLING_DEPTH

#Check that the sampling depth is a positive integer.
if ! [[ "${SAMPLING_DEPTH}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR! Sampling depth must be a positive integer."
    exit 1
fi




########################################################
## STEP 12: Core metrics for ASVs (for exploratory purposes,
##          for a more complete analysis go to the R scripts).
########################################################
log "STEP 12: Core metrics for ASV-level exploratory analysis."

#ASV_table.qza is RAW (ASV_table_rarefied.qza has rarefaction already applied!).
#QIIME 2 will rarefy the RAW table internally to SAMPLING_DEPTH.

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${TREE_QZA_DIR}/Rooted_tree.qza" \
  --i-table "${QZA_DIR}/ASV_table.qza" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${METADATA}" \
  --p-n-jobs-or-threads "${THREADS}" \
  --output-dir "${ASV_CORE_DIR}" \
  2>&1 | tee -a "${LOGFILE}"


########################################################
## Create folders for exports and QZVs.
########################################################

#Create a folder for exported ASV-level diversity tables
#and distance matrices (TSV and other standard formats).
ASV_CORE_EXP="${ASV_CORE_DIR}/Core_metrics_Qiime2_EXPORT"
mkdir -p "${ASV_CORE_EXP}"

#Create a folder for ASV-level visualizations (.qzv).
ASV_CORE_QZV="${ASV_CORE_DIR}/Core_metrics_QZV"
mkdir -p "${ASV_CORE_QZV}"


########################################################
##Export core metrics to universal formats (.tsv)
########################################################

### Alpha-diversity
## Shannon
qiime tools export \
  --input-path "${ASV_CORE_DIR}/shannon_vector.qza" \
  --output-path "${ASV_CORE_EXP}/shannon" \
  2>&1 | tee -a "${LOGFILE}"

## Observed features
qiime tools export \
  --input-path "${ASV_CORE_DIR}/observed_features_vector.qza" \
  --output-path "${ASV_CORE_EXP}/observed" \
  2>&1 | tee -a "${LOGFILE}"

## Faith's Phylogenetic Diversity
qiime tools export \
  --input-path "${ASV_CORE_DIR}/faith_pd_vector.qza" \
  --output-path "${ASV_CORE_EXP}/faith" \
  2>&1 | tee -a "${LOGFILE}"


### Beta-diversity
## Bray-Curtis
qiime tools export \
  --input-path "${ASV_CORE_DIR}/bray_curtis_distance_matrix.qza" \
  --output-path "${ASV_CORE_EXP}/bray" \
  2>&1 | tee -a "${LOGFILE}"

## Jaccard
qiime tools export \
  --input-path "${ASV_CORE_DIR}/jaccard_distance_matrix.qza" \
  --output-path "${ASV_CORE_EXP}/jaccard" \
  2>&1 | tee -a "${LOGFILE}"

## Weighted UniFrac
qiime tools export \
  --input-path "${ASV_CORE_DIR}/weighted_unifrac_distance_matrix.qza" \
  --output-path "${ASV_CORE_EXP}/weighted" \
  2>&1 | tee -a "${LOGFILE}"

## Unweighted UniFrac
qiime tools export \
  --input-path "${ASV_CORE_DIR}/unweighted_unifrac_distance_matrix.qza" \
  --output-path "${ASV_CORE_EXP}/unweighted" \
  2>&1 | tee -a "${LOGFILE}"


########################################################
## Alpha-diversity comparisons between groups.
########################################################

### Faith's Phylogenetic Diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${ASV_CORE_DIR}/faith_pd_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${ASV_CORE_QZV}/faith-pd-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"

### Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${ASV_CORE_DIR}/evenness_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${ASV_CORE_QZV}/evenness-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"

### Number of observed features
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${ASV_CORE_DIR}/observed_features_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${ASV_CORE_QZV}/observed-features-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"

### Shannon diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${ASV_CORE_DIR}/shannon_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${ASV_CORE_QZV}/shannon-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"


########################################################
## Beta-diversity group significance: comparing groups.
########################################################
#Here, sex is used as the grouping variable.
#Replace "sex" with another metadata column when needed.
#
#IMPORTANT:
#The metadata column must contain at least two groups,
#and the samples must match those in the distance matrix.

### Weighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix "${ASV_CORE_DIR}/weighted_unifrac_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${ASV_CORE_QZV}/weighted-unifrac-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"

### Unweighted UniFrac
qiime diversity beta-group-significance \
  --i-distance-matrix "${ASV_CORE_DIR}/unweighted_unifrac_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${ASV_CORE_QZV}/unweighted-unifrac-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"

### Bray-Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix "${ASV_CORE_DIR}/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${ASV_CORE_QZV}/bray-curtis-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"

### Jaccard
qiime diversity beta-group-significance \
  --i-distance-matrix "${ASV_CORE_DIR}/jaccard_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${ASV_CORE_QZV}/jaccard-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"


########################################################
## ANCOM-BC2: comparing ASV abundances between groups.
## IMPORTANT: Use RAW (non-rarefied) data.
########################################################

#The raw, non-rarefied ASV table.
RAW_ASV_QZA="${QZA_DIR}/ASV_table.qza"

#Create a temporary directory for ANCOM-BC2 intermediate files.
ANCOM_DIR="$(mktemp -d)"

#Step 1: Filter low-abundance and low-prevalence ASVs.
#Here, features must have at least 50 total counts and occur
#in at least 2 samples.

#IMPORTANT:
#Filtering changes the set of features tested. Choose thresholds
#appropriate for your study and report them in your methods.

qiime feature-table filter-features \
  --i-table "${RAW_ASV_QZA}" \
  --p-min-frequency 50 \
  --p-min-samples 2 \
  --o-filtered-table "${ANCOM_DIR}/ANCOMBC2_ASV_Filtered.qza" \
  2>&1 | tee -a "${LOGFILE}"


#Step 2: Run ANCOM-BC2 to compare ASV abundance between sex groups.
#The model uses sex as the fixed effect.
qiime composition ancombc2 \
  --i-table "${ANCOM_DIR}/ANCOMBC2_ASV_Filtered.qza" \
  --m-metadata-file "${METADATA}" \
  --p-fixed-effects-formula 'sex' \
  --o-ancombc2-output "${ANCOM_DIR}/ANCOMBC2_ASV_results.qza" \
  2>&1 | tee -a "${LOGFILE}"


#Step 3: Copy the ANCOM-BC2 results into the ASV output directory.
cp "${ANCOM_DIR}/ANCOMBC2_ASV_results.qza" \
   "${ASV_CORE_DIR}/ANCOMBC2_ASV_results.qza"


#Step 4: Create a differential-abundance barplot.
#Taxonomic labels will be limited to the assignments available
#in Taxonomy.qza (currently up to genus in your dataset).
qiime composition da-barplot \
  --i-data "${ANCOM_DIR}/ANCOMBC2_ASV_results.qza" \
  --i-taxonomy "${QZA_DIR}/Taxonomy.qza" \
  --o-visualization "${ASV_CORE_QZV}/ANCOMBC2_ASV_barplot.qzv" \
  2>&1 | tee -a "${LOGFILE}"


#Remove the temporary directory after successful completion.
rm -rf "${ANCOM_DIR}"

log "FINISHED STEP 12: Core metrics for ASV-level exploratory analysis."

##IMPORTANT!!!! Because Taxonomy.qza only contains labels up to genus the ASV results will
# only have genus labels!












########################################################
## STEP 13: Core metrics for genus-level data.
##          Exploratory analysis; continue with R scripts
##          for more comprehensive statistical analyses.
########################################################
log "STEP 13: Core metrics for genus-level exploratory analysis."



########################################################
## Define input and output paths.
########################################################

# Raw, non-rarefied genus abundance table.
GENUS_TABLE="${QZA_DIR}/Genus_taxa_tables/table-level7-genus-RAW.qza"

# Genus output directory.
GENUS_CORE_DIR="${GENUS_CORE_DIR}"

# Export and visualization directories.
GENUS_CORE_EXP="${GENUS_CORE_EXP}"
GENUS_CORE_QZV="${GENUS_CORE_QZV}"



########################################################
## Check required input files and directories.
########################################################

check_file "${GENUS_TABLE}" "Raw genus abundance table"
check_file "${METADATA}" "QIIME 2 sample metadata"



########################################################
## Run genus-level core metrics.
########################################################
#IMPORTANT:
# The genus table is derived by collapsing ASVs into genera.
# The original ASV phylogenetic tree is not directly suitable
# for genus-level metrics!!
#
#Therefore, use non-phylogenetic core metrics.
#
#Outputs include:
#   - Shannon diversity
#   - Observed genera
#   - Pielou's evenness
#   - Bray-Curtis distance matrix
#   - Jaccard distance matrix
#
#Faith's PD and UniFrac are not calculated here!!
#
#QIIME 2 requires the output directory not to exist so 
# previous genus results were removed in the initial setup.

log "Running non-phylogenetic genus-level core metrics."

qiime diversity core-metrics \
  --i-table "${GENUS_TABLE}" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${METADATA}" \
  --p-n-jobs "${THREADS}" \
  --output-dir "${GENUS_CORE_DIR}" \
  2>&1 | tee -a "${LOGFILE}"



########################################################
## Create directories for exports and visualizations.
########################################################
mkdir -p "${GENUS_CORE_EXP}"
mkdir -p "${GENUS_CORE_QZV}"


########################################################
## Export genus-level core metrics.
########################################################

### Alpha-diversity: Shannon
qiime tools export \
  --input-path "${GENUS_CORE_DIR}/shannon_vector.qza" \
  --output-path "${GENUS_CORE_EXP}/shannon" \
  2>&1 | tee -a "${LOGFILE}"

### Alpha-diversity: observed genera
qiime tools export \
  --input-path "${GENUS_CORE_DIR}/observed_features_vector.qza" \
  --output-path "${GENUS_CORE_EXP}/observed" \
  2>&1 | tee -a "${LOGFILE}"

### Alpha-diversity: Pielou's evenness
qiime tools export \
  --input-path "${GENUS_CORE_DIR}/evenness_vector.qza" \
  --output-path "${GENUS_CORE_EXP}/evenness" \
  2>&1 | tee -a "${LOGFILE}"

### Beta-diversity: Bray-Curtis
qiime tools export \
  --input-path "${GENUS_CORE_DIR}/bray_curtis_distance_matrix.qza" \
  --output-path "${GENUS_CORE_EXP}/bray" \
  2>&1 | tee -a "${LOGFILE}"

### Beta-diversity: Jaccard
qiime tools export \
  --input-path "${GENUS_CORE_DIR}/jaccard_distance_matrix.qza" \
  --output-path "${GENUS_CORE_EXP}/jaccard" \
  2>&1 | tee -a "${LOGFILE}"



########################################################
## Alpha-diversity group significance.
########################################################
#These visualizations compare alpha-diversity values across
# metadata groups. QIIME 2 can display the available metadata
# columns for exploratory comparisons.

### Shannon diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${GENUS_CORE_DIR}/shannon_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${GENUS_CORE_QZV}/shannon-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"

### Observed genera
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${GENUS_CORE_DIR}/observed_features_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${GENUS_CORE_QZV}/observed-features-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"

### Pielou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${GENUS_CORE_DIR}/evenness_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${GENUS_CORE_QZV}/evenness-group-significance.qzv" \
  2>&1 | tee -a "${LOGFILE}"



########################################################
## Beta-diversity group significance.
########################################################
#The example grouping variable is "sex".
#IMPORTANT!! Change it for other variables whn looking at other variables!!.
#
#PERMANOVA tests whether community composition differs
# among groups but its best to check group sizes and dispersion before
# interpreting the results.

### Bray-Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix "${GENUS_CORE_DIR}/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${GENUS_CORE_QZV}/bray-curtis-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"

### Jaccard
qiime diversity beta-group-significance \
  --i-distance-matrix "${GENUS_CORE_DIR}/jaccard_distance_matrix.qza" \
  --m-metadata-file "${METADATA}" \
  --m-metadata-column sex \
  --o-visualization "${GENUS_CORE_QZV}/jaccard-sex.qzv" \
  --p-pairwise \
  2>&1 | tee -a "${LOGFILE}"



########################################################
## ANCOM-BC2: genus-level differential abundance.
## Use raw, non-rarefied genus counts!
########################################################
log "Running genus-level ANCOM-BC2."

#Create a temporary directory for intermediate files.
ANCOM_GENUS_TMP="$(mktemp -d)"

#Filter genera by total abundance and prevalence.
qiime feature-table filter-features \
  --i-table "${GENUS_TABLE}" \
  --p-min-frequency 50 \
  --p-min-samples 2 \
  --o-filtered-table "${ANCOM_GENUS_TMP}/Genus_Filtered.qza" \
  2>&1 | tee -a "${LOGFILE}"

#Run ANCOM-BC2 using sex as the fixed effect.
qiime composition ancombc2 \
  --i-table "${ANCOM_GENUS_TMP}/Genus_Filtered.qza" \
  --m-metadata-file "${METADATA}" \
  --p-fixed-effects-formula 'sex' \
  --o-ancombc2-output "${ANCOM_GENUS_TMP}/ANCOMBC2_Genus_results.qza" \
  2>&1 | tee -a "${LOGFILE}"

#Save the results permanently in the genus output directory.
cp "${ANCOM_GENUS_TMP}/ANCOMBC2_Genus_results.qza" \
   "${GENUS_CORE_DIR}/ANCOMBC2_Genus_results.qza"

#Create the differential-abundance visualization.
qiime composition da-barplot \
  --i-data "${GENUS_CORE_DIR}/ANCOMBC2_Genus_results.qza" \
  --o-visualization "${GENUS_CORE_QZV}/ANCOMBC2_Genus_barplot.qzv" \
  2>&1 | tee -a "${LOGFILE}"

#Validate the outputs.
check_file "${GENUS_CORE_DIR}/ANCOMBC2_Genus_results.qza" \
  "Genus-level ANCOM-BC2 results"

check_file "${GENUS_CORE_QZV}/ANCOMBC2_Genus_barplot.qzv" \
  "Genus-level ANCOM-BC2 barplot"

#Remove temporary files only after successful completion.
rm -rf "${ANCOM_GENUS_TMP}"

log "Finished genus-level ANCOM-BC2."

log "FINISHED STEP 13: Core metrics for genus-level exploratory analysis."



########################################################
## Creating a QZA folder for better structure
########################################################
#Define separate QZA directories
ASV_CORE_QZA="${ASV_CORE_DIR}/Core_metrics_QZA"
GENUS_CORE_QZA="${GENUS_CORE_DIR}/Core_metrics_QZA"
mkdir -p "${ASV_CORE_QZA}" "${GENUS_CORE_QZA}" 


mv "${ASV_CORE_DIR}"/*.qza "${ASV_CORE_QZA}"
mv "${GENUS_CORE_DIR}"/*.qza "${GENUS_CORE_QZA}"


########################################################
## Moving some lost QZV emperor files
########################################################
#Some of the emperor QZV files are still in the MAIN core dir
# so this is just to maek sure it all ends clean. 
mv "${ASV_CORE_DIR}"/*.qzv "${ASV_CORE_QZV}"
mv "${GENUS_CORE_DIR}"/*.qzv "${GENUS_CORE_QZV}"



########################################################
## Summary
########################################################
log "==========================================================================================================="
log "Script 7_Core-metrics.sh finished running!"
log " "
log "CORE METRICS AND DIVERSITY ANALYSIS COMPLETE."
log " "
log "ASV-LEVEL RESULTS:"
log "  Core metrics artifacts:"
log "    - ${ASV_CORE_DIR}/"
log "      (ASV-level rarefied diversity metrics and distance matrices)"
log " "
log "  Exported results:"
log "    - ${ASV_CORE_EXP}/"
log " "
log "  Visualizations:"
log "    - ${ASV_CORE_QZV}/"
log " "
log "  ANCOM-BC2 results, if generated in STEP 12:"
log "    - ${ASV_CORE_DIR}/ANCOMBC2_ASV_results.qza"
log "    - ${ASV_CORE_QZV}/ANCOMBC2_ASV_barplot.qzv"
log " "
log "GENUS-LEVEL RESULTS:"
log "  Core metrics artifacts:"
log "    - ${GENUS_CORE_DIR}/"
log "      (Shannon, observed genera, evenness, Bray-Curtis and Jaccard)"
log " "
log "  Exported results:"
log "    - ${GENUS_CORE_EXP}/"
log " "
log "  Visualizations:"
log "    - ${GENUS_CORE_QZV}/"
log " "
log "==========================================================================================================="
log "NEXT STEPS:"
log "  1. Inspect the .qzv visualizations at https://view.qiime2.org."
log "  2. Check the sampling depth and confirm that enough samples were retained."
log "  3. Check the group sizes in the metadata before interpreting significance tests."
log "  4. Use raw, non-rarefied counts for appropriate differential-abundance analyses."
log "  5. Continue with the comprehensive downstream R-based statistical analyses."
log "==========================================================================================================="