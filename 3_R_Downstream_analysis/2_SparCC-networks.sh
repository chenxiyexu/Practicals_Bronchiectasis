#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob #Modify filename expansion (globbing). 
#
#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 
#
#This part of the downstream analysis is for the correlation network analysis using SparCC/Fastspar.
#Fastspar in a fast implementation of the SparCC algorithm, which is designed to estimate correlations 
# in compositional data (like microbiome data).

#Basic variables. 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Global_functions_DOWNSTREAM.sh"
trap "error_check ${LINENO}" ERR





#########################################################
### FastSpar code for microbiome networks analysis!
#########################################################

################# Files needed for the analysis.
QIIME_DIR="$SCRIPT_DIR/../2_Qiime2"

#Non-rarefied genus abundance tables.
#We use non-rarefied data because this is a compositional analysis.

#Data without sex separation.
INPUT_WITHOUT="${QIIME_DIR}/Final_results/FINAL_tables/RAW/GENUS_table_RAW.tsv"

#Data separated by patient sex.
INPUT_FEMALE="${QIIME_DIR}/Final_results/FINAL_tables/SEX_separated/FEMALE_only/GENUS_table_RAW_FEMALE.tsv"

INPUT_MALE="${QIIME_DIR}/Final_results/FINAL_tables/SEX_separated/MALE_only/GENUS_table_RAW_MALE.tsv"


################# Files created as a result of this analysis.
NETWORKS_DIR="${SCRIPT_DIR}/Correlation_networks"
BASE_DIR="${NETWORKS_DIR}/SparCC_results"

#Main output directories.
WITHOUT_DIR="${BASE_DIR}/WITHOUT_separation"
SEX_DIR="${BASE_DIR}/SEX_separation"

FEMALE_DIR="${SEX_DIR}/FEMALE_only"
MALE_DIR="${SEX_DIR}/MALE_only"

mkdir -p "$WITHOUT_DIR" "$FEMALE_DIR" "$MALE_DIR"



################# Other variable configurations.
PERMUTATIONS="${PERMUTATIONS:-1000}"
ITERATIONS="${ITERATIONS:-500}"
THREADS="${THREADS:-1}"

#Number of bootstraps to run in parallel.
#Each individual bootstrap correlation uses --threads 1.
JOBS="${JOBS:-$THREADS}"

log "Starting downstream taxa correlation analysis using FastSpar!"



#########################################################
### STEP 0: Checking requirements
#########################################################
log "STEP 0: Checking requirements."

#Verify that the correct conda environment is activated.
for bin in fastspar fastspar_bootstrap fastspar_pvalues; do
    if ! command -v "$bin" >/dev/null 2>&1; then
        echo "ERROR: Could not find '$bin' in the PATH." >&2
        echo "Please activate the correct environment or install FastSpar!" >&2
        exit 1111
    fi
done


#Check that all input abundance tables exist and are not empty.
for INPUT in "$INPUT_WITHOUT" "$INPUT_FEMALE" "$INPUT_MALE"; do
    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: Input abundance table is missing or empty: $INPUT" >&2
        exit 1111
    fi
done

log "All three input abundance tables are available."

log "FINISHED STEP 0: Checking requirements."


#########################################################
### FUNCTION: run_fastspar
###
###This function independently performs the complete
### FastSpar analysis for one abundance table.
###
###Arguments:
###  1. Input abundance table
###  2. Output analysis directory
###
###The function generates:
###   - Original abundance matrix
###   - Filtered abundance matrix
###   - Correlation matrix
###   - Covariance matrix
###   - Bootstrap abundance tables
###   - Bootstrap correlation matrices
###   - Bootstrap covariance matrices
###   - Permutation-based p-values
#########################################################

run_fastspar() {
    local INPUT_TSV="$1"
    local ANALYSIS_DIR="$2"

    #####################################################
    ## Output directories
    #####################################################
    local ABUNDANCE_MATRIX_DIR="${ANALYSIS_DIR}/Matrices"

    local FILTERED_TSV="${ABUNDANCE_MATRIX_DIR}/abundance_matrix_genus_filtered.tsv"

    local BOOT_DIR="${ANALYSIS_DIR}/Bootstraps"
    local BOOT_CORR_DIR="${ANALYSIS_DIR}/Bootstrap_correlations"

    local CORR_OUT="${ANALYSIS_DIR}/correlation.tsv"
    local COV_OUT="${ANALYSIS_DIR}/covariance.tsv"
    local PVAL_OUT="${ANALYSIS_DIR}/pvalues.tsv"

    mkdir -p "$ABUNDANCE_MATRIX_DIR" "$BOOT_DIR" "$BOOT_CORR_DIR"


    #####################################################
    ## STEP 0: Checking input and preparing directories
    #####################################################
    log "Preparing FastSpar analysis: ${ANALYSIS_DIR}"

    if [[ ! -s "$INPUT_TSV" ]]; then
        echo "ERROR: Input abundance table missing or empty: $INPUT_TSV" >&2
        return 1111
    fi

    #Copy original abundance matrix for convenience.
    local ORIGINAL_TSV="${ABUNDANCE_MATRIX_DIR}/OG_GENUS_table_RAW.tsv"

    cp "$INPUT_TSV" "$ORIGINAL_TSV"

    #Use the copied file for the remaining analysis.
    INPUT_TSV="$ORIGINAL_TSV"

    log "Original abundance matrix copied to: ${ORIGINAL_TSV}"


    #####################################################
    ### STEP 1: Filter low abundance taxa
    #####################################################
    log "STEP 1: Filtering low abundance taxa."

    #This filtering approach checks that general prevalence
    # is above 20% of the samples.

    awk 'BEGIN { FS=OFS="\t" }

    #Skip BIOM comment line.
    NR==1 && $0 ~ /^# Constructed from biom file/ {
        next
    }

    #Process the actual header.
    !header_seen && $1=="#OTU ID" {

        print

        n_samples=NF-1

        min_prevalence=int(n_samples*0.20)

        if (n_samples*0.20 > min_prevalence)
            min_prevalence++

        header_seen=1

        next
    }

    #Filter taxa by prevalence.
    header_seen {

        prevalence=0

        for (i=2; i<=NF; i++) {

            if (($i+0)>0)
                prevalence++

        }

        if (prevalence>=min_prevalence)
            print
    }

    END {

        if (!header_seen) {

            print "ERROR: #OTU ID header not found" > "/dev/stderr"

            exit 1
        }

    }' "$INPUT_TSV" > "$FILTERED_TSV"


    if [[ $? -ne 0 ]]; then
        echo "ERROR: Filtering failed for $INPUT_TSV" >&2
        return 1111
    fi


    #Count retained taxa.
    local N_TAXA

    N_TAXA=$(awk -F '\t' '
    $1=="#OTU ID" {next}
    /^# Constructed from biom file/ {next}
    NF>0 {n++}
    END {print n+0}
    ' "$FILTERED_TSV")


    if [[ "$N_TAXA" -lt 2 ]]; then
        echo "ERROR: Only $N_TAXA taxa remain after filtering." >&2
        return 1111
    fi

    log "Number of taxa after filtering: ${N_TAXA}"

    log "FINISHED STEP 1: Filtering low abundance taxa."


    #####################################################
    ### STEP 2: Main correlations using SparCC/FastSpar
    #####################################################
    log "STEP 2: Main correlations using SparCC/FastSpar."

    fastspar \
        --otu_table "$FILTERED_TSV" \
        --correlation "$CORR_OUT" \
        --covariance "$COV_OUT" \
        --iterations "$ITERATIONS" \
        --threads "$THREADS"

    if [[ ! -s "$CORR_OUT" || ! -s "$COV_OUT" ]]; then
        echo "ERROR: FastSpar failed to generate correlation/covariance matrices." >&2
        return 1111
    fi

    log "FINISHED STEP 2: Main correlations using SparCC/FastSpar."



    #####################################################
    ### STEP 3: Generate bootstrap tables
    #####################################################
    log "STEP 3: Generate bootstrap tables using FastSpar."

    fastspar_bootstrap \
        --otu_table "$FILTERED_TSV" \
        --number "$PERMUTATIONS" \
        --prefix "$BOOT_DIR/bootstrap" \
        --threads "$THREADS"

    shopt -s nullglob

    local boot_files=("$BOOT_DIR"/bootstrap_*.tsv)
    shopt -u nullglob


    if [[ "${#boot_files[@]}" -eq 0 ]]; then
        echo "ERROR: FastSpar did not generate bootstrap files." >&2
        return 1111
    fi


    if [[ "${#boot_files[@]}" -ne "$PERMUTATIONS" ]]; then
        echo "WARNING: Expected $PERMUTATIONS bootstraps, but found ${#boot_files[@]}!" >&2
    fi

    log "Generated ${#boot_files[@]} bootstrap abundance tables."

    log "FINISHED STEP 3: Generate bootstrap tables."



    #####################################################
    ### STEP 4: Bootstrap correlations
    ###
    ### Each bootstrap is analysed independently.
    ### Up to JOBS analyses can run simultaneously.
    #####################################################
    log "STEP 4: Correlations using bootstraps."

    run_one_bootstrap() {

        local f="$1"
        local base

        base=$(basename "$f")
        fastspar \
            --otu_table "$f" \
            --correlation "$BOOT_CORR_DIR/corr_${base}" \
            --covariance "$BOOT_CORR_DIR/cov_${base}" \
            --iterations "$ITERATIONS" \
            --threads 1
    }


    export -f run_one_bootstrap

    export BOOT_CORR_DIR ITERATIONS


    if ! printf '%s\n' "${boot_files[@]}" |
        xargs -P "$JOBS" -I{} bash -c 'run_one_bootstrap "$@"' _ {}; then
        echo "ERROR: One or more bootstrap correlation jobs failed!" >&2
        return 1111
    fi

    local n_cor

    n_cor=$(find "$BOOT_CORR_DIR" -maxdepth 1 -type f \
        -name 'corr_bootstrap_*.tsv' | wc -l)
    if [[ "$n_cor" -ne "${#boot_files[@]}" ]]; then
        echo "ERROR: Generated $n_cor correlation files but expected ${#boot_files[@]}!" >&2
        return 1111
    fi

    log "Generated ${n_cor} bootstrap correlation matrices."

    log "FINISHED STEP 4: Correlations using bootstraps."


    #####################################################
    ### STEP 5: Obtaining p-values
    #####################################################
    log "STEP 5: Obtaining p-values."

    if ! fastspar_pvalues \
        --otu_table "$FILTERED_TSV" \
        --correlation "$CORR_OUT" \
        --prefix "$BOOT_CORR_DIR/corr_bootstrap_" \
        --permutations "$PERMUTATIONS" \
        --outfile "$PVAL_OUT"; then

        echo "ERROR: FastSpar failed to calculate p-values!" >&2
        return 1111
    fi


    if [[ ! -s "$PVAL_OUT" ]]; then
        echo "ERROR: P-value output file is missing or empty: $PVAL_OUT" >&2
        return 1111
    fi

    log "FINISHED STEP 5: Obtaining p-values."


    #####################################################
    ### STEP 6: Analysis summary
    #####################################################
    log "FastSpar analysis completed successfully: ${ANALYSIS_DIR}"

    log "Filtered matrix: ${FILTERED_TSV}"
    log "Correlations: ${CORR_OUT}"
    log "Covariances: ${COV_OUT}"
    log "P-values: ${PVAL_OUT}"

    log "Number of taxa retained: ${N_TAXA}"
    log "Number of bootstraps: ${#boot_files[@]}"

    return 0
}



#########################################################
### STEP 1: FastSpar WITHOUT sex separation
#########################################################
log "=========================================================="
log "STARTING ANALYSIS WITHOUT SEX SEPARATION"
log "=========================================================="

run_fastspar \
    "$INPUT_WITHOUT" \
    "$WITHOUT_DIR" || exit 1111

log "FINISHED ANALYSIS WITHOUT SEX SEPARATION."



#########################################################
### STEP 2: FastSpar FEMALE only
#########################################################
log "=========================================================="
log "STARTING FEMALE-ONLY ANALYSIS"
log "=========================================================="

run_fastspar \
    "$INPUT_FEMALE" \
    "$FEMALE_DIR" || exit 1111

log "FINISHED FEMALE-ONLY ANALYSIS."



#########################################################
### STEP 3: FastSpar MALE only
#########################################################
log "=========================================================="
log "STARTING MALE-ONLY ANALYSIS"
log "=========================================================="

run_fastspar \
    "$INPUT_MALE" \
    "$MALE_DIR" || exit 1111

log "FINISHED MALE-ONLY ANALYSIS."



#########################################################
### FINAL SUMMARY
#########################################################
log "==========================================================================================================="
log "FASTSPAR-BASED GENUS-LEVEL CORRELATION NETWORK ANALYSIS COMPLETE."
log ""
log "ANALYSES PERFORMED:"
log ""
log "  1. WITHOUT SEX SEPARATION:"
log "     - ${WITHOUT_DIR}/"
log ""
log "  2. FEMALE-ONLY:"
log "     - ${FEMALE_DIR}/"
log ""
log "  3. MALE-ONLY:"
log "     - ${MALE_DIR}/"
log ""
log "EACH ANALYSIS CONTAINS:"
log "  - Matrices/OG_GENUS_table_RAW.tsv"
log "  - Matrices/abundance_matrix_genus_filtered.tsv"
log "  - Bootstraps/bootstrap_*.tsv"
log "  - Bootstrap_correlations/corr_bootstrap_*.tsv"
log "  - Bootstrap_correlations/cov_bootstrap_*.tsv"
log "  - correlation.tsv"
log "  - covariance.tsv"
log "  - pvalues.tsv"
log ""
log "ANALYSIS CONFIGURATION:"
log "  - Number of bootstrap/permutation tables: ${PERMUTATIONS}"
log "  - FastSpar iterations per correlation analysis: ${ITERATIONS}"
log "  - FastSpar threads per main analysis: ${THREADS}"
log "  - Maximum parallel bootstrap correlation jobs: ${JOBS}"
log ""
log "==========================================================================================================="
log "NEXT STEPS:"
log "  1. Inspect the filtered abundance matrices to confirm that the retained taxa"
log "     and sample columns are correct."
log "  2. Compare correlation coefficients between the three independent analyses."
log "  3. Use the corresponding pvalues.tsv files to assess statistical significance."
log "  4. Apply an appropriate multiple-testing correction (e.g., Benjamini-Hochberg FDR)."
log "  5. Prepare the resulting networks for downstream analysis and visualization in R."
log "==========================================================================================================="
log "All FastSpar downstream correlation network analyses completed successfully."
log "==========================================================================================================="
