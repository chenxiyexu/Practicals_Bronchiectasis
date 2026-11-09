#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob #Modify filename expansion (globbing). 
#
#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 
#
#Special thanks to Pneumobiome from IRYCIS for their support and for allowing me to 
# research such an interesting topic, and above all special thanks to Sandra Mingo Ramírez for 
# letting me use her scripts as a base and for her help in the development of this pipeline.
#
#This part of the downstream analysis is for the correlation network analysis using SparCC/Fastspar.
#Fastspar in a fast implementation of the SparCC algorithm, which is designed to estimate correlations 
# in compositional data (like microbiome data).

#Basic variables. 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Global_functions_DOWNSTREAM.sh"
trap "error_check ${LINENO}" ERR




#########################################################
### FastSpar code for microbiome networkds analysis!
#########################################################


#################Files needed for the analysis. 
QIIME_DIR="$SCRIPT_DIR/../2_Qiime2"

#Inside the Qiime folder find the NON-RAREFIED ABUNDANCE TABLE .tsv for genus!
#We use non-rarefied data because this is a compositional analysis!
INPUT_TSV="$QIIME_DIR/Final_results/FINAL_tables/RAW/GENUS_table_RAW.tsv"


#################Files created as a result of this analysis. 
BASE_DIR="$SCRIPT_DIR/Correlation_networks"
mkdir -p "$BASE_DIR"

ABUNDANCE_MATRIX_DIR="$BASE_DIR/Matrices"
mkdir -p "$ABUNDANCE_MATRIX_DIR"

FILTERED_TSV="$ABUNDANCE_MATRIX_DIR/abundance_matrix_genus_filtered.tsv"

BOOT_DIR="$BASE_DIR/Bootstraps"
BOOT_CORR_DIR="$BASE_DIR/Bootstrap_correlations"
mkdir -p "$BOOT_DIR" "$BOOT_CORR_DIR"

CORR_OUT="$BASE_DIR/correlation.tsv"
COV_OUT="$BASE_DIR/covariance.tsv"
PVAL_OUT="$BASE_DIR/pvalues.tsv"


#################Other variable configurations. 
PERMUTATIONS="${PERMUTATIONS:-1000}"
ITERATIONS="${ITERATIONS:-500}"
THREADS="${THREADS:-1}"


#Number of bootstraps to run in parallel (independent of --threads of fastspar,
# which controls the internal threads of EACH call). 
JOBS="${JOBS:-$THREADS}"

log "Starting downstream taxa correlation analysis using FastSpar!"


#########################################################
### STEP 0: Checking requirements
#########################################################
log "STEP 0: Checking requirements."

#Verify the correct conda env is activated (which includes Fastspar).
for bin in fastspar fastspar_bootstrap fastspar_pvalues; do
    if ! command -v "$bin" >/dev/null 2>&1; then
        echo "ERROR: could not find '$bin' in the PATH." >&2
        echo "Please activate the correct environment or download Fastspar!" >&2
        exit 1111
    fi
done

#Just checking if file exists (wth -f)! 
if [[ ! -f "$INPUT_TSV" ]]; then
    echo "ERROR: input raw/non-rarefied abundance table does not exist at: $INPUT_TSV !" >&2
    exit 1111
fi

#Check that the input file is not empty (with -s).
if [[ ! -s "$INPUT_TSV" ]]; then
    echo "ERROR: input abundance table is empty: $INPUT_TSV" >&2
    exit 1111
fi

#Creating output directories.
echo "Copying original genus raw/non-rarefied abundance matrix into $ABUNDANCE_MATRIX_DIR for convenience!"
cp "$INPUT_TSV" "$ABUNDANCE_MATRIX_DIR/OG_GENUS_table_RAW.tsv"

#Changing to the file copied at the closer directory (in case we need to use or change this file 
# for whatever reason)!
INPUT_TSV="$ABUNDANCE_MATRIX_DIR/OG_GENUS_table_RAW.tsv"

log "FINISHED STEP 0: Checking requirements."


#########################################################
### STEP 1: Filter low abundance taxa
#########################################################
#This is done in order to avoid problems with SparCC/FastSpar, which is not designed 
# to handle taxa that are absent in most samples.

log "STEP 1: Filtering low abundance taxa."

: ' comentario! 
This filter checks that total abundance is > 1. 
Removed single comillas before BEGIN and afer } for
multiline comment.

awk BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
s=0
for(i=2;i<=NF;i++) s+=($i+0)
if(s>1) print
} "$INPUT_TSV" > "$FILTERED_TSV"
'

#This filtering approach checks that general prevalence is above 10% of the samples. 
awk 'BEGIN { FS=OFS="\t" }

#Skip the BIOM comment line
NR==1 && $0 ~ /^# Constructed from biom file/ { next }

#Process the actual header
!header_seen && $1=="#OTU ID" {
    print
    n_samples=NF-1
    min_prevalence=int(n_samples*0.10)
    if (n_samples*0.10 > min_prevalence) min_prevalence++
    header_seen=1
    next
}

#Filter taxa by prevalence
header_seen {
    prevalence=0
    for (i=2; i<=NF; i++) {
        if (($i+0)>0) prevalence++
    }
    if (prevalence>=min_prevalence) print
}

END {
    if (!header_seen) {
        print "ERROR: #OTU ID header not found" > "/dev/stderr"
        exit 1
    }
}' "$INPUT_TSV" > "$FILTERED_TSV"

if [[ $? -ne 0 ]]; then
    echo "ERROR: filtering failed." >&2
    exit 1
fi

N_TAXA=$(($(wc -l < "$FILTERED_TSV") - 1))

if [[ "$N_TAXA" -lt 2 ]]; then
    echo "ERROR: after filtering low abundance taxa there are only $N_TAXA taxa left (at least 2 or more are needed!) in $FILTERED_TSV" >&2
    exit 1111
fi

echo "Number of taxa after filtering: $N_TAXA"

log "FINISHED STEP 1: Filter low abundance taxa."


#########################################################
### STEP 2: Main correlations using SparCC/Fastspar
#########################################################
log "STEP 2: Main correlations using SparCC/Fastspar."

#Correlation is in Pearson correlation from variance of log
# ratios. It ranges from -1 to 1. 
fastspar \
--otu_table "$FILTERED_TSV" \
--correlation "$CORR_OUT" \
--covariance "$COV_OUT" \
--iterations "$ITERATIONS" \
--threads "$THREADS"

log "FINISHED STEP 2: Main correlations using SparCC/Fastspar."


#########################################################
### STEP 3: Generate bootstrap tables using SparCC/Fastspar
#########################################################
log "STEP 3: Generate bootstrap tables using SparCC/Fastspar."

fastspar_bootstrap \
--otu_table "$FILTERED_TSV" \
--number "$PERMUTATIONS" \
--prefix "$BOOT_DIR/bootstrap" \
--threads "$THREADS"

shopt -s nullglob
boot_files=("$BOOT_DIR"/bootstrap_*.tsv)
shopt -u nullglob

if [[ "${#boot_files[@]}" -eq 0 ]]; then
    echo "ERROR: fastspar_bootstrap did NOT generate any files at: $BOOT_DIR" >&2
    exit 1111
fi

if [[ "${#boot_files[@]}" -ne "$PERMUTATIONS" ]]; then
    echo "WARNING: we expected $PERMUTATIONS bootstraps, but found ${#boot_files[@]}!" >&2
fi

log "FINISHED STEP 3: Generate bootstrap tables using SparCC/Fastspar."


#########################################################
### STEP 4: Correlations using bootstraps (in parallel, 
###         each with --threads 1 to avoid overloading 
###         CPUs; up to $JOBS simultaneously).
#########################################################
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
    echo "ERROR: one or more bootstrap correlation jobs failed!" >&2
    exit 1
fi

n_cor=$(find "$BOOT_CORR_DIR" -maxdepth 1 -type f \
    -name 'corr_bootstrap_*.tsv' | wc -l)

if [[ "$n_cor" -ne "$PERMUTATIONS" ]]; then
    echo "ERROR: we generated $n_cor correlation bootstrap files but expected $PERMUTATIONS!" >&2
    exit 1
fi

log "FINISHED STEP 4: Correlations using bootstraps."


#########################################################
### STEP 5: Obtaining p-values
#########################################################
#Fastspar uses a permutation based approach and the resulting p-value 
# is how often a more extreme correlation is observed for randomly 
# permutated data. This ensures it is robust and not overestimated. 

log "STEP 5: Obtaining p-values."

if ! fastspar_pvalues \
    --otu_table "$FILTERED_TSV" \
    --correlation "$CORR_OUT" \
    --prefix "$BOOT_CORR_DIR/corr_bootstrap_" \
    --permutations "$PERMUTATIONS" \
    --outfile "$PVAL_OUT"; then
    echo "ERROR: FastSpar failed to calculate p-values!" >&2
    exit 1
fi

if [[ ! -s "$PVAL_OUT" ]]; then
    echo "ERROR: p-value output file is missing or empty: $PVAL_OUT" >&2
    exit 1
fi

echo "Pipeline completed."
echo "Filtered matrix:  $FILTERED_TSV"
echo "Correlations:     $CORR_OUT"
echo "Covariances:      $COV_OUT"
echo "P-values:         $PVAL_OUT"

log "FINISHED STEP 5: Obtaining p-values."
log "FINISHED Downstream taxa correlation analysis using FastSpar!"


#########################################################
## Final logs
#########################################################
log "==========================================================================================================="
log "FASTSPAR-BASED GENUS-LEVEL CORRELATION NETWORK ANALYSIS COMPLETE."
log "Key outputs:"
log "  Input abundance matrix:"
log "    - ${ABUNDANCE_MATRIX_DIR}/OG_GENUS_table_RAW.tsv"
log "      (Original genus-level raw/non-rarefied abundance table)"
log ""
log "  Filtered abundance matrix:"
log "    - ${FILTERED_TSV}"
log "      (Taxa retained after filtering rows with a total abundance <= 1)"
log ""
log "  Main correlation and covariance matrices:"
log "    - ${CORR_OUT}"
log "      (Genus-level correlation coefficients estimated using FastSpar)"
log "    - ${COV_OUT}"
log "      (Genus-level covariance estimates produced by FastSpar)"
log ""
log "  Bootstrap outputs:"
log "    - ${BOOT_DIR}/bootstrap_*.tsv"
log "      (Bootstrap abundance tables generated from the filtered matrix)"
log "    - ${BOOT_CORR_DIR}/corr_bootstrap_*.tsv"
log "      (Correlation matrices estimated from bootstrap abundance tables)"
log "    - ${BOOT_CORR_DIR}/cov_bootstrap_*.tsv"
log "      (Covariance matrices estimated from bootstrap abundance tables)"
log ""
log "  Statistical significance:"
log "    - ${PVAL_OUT}"
log "      (Permutation-based p-values for the estimated correlations)"
log ""
log "  Analysis configuration:"
log "    - Number of bootstrap/permutation tables: ${PERMUTATIONS}"
log "    - FastSpar iterations per correlation analysis: ${ITERATIONS}"
log "    - FastSpar threads per main analysis: ${THREADS}"
log "    - Maximum parallel bootstrap correlation jobs: ${JOBS}"
log "    - Number of taxa retained after filtering: ${N_TAXA}"
log "==========================================================================================================="
log "NEXT STEPS:"
log "  1. Inspect the filtered abundance matrix to confirm that the retained taxa"
log "     and sample columns are correct."
log "  2. Examine ${CORR_OUT} to identify positive and negative genus-level associations."
log "  3. Use ${PVAL_OUT} to assess the statistical significance of the estimated correlations."
log "  4. Apply an appropriate multiple-testing correction (e.g., Benjamini-Hochberg FDR)"
log "     before interpreting statistically significant associations."
log "  5. Combine correlation coefficients and adjusted p-values to construct and visualize"
log "     the genus-level microbial association network."
log "==========================================================================================================="
log "All FastSpar downstream correlation network analysis steps completed successfully."
log "==========================================================================================================="