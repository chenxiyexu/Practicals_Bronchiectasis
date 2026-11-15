#!/usr/bin/env bash

#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#We create a folder for the logs.
mkdir -p "${SCRIPT_DIR}/Qiime2_logs"
LOG_DIR="${SCRIPT_DIR}/Qiime2_logs"
LOGFILE="${LOG_DIR}/Qiime2_pipeline_log.txt"


#########################################################
### This log function writes down the time of each step
### 	to the file Qiime2_pipeline_log.txt. 
#########################################################
function log {
    local message="$1"
    local timestamp

    timestamp="$(date "+%d-%m-%Y %H:%M:%S")"
    
    echo "[${timestamp}] ${message}" | tee -a "$LOGFILE"
}


#########################################################
### This error_check function is called anytime a command
### 	fails. 	
#########################################################

#BASH_SOURCE[0] is the current file/function where this code is executing.
#BASH_SOURCE[1] the caller file or function. 

function error_check {
    local exit_code="$1"
    local line_number="$2"
    local failed_command="$3"
    local script="${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}"

    log "ERROR! Command failed with exit code ${exit_code}"
    log "Script: ${script}"
    log "Line: ${line_number}"
    log "Command: ${failed_command}"
    log "Exiting."

    exit 1111
}


#########################################################
### This check_file function is called to verify a file
### 	was created alright. 	
#########################################################
function check_file {
    local file="$1"
    local description="$2"

    if [[ -f "$file" ]]; then
        log "SUCCESS! ${description}: ${file}"
    else
        log "ERROR! ${file} was not created!"
        exit 1111
    fi
}



#########################################################
## This split_by_sex function separates samples using
##      sex as a variable to study the networks independently.
#########################################################
function split_by_sex() {

    local INPUT_TABLE="$1"
    local METADATA_TSV="$2"
    local SEX_VALUE="$3"
    local OUTPUT_TABLE="$4"

    local TEMP_OUTPUT="${OUTPUT_TABLE}.tmp"

    log "Separating samples: ${SEX_VALUE}"


    #####################################################
    ## Validate arguments and input files
    #####################################################
    if [[ "$#" -ne 4 ]]; then
        echo "ERROR: split_by_sex requires exactly 4 arguments." >&2
        return 1111
    fi

    if [[ ! -s "$INPUT_TABLE" ]]; then
        echo "ERROR: Input table missing or empty: $INPUT_TABLE" >&2
        return 1111
    fi

    if [[ ! -s "$METADATA_TSV" ]]; then
        echo "ERROR: Metadata missing or empty: $METADATA_TSV" >&2
        return 1111
    fi

    if [[ -z "$SEX_VALUE" || -z "$OUTPUT_TABLE" ]]; then
        echo "ERROR: Sex value or output path is empty." >&2
        return 1111
    fi

    mkdir -p "$(dirname "$OUTPUT_TABLE")" || return 1111


    #####################################################
    ## Separate samples using metadata
    #####################################################

    if ! awk -F '\t' -v sex_value="${SEX_VALUE}" -v metadata_file="${METADATA_TSV}" '

    BEGIN {
        OFS = "\t"
        sex_value = tolower(sex_value)
    }

    #####################################################
    ## PART 1: Read metadata
    #####################################################

    FILENAME == metadata_file {

        for (i=1; i<=NF; i++)
            gsub(/\r$/, "", $i)

        if (FNR == 1) {

            for (i=1; i<=NF; i++) {
                column = tolower($i)

                if (column == "sample-id" || column == "#sampleid" || column == "sampleid")
                    id_col = i

                if (column == "sex")
                    sex_col = i
            }

            if (!id_col || !sex_col) {
                print "ERROR: SampleID or sex column not found." \
                    > "/dev/stderr"
                failed = 1
                exit 1111
            }

            next
        }

        sample_id = $id_col
        sample_sex = tolower($sex_col)

        if (sample_id == "")
            next

        if (sample_id in metadata_ids) {
            print "ERROR: Duplicate sample ID in metadata: " sample_id \
                > "/dev/stderr"
            failed = 1
            exit 1111
        }

        metadata_ids[sample_id] = 1

        if (sample_sex == "") {
            missing_sex++
            next
        }

        if (sample_sex == sex_value) {
            selected[sample_id] = 1
            metadata_count++
        }

        next
    }


    #####################################################
    ## PART 2: Process abundance table
    #####################################################
    FILENAME != metadata_file {

        for (i=1; i<=NF; i++)
            gsub(/\r$/, "", $i)

        # Preserve BIOM comment
        if ($0 ~ /^# Constructed from biom file/) {
            print
            next
        }

        # Read abundance header
        if (!header_seen && $1 == "#OTU ID") {

            printf "%s", $1

            n_selected = 0
            n_unmatched = 0

            for (i=2; i<=NF; i++) {

                if ($i in table_ids) {
                    print "ERROR: Duplicate sample ID in abundance table: " $i \
                        > "/dev/stderr"
                    failed = 1
                    exit 1111
                }

                table_ids[$i] = 1

                if ($i in selected) {

                    columns[++n_selected] = i

                    matched[$i] = 1

                    printf "\t%s", $i

                } else {

                    if (!($i in metadata_ids)) {
                        print "WARNING: Sample absent from metadata: " $i \
                            > "/dev/stderr"
                        n_unmatched++
                    }
                }
            }

            print ""

            header_seen = 1

            if (n_selected == 0) {
                print "ERROR: No matching samples found." \
                    > "/dev/stderr"
                failed = 1
                exit 1111
            }

            next
        }

        # Extract taxa
        if (header_seen) {

            printf "%s", $1

            for (j=1; j<=n_selected; j++)
                printf "\t%s", $(columns[j])

            print ""

            taxa_count++
        }
    }


    #####################################################
    ## PART 3: Final validation
    #####################################################
    END {

        if (failed)
            exit 1111

        if (!header_seen) {
            print "ERROR: #OTU ID header not found." \
                > "/dev/stderr"
            exit 1
        }

        if (taxa_count == 0) {
            print "ERROR: No taxa found." \
                > "/dev/stderr"
            exit 1111
        }

        # Check selected metadata samples absent from table
        for (sample in selected) {

            if (!(sample in matched)) {
                print "WARNING: Selected metadata sample absent from abundance table: " sample \
                    > "/dev/stderr"
            }
        }

        printf "Sex: %s | Samples: %d | Taxa: %d\n", \
            sex_value, n_selected, taxa_count > "/dev/stderr"

        printf "Samples missing from metadata: %d\n", \
            n_unmatched > "/dev/stderr"

        printf "Samples without sex information: %d\n", \
            missing_sex > "/dev/stderr"
    }

    ' "$METADATA_TSV" "$INPUT_TABLE" > "$TEMP_OUTPUT"; then

        rm -f "$TEMP_OUTPUT"

        echo "ERROR: Failed to separate ${SEX_VALUE} samples." >&2

        return 1111
    fi


    #####################################################
    ## PART 4: Validate and save output
    #####################################################
    if [[ ! -s "$TEMP_OUTPUT" ]]; then
        rm -f "$TEMP_OUTPUT"
        echo "ERROR: Generated output is empty." >&2
        return 1111
    fi

    if ! mv "$TEMP_OUTPUT" "$OUTPUT_TABLE"; then
        rm -f "$TEMP_OUTPUT"
        echo "ERROR: Could not save output: $OUTPUT_TABLE" >&2
        return 1111
    fi

    log "Successfully generated: ${OUTPUT_TABLE}"

    return 0
}
