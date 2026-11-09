#!/usr/bin/env bash

#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#We create a folder for the logs.
mkdir -p "${SCRIPT_DIR}/Downstream_logs"
LOG_DIR="${SCRIPT_DIR}/Downstream_logs"
LOGFILE="${LOG_DIR}/Downstream_log.txt"


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
