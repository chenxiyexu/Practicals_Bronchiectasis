#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 


#Import all used libraries.
import os
import subprocess
import sys
from datetime import datetime
import pandas as pd



#######################################################################
## AUTOMATED SEQUENCE DOWNLOAD
#######################################################################

#This first script is used for automated sequences download from NCBI using sratools. 

#The only neccesary thing/file is a .csv with the accesion SRA code and the sample IDS 
# (later used to create the manifest.tsv that Qiime2 requires to import data).

#Download all SRA sequences from a list.
SRA_filename = input(f"Please input the name of the file with the SRA metadata (.csv) "
                      "(with both run code and sample IDs!): \n")
print("\n")

if not os.path.exists(SRA_filename):
    print(f"ERROR: The file {SRA_filename} was NOT found! Please try again.\n")
    sys.exit()


#Read and check SRA run and sample-id columns.
run_column = input(f"""Please input the column name of the SRA run codes: \n""")
print("\n")
sample_id_column = input(f"""Please input the column name of the sample IDs: \n""")
print("\n")


try:
    data = pd.read_csv(SRA_filename, sep = ",", header = 0)
    run_list = data[run_column].tolist()
    sample_id_list = data[sample_id_column].tolist()

except Exception as e:
    print(f"ERROR: {e}")
    sys.exit()
    
#Summary of the SRA codes file.
if len(run_list) == 0:
    print(f"No codes were found inside the provided list! Exiting the program.")
    sys.exit()
else:
    print(f"A total of {len(run_list)} SRA accession codes were found!\n"
        "Proceeding to download!"
        )


#We create a folder names "Sequences" and download all sequences inside it for a cleaner
# result. 
output_dir = "Sequences"
os.makedirs(output_dir, exist_ok = True)


logs = []
#We validate each SRA code and download it.
for i, sra in enumerate(run_list, 1):
    print(f"{i}: {sra} is being downloaded...")

    cmd = ["fasterq-dump", sra, "--outdir", output_dir]
    check_SRA_validity = subprocess.run(cmd, capture_output = True, text = True )

    if check_SRA_validity.returncode == 0:
        print(f"{i}: {sra} was successfully downloaded!\n")
        logs.append((i, sra, "SUCCESS"))
    else:
        print(f"ERROR {i}: {sra} failed to download! {check_SRA_validity.returncode}\n")
        logs.append((i, sra, "ERROR"))

timestamp = datetime.now()

#Create a new Logs directory.
os.makedirs("./Logs", exist_ok = True)

#Writing the logs to check for errors. 
with open("./Logs/SRA_download_LOGS.txt", "w") as log:
    log.write(f"# {timestamp}\n\n")
    log.write("Number\tSRA accession\tStatus")
    log.write("\n")
    for number, sra_code, status in logs:
        log.write(f"{number}\t{sra_code}\t{status}\n")

print(f"\n\nFinished downloading all available sequences!")
print(f"Please check the log file for any errors!\n")



#######################################################################
## CREATE A MANIFEST.TSV FOR QIIME2 TO USE
#######################################################################

#Create a MANIFEST file for Qiime2 for all the downloaded sequences. 
print(f"With the downloaded sequences we will create now a manifest.tsv for Qiime2 " 
       "(necessary for importing data).\n")



#######################################################################
## FOR PAIRED-END SEQUENCES
#######################################################################

#We can adapt parts of the code for other file types. 
print(f"In this case the manifest corresponds to an experiment using paired-end data")

with open("Manifest_Qiime2.tsv", "w") as m_line:  
    m_line.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for i in range(len(run_list)):
       m_line.write(f"{sample_id_list[i]}\t$PWD/../1_SRA_Download/Sequences/{run_list[i]}_1.fastq\t$PWD/../1_SRA_Download/Sequences/{run_list[i]}_2.fastq\n") 

print(f"\nFinished creating the manifest file (Manifest_Qiime2.tsv) for Qiime2!\n")
