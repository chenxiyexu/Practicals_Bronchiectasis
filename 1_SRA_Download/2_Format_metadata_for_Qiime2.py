#These scripts were written by Chen Xi Ye Xu as part of her Master´s Thesis 
# for UAM - Master of Bioinformatics and Computational Biology. 

#Load the neccesary libraries.
import pandas as pd


def format_qiime2_metadata():
    """ 
    This function takes a metadata file (CSV or TSV) and formats it for QIIME 2 by 
    renaming the specified sample ID column to 'sample-id', reordering the columns
    and saving it as a tab-separated file.
    
    Input:
        - input_file (str): path to the input metadata file (CSV or TSV).

    Output:
        - metadata_Qiime_adapted.tsv (tab-separated file/ .tsv): formatted for QIIME 2.
    
    """


    #1)Ask for name of the file. 
    input_file = input("Enter your metadata file name (with file extension like .txt or .csv): ").strip()

    
    #2)Load the file (it automatically detects csv or tsv)
    try:
        if input_file.endswith(".csv"):
            df = pd.read_csv(input_file)
        else:
            df = pd.read_csv(input_file, sep = None, engine = "python")
    except Exception as e:
        print(f"ERROR reading file: {e}")
        return


    #3)Display available columns from the table
    print("\nAvailable columns in your metadata file[cite: 1]:")
    for col in df.columns:
        print(f"   - {col}")

        
    #4)Ask which column matches the sample IDs used in the manifest
    target_col = input("\nEnter the exact name of the column you want to use as your sample-id (for example: id)[cite: 1]: ").strip()
    
    if target_col not in df.columns:
        print(f"\nERROR: Column {target_col} not found in the file.")
        return


    #5)Rename the chosen column to "sample-id" (QIIME 2 requires it)
    df = df.rename(columns = {target_col: "sample-id"})

    
    #6)Reorder columns so 'sample-id' is strictly the very first column
    cols = ["sample-id"] + [col for col in df.columns if col != "sample-id"]
    df = df[cols]

    
    #7)Save as a tab-separated file (.tsv) as required by QIIME 2
    output_file = "metadata_Qiime_adapted.tsv"
    df.to_csv(output_file, sep = "\t", index = False)

    
    print(f"\nSuccess! Formatted metadata saved as: {output_file}.")
    print("You can now enter metadata_Qiime_adapted.tsv when running Qiime2 bash scripts.")


if __name__ == "__main__":
    format_qiime2_metadata()
