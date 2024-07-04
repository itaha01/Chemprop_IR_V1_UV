import pandas as pd
import os

"""
_________________________________________________________________________________
This python script applies normalization to all ornl_aisd_ex_X_trimmed_smoothed.csv  
This produces the ornl_aisd_ex_X_trimmed_smoothed_norm.csv files.
_________________________________________________________________________________
"""

# Normalization for one SMILES
def norm_scale_spectra(intensities_list):
    for i in range(len(intensities_list)):
        intensities = intensities_list.iloc[i, 1:]  # Skip the first column containing SMILES
        norm_intensities = (intensities) / (sum(intensities))
        intensities_list.iloc[i, 1:] = norm_intensities  # Change the intensities to min-max and keep the SMILES identical
    return intensities_list

# Function to load, apply min-max normalization and save.
def process_file(input_file, output_file):
    intensities = pd.read_csv(input_file, index_col=0)
    intensities_norm_normalized = norm_scale_spectra(intensities)
    intensities_norm_normalized.to_csv(output_file)

# Main execution
if __name__ == "__main__":
    input_dir = '/home/itaha/chemprop-IR/Trimmed_Smoothed_ORNL_6nm/'
    output_dir = '/home/itaha/chemprop-IR/Trimmed_Smoothed_ORNL_Norm_6nm/'
    for i in range(1, 1001):  # For X=1,1000
        input_file = os.path.join(input_dir, f"ornl_aisd_ex_{i}_trimmed_smoothed.csv")
        output_file = os.path.join(output_dir, f"ornl_aisd_ex_{i}_trimmed_smoothed_norm.csv")
        if os.path.exists(input_file):
            process_file(input_file, output_file)
            print(f"The file '{input_file}' has been processed as '{output_file}'.")
        else:
            print(f"The file '{input_file}' could not be found.")
    print(f"Finished processing all of the ornl_files.")
