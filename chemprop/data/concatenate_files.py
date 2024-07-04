import os
import pandas as pd

"""
_________________________________________________________________________________
This python script concatenates multiple ornl files of the same header, specifically
ornl_aisd_ex_X_trimmed_smoothed_min_max.csv files
_________________________________________________________________________________
"""

# Define directory containing the smoothed ORNL files
data_dir = '/home/itaha/chemprop-IR/Trimmed_Smoothed_Min_Max_ORNL_6nm/'
n_files = 999 
output_file = "/home/itaha/chemprop-IR/chemprop/data/6nm_dataset/10.5M_smoothed_6nm_min_max.csv"
combined_df = None

for i in range(1, n_files + 1):
  # Create file path for current file
  file_path = os.path.join(data_dir, f"ornl_aisd_ex_{i}_trimmed_smoothed_min_max.csv")

  # Check if file exists
  if os.path.exists(file_path):
    current_df = pd.read_csv(file_path)

    # If it's the first file, initialize combined_df
    if combined_df is None:
      combined_df = current_df.copy()
    # Otherwise, concatenate current_df to combined_df
    else:
      combined_df = pd.concat([combined_df, current_df])
  else:
    print(f"The file '{file_path}' could not be found.")

# Check if any files were concatenated
if combined_df is not None:
  combined_df.to_csv(output_file, index=False)
  print(f"Successfully concatenated {n_files} files and saved to '{output_file}'.")
else:
  print(f"No files were found in '{data_dir}'.")
