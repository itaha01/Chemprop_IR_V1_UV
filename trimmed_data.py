import csv
import os
import re

"""
________________________________________________________________________________
This python script extracts the relevant information from all ornl_aisd_ex_X.csv files 
to create trimmed versions called ornl_aisd_ex_X_trimmed.csv
_________________________________________________________________________________
"""
# Data paths to input and output data. Modify to suit your locations
input_path = "/mnt/c/Users/Ibrahim/Desktop/AstraZeneca/ornl_aisd_ex_csv"
output_path = "Trimmed_ORNL"

# The necessary columns
smiles_pos = 2  # SMILES are the second column
intensities_start_pos = 7 #The intensities start from the 7th column
intensities_end_pos = 106 #The intensities end at the 106th column

#Extract the files and sort them numerically
ornl_files = [f for f in os.listdir(input_path) if f.startswith("ornl_aisd_ex_") and f.endswith(".csv")]
ornl_files.sort(key=lambda f: int(re.search(r'\d+', f).group()))

# Go through all ornl_files
for filename in ornl_files:
    file_number = filename.split("_")[-1].split(".")[0]

    #Create the new trimmed file
    ornl_trimmed_file = os.path.join(output_path, f"ornl_aisd_ex_{file_number}_trimmed.csv")

    # Open both of the files. Write for the new file and read for the old file
    with open(os.path.join(input_path, filename), 'r') as ornl_file, \
         open(ornl_trimmed_file, 'w', newline='') as ornl_trimmed_file:
      reader = csv.reader(ornl_file)
      writer = csv.writer(ornl_trimmed_file)

      # Extract the header to maintain it. 
      header_row = next(reader)
      trimmed_header = [header_row[smiles_pos - 1]] + header_row[intensities_start_pos-1:intensities_end_pos]

      # Write the trimmed header to the new file
      writer.writerow(trimmed_header)

      # Rowwise extract the smiles and corresponding intensities
      for row in reader:
        smiles = row[smiles_pos - 1]  # Access SMILES from the second column
        intensities = row[intensities_start_pos - 1:intensities_end_pos]
        final_row = [smiles] + intensities

        # Write the final row to the new file
        writer.writerow(final_row)

      print(f"The file '{ornl_file}'. has been trimmed to '{ornl_trimmed_file}'.")

print(f"All of the files in '{input_path}' have been trimmed to '{output_path}'.")
