import os
target_dir = "/home/itaha/Desktop/10.13139_OLCF_1907919/All_data"

# Loop through folders starting with "mol_"
for folder in os.listdir(target_dir):
  if folder.startswith("mol_"):
    folder_path = os.path.join(target_dir, folder)
    
    # Keep track of desired files
    keep_files = ["smiles.pdb", "EXC.DAT"]

    # Loop through files within the folder
    for filename in os.listdir(folder_path):
      # Delete if not in the keep list
      if filename not in keep_files:
        file_path = os.path.join(folder_path, filename)
        os.remove(file_path)
        #print(f"Deleted: {file_path}")  # Optional: Print deleted filenames
