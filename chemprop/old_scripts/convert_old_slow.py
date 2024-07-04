import os
import glob
import csv
import random
from scipy.ndimage import gaussian_filter1d
from rdkit import Chem

def process_mol_file(pre_smiles_file):
    # Convert from pdb to SMILES
    pre_smiles = Chem.MolFromPDBFile(pre_smiles_file)
    if pre_smiles is None:
        return None
    smiles = Chem.MolToSmiles(pre_smiles)

    # Extract the intensities and their coresponding wavelengths
    pre_exc_file = os.path.join(os.path.dirname(pre_smiles_file), 'EXC.DAT')
    with open(pre_exc_file, 'r') as exc:
        lines = exc.readlines()
    peaks = [0] * (401 - 100)

    # Read the data and skip to the first value
    for line in lines:
        if not line.strip() or line.startswith('='):
            continue

        parts = line.split()
        # Convert into nm and take range 100-400 nm
        try:
            wavelen = round(1240 / float(parts[0]))
            if (100 <= wavelen <= 400):
                intensity = float(parts[1])
                # Adjust on wavelength value to column number
                peaks[wavelen - 100] = intensity
        except ValueError:
            continue

    # Apply Gaussian filter for smoothing (sigma to be adjusted)
    sigma = 5  # Adjust this parameter for desired smoothing strength
    smoothed_peaks = gaussian_filter1d(peaks, sigma=sigma)

    data_row = [smiles] + smoothed_peaks.tolist()
    return data_row

def create_dataset(all_mol_folder, pre_smiles_file):
    # Get the list of all 'mol_*' directories in the directory
    mol_directories = sorted(glob.glob(os.path.join(all_mol_folder, 'mol_*')))
    smiles_files = []

    for mol_directory in mol_directories:
        try:
            smiles_file = os.path.join(mol_directory, 'smiles.pdb')
            # Check if the smiles.pdb file can be read successfully
            Chem.MolFromPDBFile(smiles_file)  # Raise OSError if file is bad
            smiles_files.append(smiles_file)
        except OSError as e:
            print(f"Skipping folder '{mol_directory}' due to error: {e}")

    # Get the 'smiles.pdb' file from each directory (using valid files)
    # smiles_files = [os.path.join(mol_directory, 'smiles.pdb') for mol_directory in mol_directories]

    with open(pre_smiles_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['smiles'] + [str(j) for j in range(100, 401)]
        writer.writerow(header)

        for smiles_file in smiles_files:
            row = process_mol_file(smiles_file)
            if row is not None:
                writer.writerow(row)

# Example usage:
create_dataset('/home/itaha/MVEX03/dftb_gdb9_electronic_excitation_spectrum',
               '/home/itaha/chemprop-IR/chemprop/data/Small_filtered.csv')