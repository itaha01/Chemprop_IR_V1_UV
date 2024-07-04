from rdkit import Chem
import os
import glob
import csv
import random

def process_mol_file(pre_smiles_file):
    # Convert from pdb to SMILES
    pre_smiles = Chem.MolFromPDBFile(pre_smiles_file)
    if pre_smiles is None:
        return None
    smiles = Chem.MolToSmiles(pre_smiles)

    # Extract the intensities and their corresponding wavelengths
    pre_exc_file = os.path.join(os.path.dirname(pre_smiles_file), 'EXC_smoothed.csv')
    if not os.path.exists(pre_exc_file):
        return None

    peaks = {}  # Use a dictionary to store intensities with wavelengths as keys

    with open(pre_exc_file, 'r') as exc:
        reader = csv.reader(exc)
        for line in reader:
            if len(line) == 2:
                wavelength = int(line[0])
                intensity = float(line[1])
                if 150 <= wavelength <= 448:  # Consider wavelengths within the desired range
                    peaks[wavelength] = intensity

    # Create the data row
    data_row = [smiles] + [peaks.get(wavelength, 0) for wavelength in range(150, 449,2)]
    #print(peaks)
    return data_row

def create_dataset(all_mol_folder, pre_smiles_file):
    # Get the list of all 'mol_*' directories in the directory
    mol_directories = sorted(glob.glob(os.path.join(all_mol_folder, 'mol_*')))
    mol_directories = mol_directories[1:500000]
    #Randomize to take x randomly placed molecules
    #mol_directories = random.sample(mol_directories,500000)

    # Get the 'smiles.pdb' file from each directory
    smiles_files = [os.path.join(mol_directory, 'smiles.pdb') for mol_directory in mol_directories]

    with open(pre_smiles_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['smiles'] + [str(j) for j in range(150,449,2)]
        writer.writerow(header)

        for smiles_file in smiles_files:
            row = process_mol_file(smiles_file)
            if row is not None:
                writer.writerow(row)

create_dataset('/home/itaha/Desktop/10.13139_OLCF_1907919/Subdata_500','/home/itaha/chemprop-IR/chemprop/data/500K_smoothed.csv')
