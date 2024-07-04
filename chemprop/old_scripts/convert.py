from rdkit import Chem
import os
import glob
import csv
import random

def process_mol_file(pre_smiles_file):
    #Convert from pdb to SMILES
    pre_smiles = Chem.MolFromPDBFile(pre_smiles_file)
    if pre_smiles is None:
        return None
    smiles = Chem.MolToSmiles(pre_smiles)

    #Extract the intensities and their coresponding wavelengths
    pre_exc_file = os.path.join(os.path.dirname(pre_smiles_file), 'EXC.DAT')
    with open(pre_exc_file, 'r') as exc:
        lines = exc.readlines()
    peaks = [0]*(401-100)

    #Read the data and skip to the first value
    for line in lines:
        if not line.strip() or line.startswith('='):
            continue

        parts = line.split()
        #Convert into nm and take range 100-400 nm
        try:
            wavelen = round(1240/float(parts[0]))
            if (100<=wavelen<=400):
                intensity = float(parts[1])
                #Adjust on wavelength value to column number
                peaks[wavelen-100] = intensity
        except ValueError:
            continue

    data_row = [smiles] + peaks
    return data_row

def create_dataset(all_mol_folder, pre_smiles_file):
    # Get the list of all 'mol_*' directories in the directory
    mol_directories = sorted(glob.glob(os.path.join(all_mol_folder, 'mol_*')))
    #Take the first x number of molecules in order
    #mol_dirs = mol_dirs[1:250000]
    #Randomize to take x randomly placed molecules
    mol_directories = random.sample(mol_directories,250000)

    # Get the 'smiles.pdb' file from each directory
    smiles_files = [os.path.join(mol_directory, 'smiles.pdb') for mol_directory in mol_directories]

    with open(pre_smiles_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['smiles'] + [str(j) for j in range(100,401)]
        writer.writerow(header)

        for smiles_file in smiles_files:
            row = process_mol_file(smiles_file)
            if row is not None:
                writer.writerow(row)

create_dataset('/home/itaha/Desktop/10.13139_OLCF_1907919/All_data','/home/itaha/chemprop-IR/chemprop/data/Medium_data_4.csv')
