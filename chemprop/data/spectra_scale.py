"""Min-max normalization to scale the spectra"""
import numpy as np
import pandas as pd
import sys

"""
_________________________________________________________________________________
This python script applies min-max normalization to the predicted .csv spectra which  
This produces the scaled_spectra.csv file containing all of the min_max normalized predicted spectra.
_________________________________________________________________________________
"""

# Min-max normalization for one SMILES
def min_max_scale_spectra(intensities_list):
    for i in range(len(intensities_list)):
        intensities = intensities_list.iloc[i, 1:]  # Skip the first column containing SMILES
        min_max_intensities = (intensities - intensities.min()) / (intensities.max() - intensities.min())
        intensities_list.iloc[i, 1:] = min_max_intensities  # Change the intensities to min-max and keep the SMILES identical
    return intensities_list

if __name__ == "__main__":
    raw_spectra = pd.read_csv(sys.argv[1],index_col=0)
    scaled_spectra = min_max_scale_spectra(raw_spectra)
    scaled_spectra.to_csv('/home/itaha/chemprop-IR/chemprop/data/scaled_spectra.csv')
