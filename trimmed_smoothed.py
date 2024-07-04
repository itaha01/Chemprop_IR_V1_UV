import os
import re
import shutil
import numpy as np
import scipy
import pandas as pd

from scipy.signal import find_peaks
from scipy import signal
"""
_________________________________________________________________________________
This python script converts from eV to nm and applies Gaussian smoothing to
all ornl_aisd_ex_X_trimmed.csv files while using a 6 nm discretization and only 
looking at the range between 150 and 450 nm. This produces the ornl_aisd_ex_X_trimmed_smoothed.csv files.
_________________________________________________________________________________
"""

# Parameters
w_nm = 10.0
w_max = 451
w_min = 150
w_step = 6

# Convert from eV to nm with equation lambda = hc/E
def eV_to_nm(eV):
  planck = 4.1357 * 1e-15  # eV s (This is h)
  c_speed = 299792458  # m / s (This is c)
  m_to_nm = 1e+9 # Convert from m to nm
  nm = 1 / eV * planck * c_speed * m_to_nm
  return nm

# Gaussian smoothing function
def gauss(a, m, x, w):
  return a * np.exp(-(np.log(2) * ((m - x) / w) ** 2))

# Extract the SMILES, intensities and eV from ornl files (.csv format only)
def read_csv(csv_file):
  df = pd.read_csv(csv_file)
  smiles = df['smiles'].values
  intensities_list = df.values[:, 51:101].tolist()  # Intensities are in columns 52 to 101
  eV_list = df.values[:, 1:51].tolist()  # eV values are in columns 2 to 51
  return smiles, eV_list, intensities_list

# Apply Gaussian smoothing to the intensities by using a Gaussian filter
def spectra_smoothing(eV_list, intensities_list):
  nm_list = [eV_to_nm(eV) for eV in eV_list]
  nm_list.sort()
  gauss_intensities = [] #This is a list of all the gaussian smoothed elements
  nm_range = np.arange(w_min, w_max, w_step)
  for i, wn in enumerate(nm_list):
    gauss_intensity = gauss(intensities_list[i], nm_range, wn, w_nm)
    gauss_intensities.append(gauss_intensity)
  gauss_intensities_final = np.sum(gauss_intensities, axis=0) #Sum the gaussian smoothed elements column wise
  return gauss_intensities_final

# Process the file and write all smoothed data to a single output file
def process_file(input_file, output_file):
  if os.path.exists(input_file):
    with open(output_file, "w") as output_file:  # Open in write mode to create a new file
      smiles, eV_list, intensities_list = read_csv(input_file)
      header_line = "smiles," + ",".join([str(w_min + i * w_step) for i in range(len(spectra_smoothing(eV_list[0], intensities_list[0])))]) + "\n"
      output_file.write(header_line)
      for smile, ev_data, intensity_data in zip(smiles, eV_list, intensities_list):
        smooth_spectra = spectra_smoothing(ev_data, intensity_data)
        output_line = smile + "," + ",".join([str(intensity) for intensity in smooth_spectra]) + "\n"
        output_file.write(output_line)
    print(f"The file '{input_file}' has been processed as '{output_file}'.")
  else:
    print(f"The file '{input_file}' could not be found.")

# Main execution
if __name__ == "__main__":
  input_path = '/home/itaha/chemprop-IR/Trimmed_ORNL/'
  output_path = '/home/itaha/chemprop-IR/Trimmed_Smoothed_ORNL_6nm/'
  for i in range(1, 1001):  # For X=1,1000
    input_file = os.path.join(input_path, f"ornl_aisd_ex_{i}_trimmed.csv")
    output_file = os.path.join(output_path, f"ornl_aisd_ex_{i}_trimmed_smoothed.csv")
    process_file(input_file, output_file)
  print(f"Finished processing all of the ornl_files.")
