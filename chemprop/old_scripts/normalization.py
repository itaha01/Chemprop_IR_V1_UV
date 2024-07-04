import csv

def normalize_csv(input_csv_file, output_csv_file):
    """
    Normalizes the values in each row of a CSV file to sum to 1 and saves to a new file.

    Args:
        input_csv_file (str): The path to the input CSV file.
        output_csv_file (str): The path to the output CSV file.
    """

    with open(input_csv_file, 'r', newline='') as infile, \
         open(output_csv_file, 'w', newline='') as outfile:

        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        header = next(reader)  # Read the header row
        writer.writerow(header)  # Write the header to the output file

        for row in reader:
            smiles = row[0]
            intensities = [float(x) for x in row[1:]]
            total_intensity = sum(intensities)
            normalized_intensities = [x / total_intensity for x in intensities]
            writer.writerow([smiles] + normalized_intensities)

# Example usage (replace with your actual filepaths)
input_csv_file = '/home/itaha/chemprop-IR/chemprop/data/scaled_spectra.csv'
output_csv_file = '/home/itaha/chemprop-IR/chemprop/data/scaled_norm_spectra.csv'
normalize_csv(input_csv_file, output_csv_file)
