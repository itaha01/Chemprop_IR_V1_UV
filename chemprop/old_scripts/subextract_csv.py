#Extracts the first num_rows of the original_file to create a subset called new_file.
def subextract_csv(original_file, new_file, num_rows):
  with open(original_file, 'r') as original, open(new_file, 'w') as new:
    row_count = 0
    for line in original:
      if row_count < num_rows:  # Include up to (but not exceeding) the specified number of rows
        new.write(line)
      row_count += 1

original_file = "/home/itaha/chemprop-IR/chemprop/data/500K_smoothed.csv"
new_file = "/home/itaha/chemprop-IR/chemprop/data/50K_smoothed_2nm.csv"
num_rows = 50000 + 1  # + 1 to include the header

subextract_csv(original_file, new_file, num_rows)

print(f"Created new file: {new_file}")
