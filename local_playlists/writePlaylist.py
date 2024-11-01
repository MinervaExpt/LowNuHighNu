import os

playlist = "13C"

# Define the root directory and output file
root_dir = "/pnfs/minerva/persistent/users/drut1186/LE_Tuples/P4/minerav{0}/grid/central_value/minerva/ana/v10r8p9/00/01/32/".format(playlist)
output_file = "minerva_{0}.txt".format(playlist)

# Open the output file for writing
with open(output_file, 'w') as f:
  # Walk through all directories and files in the root directory
  for dirpath, _, filenames in os.walk(root_dir):
    # Find all .root files in each directory
    for filename in filenames:
      if filename.endswith(".root"):
        # Write the full path to the output file
        filepath = os.path.join(dirpath,filename)
        filepath = filepath.lstrip("/pnfs")
        filepath_xrootd = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/" + filepath
        print("Adding this file: {0}".format(filepath_xrootd))
        f.write(filepath_xrootd + "\n")

print(f"File paths have been saved to {output_file}")

