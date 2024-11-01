import os

playlist = "13C"

# Define the root directory and output file
root_dir_mc = "/pnfs/minerva/persistent/users/drut1186/LE_Tuples/P4/minerav{0}/grid/central_value/minerva/ana/v10r8p9/00/01/32/".format(playlist)
output_file_mc = "minerva{0}_MC.txt".format(playlist)

root_dir_data = "/pnfs/minerva/persistent/users/drut1186/LE_Tuples/P4/data/minerav{0}/grid/minerva/ana/numibeam/v10r8p9/00/00/".format(playlist)
output_file_data = "minerva{0}_Data.txt".format(playlist)

for root_dir,output_file in [(root_dir_mc,output_file_mc),(root_dir_data,output_file_data)]:
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

