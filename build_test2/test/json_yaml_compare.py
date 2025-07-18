import os
import filecmp
import sys

def compare_files(file_name, folder_path1, folder_path2):
  file1 = os.path.join(folder_path1, file_name)
  file2 = os.path.join(folder_path2, file_name)

  if os.path.isfile(file1) and os.path.isfile(file2):
    if filecmp.cmp(file1, file2):
      print("The files are equal.")
    else:
      print("The files are not equal.")
      return 1  # Return a failure code
  else:
    print("One or both files do not exist.")
    return 1  # Return a failure code

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python script.py <folder_path1> <folder_path2>")
  else:
    folder_path1 = sys.argv[1]
    folder_path2 = sys.argv[2]
    compare_files("photolysis_rate_constants.nc", folder_path1, folder_path2)
    compare_files("dose_rates.nc", folder_path1, folder_path2)
