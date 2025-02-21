import subprocess
import os
import numpy as np
from scipy.sparse import csr_matrix

def create_input_file(filename, alpha_value):
   # Create an input file with specified alpha value in namelist format
   with open(filename, 'w') as f:
       f.write("&comppar\n")
       f.write(f"alpha = {alpha_value:.4f}\n")
       f.write("/\n")

def run_fortran_simulation(input_file, output_file):
   # Execute the Fortran simulation with the specified input file and output file
   with open(output_file, "w") as outfile:
       subprocess.run(["C:\\cygwin64\\bin\\bash.exe", "-c", "C:\\cygwin64\\bin\\poisson1.exe"], input=open(input_file).read(), stdout=outfile, text=True)

def append_output_files(input_files, output_combined):
   # Open the combined output file for writing
   with open(output_combined, 'w') as f_combined:
       # Iterate over each input file
       for input_file in input_files:
           # Read the contents of the current output file and append to combined file
           with open(input_file.replace("input_", "output_"), 'r') as f_output:
               f_combined.write(f_output.read())

def write_alpha_values(alpha_values, filename):
   # Write the alpha values to a specified file
   with open(filename, 'w') as f:
       for alpha in alpha_values:
           f.write(f"{alpha:.4f}\n")


alpha_values = np.linspace(1.0, 10.1, 50)

# Write the alpha values to a new file
write_alpha_values(alpha_values, "parameters.txt")
print("Alpha values have been saved to parameters.txt.")

if __name__ == "__main__":
   # Create input files with different alpha values
   input_files = []
   for alpha in alpha_values:
       filename = f"input_{alpha:.4f}.txt"
       create_input_file(filename, alpha)
       input_files.append(filename)

   # Run simulations for each input file
   for input_file in input_files:
       output_file = input_file.replace("input_", "output_")
       run_fortran_simulation(input_file, output_file)

   # Append output from all simulations into a single combined output file
   output_combined = "snapshots.txt"
   append_output_files(input_files, output_combined)
   print(f"snapshots created.")