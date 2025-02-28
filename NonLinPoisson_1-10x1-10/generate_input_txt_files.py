import numpy as np

# Load the parameters from the .npy file
parameters = np.load("parameters.npy")

# Template for the text file
template = """&comppar
  outputfile = 'sol_{:.1f}_{:.1f}.vtk'    ! name of the output file
  delta = {:.1f}               ! delta coefficient
  gamma = {:.1f}               ! gamma coefficient
  nx = 20                   ! number of elements in x-direction
  ny = 20                   ! number of elements in y-direction
  eps = 1e-9                ! accuracy in iteration
  nPicard = 10,             ! number of Picard iterations before Newton-Raphson
  itermax = 100             ! maximum iterations
/
"""

# Loop through each row in parameters and create a file
for delta, gamma in parameters:
    filename = f"parameters_{delta:.1f}_{gamma:.1f}.txt"
    content = template.format(delta, gamma, delta, gamma)
    
    with open(filename, "w") as file:
        file.write(content)
    
    print(f"Generated: {filename}")
