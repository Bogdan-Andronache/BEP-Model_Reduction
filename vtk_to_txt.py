import vtk
import os
import re

def read_vtk_file(file_path):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    data = reader.GetOutput()
    point_data = data.GetPointData()
    if point_data.GetNumberOfArrays() == 0:
        raise ValueError(f"The VTK file {file_path} does not contain scalar data.")

    scalars = point_data.GetArray("u")
    if scalars is None:
        raise ValueError(f"The VTK file {file_path} does not contain a scalar array named 'u'.")

    scalar_values = [[scalars.GetValue(i)] for i in range(scalars.GetNumberOfTuples())]
    return scalar_values

def write_matrix_to_txt(matrix, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        for row in matrix:
            f.write(' '.join(map(str, row)) + '\n')

def main():
    input_folder = "Data/Vtk"
    output_folder = "Data/Txt"
    
    if not os.path.exists(input_folder):
        print(f"Error: Input folder '{input_folder}' does not exist.")
        return

    # Pattern to match files like sol 0.10.vtk, sol 1.00.vtk, ..., sol 10.00.vtk
    pattern = r"sol  (\d+(\.\d{1,2})?)\.vtk"
    vtk_files = sorted(
        [f for f in os.listdir(input_folder) if re.match(pattern, f)],
        key=lambda x: float(re.search(pattern, x).group(1))  # Sort by numeric value
    )

    for file_name in vtk_files:
        vtk_file = os.path.join(input_folder, file_name)
        # Remove space between 'sol' and number in the output filename
        txt_file = os.path.join(output_folder, file_name.replace(" ", "").replace(".vtk", ".txt"))
        
        try:
            scalar_values = read_vtk_file(vtk_file)
            write_matrix_to_txt(scalar_values, txt_file)
            print(f"Processed {vtk_file} -> {txt_file}")
        except ValueError as e:
            print(f"Skipping {vtk_file}: {e}")

if __name__ == "__main__":
    main()
