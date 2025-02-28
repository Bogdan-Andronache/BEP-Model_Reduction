import numpy as np

array = np.array([[i, j] for i in range(1, 11) for j in range(1, 11)])

# Save as a text file
np.savetxt("parameters.txt", array, fmt="%d", delimiter=" ")

# Save as a .npy file
np.save("parameters.npy", array)

print("Files saved: parameters.txt and parameters.npy")
