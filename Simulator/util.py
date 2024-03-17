import numpy as np
from scipy.signal import find_peaks

def all_close(matrix1, matrix2):
    return np.allclose(np.real(matrix1), np.real(matrix2), rtol=1e-8, atol=1e-5) and np.allclose(np.imag(matrix1), np.imag(matrix2), rtol=1e-8, atol=1e-5)


def hilbert_schmidt_distance(V: np.ndarray, U: np.ndarray):
    Vdag = V.transpose()
    Vdag = Vdag.conjugate()
    return np.sqrt(np.abs(1 - np.abs(np.abs(np.trace(np.matmul(Vdag, U))) ** 2) / (4 ** 2)))


def norm_2_distance(matrix1: np.ndarray, matrix2: np.ndarray):
    return np.abs(np.sum((matrix1 - matrix2) ** 2))


def find_two_largest_peaks(x, y):
    # Convert to numpy arrays for processing
    x = np.array(x)
    y = np.array(y)

    # Find all peaks, both positive and negative, by finding peaks on the original and the negated signal
    positive_peaks, _ = find_peaks(y, distance=20)
    negative_peaks, _ = find_peaks(-y, distance=20)

    # Combine the found positive and negative peaks
    all_peaks = np.concatenate((positive_peaks, negative_peaks))

    # Check if we have found at least two peaks
    if len(all_peaks) >= 2:
        # Get the absolute values of the peaks to sort them
        peak_values = y[all_peaks]
        abs_peak_values = np.abs(peak_values)

        # Sort the peaks by the absolute values in descending order and take the indices of the two largest
        largest_peaks_indices = all_peaks[np.argsort(abs_peak_values)[-2:]]

    elif len(all_peaks) == 1:
        # If there's only one peak, take it
        largest_peaks_indices = all_peaks
    else:
        # If there are no peaks, return empty arrays
        return np.array([]), np.array([])

    # Find the positions and values of the two largest absolute peaks
    peak_positions = x[largest_peaks_indices]
    peak_values = y[largest_peaks_indices]

    return list(peak_positions), list(peak_values)


def store_string_to_file(output_string, file_path):
    """
    Stores a given string to a specified file.

    :param output_string: The string to be stored in the file
    :param file_path: The path of the file where the string will be stored
    """
    with open(file_path, 'w') as file:
        file.write(output_string)



if __name__ == "__main__":
    A = np.array([[0, 1], [1, 1j]], dtype=complex)
    B = np.array([[0, 1.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)

    assert all_close(A, B)
    print(norm_2_distance(A, B))

    C = np.array([[0, 2.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)
    assert not all_close(A, C)
    print(norm_2_distance(A, C))
