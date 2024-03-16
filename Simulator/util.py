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

    # Find peaks
    peaks, _ = find_peaks(y, distance=20)  # The distance parameter may need to be adjusted

    # Sort the peaks based on the peak heights in descending order
    sorted_peak_indices = np.argsort(y[peaks])[::-1]

    # Get the indices of the two largest peaks
    if len(sorted_peak_indices) >= 2:
        # If there are two or more peaks, take the first two
        top_peaks_indices = sorted_peak_indices[:2]
    elif len(sorted_peak_indices) == 1:
        # If there's only one peak, take it
        top_peaks_indices = sorted_peak_indices
    else:
        # If there are no peaks, return empty arrays
        return np.array([]), np.array([])

    # Find the positions of the two largest peaks
    peak_positions = x[peaks][top_peaks_indices]
    peak_values = y[peaks][top_peaks_indices]

    return list(peak_positions), list(peak_values)



if __name__ == "__main__":
    A = np.array([[0, 1], [1, 1j]], dtype=complex)
    B = np.array([[0, 1.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)

    assert all_close(A, B)
    print(norm_2_distance(A, B))

    C = np.array([[0, 2.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)
    assert not all_close(A, C)
    print(norm_2_distance(A, C))
