import numpy as np


def all_close(matrix1, matrix2):
    return np.allclose(np.real(matrix1), np.real(matrix2), rtol=1e-8, atol=1e-5) and np.allclose(np.imag(matrix1), np.imag(matrix2), rtol=1e-8, atol=1e-5)


def hilbert_schmidt_distance(V: np.ndarray, U: np.ndarray):
    L = V.shape[0]
    Vdag = V.transpose()
    Vdag = Vdag.conjugate()
    return np.sqrt(1 - np.abs(np.abs(np.trace(np.matmul(Vdag, U))) ** 2) / (L ** 2))


def norm_2_distance(matrix1: np.ndarray, matrix2: np.ndarray):
    return np.abs(np.sum((matrix1 - matrix2) ** 2))


if __name__ == "__main__":
    A = np.array([[0, 1], [1, 1j]], dtype=complex)
    B = np.array([[0, 1.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)

    assert all_close(A, B)
    print(norm_2_distance(A, B))

    C = np.array([[0, 2.00000000001 + 0.000000001 * 1j], [1, 1j]], dtype=complex)
    assert not all_close(A, C)
    print(norm_2_distance(A, C))
