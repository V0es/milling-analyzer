import numpy as np
from numba.core.extending import register_jitable
from numpy.linalg import inv


@register_jitable
def npexpm(matrix: np.ndarray):
    eigvals, eigvecs = np.linalg.eig(matrix)
    diag_exp = np.diag(np.exp(eigvals))
    return eigvecs @ diag_exp @ inv(eigvecs).astype(np.complex128)