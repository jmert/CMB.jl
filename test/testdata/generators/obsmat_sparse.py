import h5sparse
import scipy.sparse as sparse
import numpy as np
import pathlib as path

Rdense = np.array([
    [0.1, 0.2, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]], dtype = np.float64)

base = path.Path(__file__).resolve().parent.parent

with h5sparse.File(base.joinpath("obsmat_sparse_pycsc.h5"), "w") as f:
    f.create_dataset("R", data = sparse.csc_matrix(Rdense))

with h5sparse.File(base.joinpath("obsmat_sparse_pycsr.h5"), "w") as f:
    f.create_dataset("R", data = sparse.csr_matrix(Rdense))
