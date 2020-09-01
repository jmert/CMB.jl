import h5sparse
import scipy.sparse as sparse
import numpy as np
import pathlib as path
import os

Rdense = np.array([
    [0.1, 0.2, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]], dtype = np.float64)

base = path.Path(__file__).resolve().parent.parent

# Generate a CSC matrix, as h5sparse generates by default. We're explicit about the
# index data types to more clearly differentiate it from the following case.
csc_mixed = base.joinpath("obsmat_sparse_pycsc_i64i32f64.h5")
with h5sparse.File(csc_mixed, "w") as f:
    f.create_dataset("R", data = sparse.csc_matrix(Rdense), dtype = np.float64,
            indptr_dtype = np.int64, indices_dtype = np.int32)

# Generate another CSC matrix, but this time force the index types to match.
csc_uniform = base.joinpath("obsmat_sparse_pycsc_i32i32f64.h5")
with h5sparse.File(csc_uniform, "w") as f:
    f.create_dataset("R", data = sparse.csc_matrix(Rdense), dtype = np.float64,
            indptr_dtype = np.int32, indices_dtype = np.int32)

# Give the default case a simple name.
csc_simple = base.joinpath("obsmat_sparse_pycsc.h5")
if not csc_simple.exists():
    os.symlink(csc_mixed, csc_simple)

# Generate a CSR matrix as well
with h5sparse.File(base.joinpath("obsmat_sparse_pycsr.h5"), "w") as f:
    f.create_dataset("R", data = sparse.csr_matrix(Rdense))
