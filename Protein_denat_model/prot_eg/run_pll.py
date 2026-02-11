"""
MPI-parallel compilation and execution of parameterized C++ simulations.

Each MPI rank receives a subset of run indices, compiles a unique executable,
and executes it with the corresponding index as an argument.
"""

import os
import sys
import numpy as np
from mpi4py import MPI

# ---------------------------------------------------------------------
# MPI initialization
# ---------------------------------------------------------------------
comm = MPI.COMM_WORLD
size = comm.Get_size()   # Total number of MPI ranks
rank = comm.Get_rank()   # Rank ID of this process

# ---------------------------------------------------------------------
# Input arguments
# ---------------------------------------------------------------------
# sys.argv[1]: total number of simulation runs
# sys.argv[2]: C++ source file to compile
n_runs = int(sys.argv[1])
source_file = sys.argv[2]

# ---------------------------------------------------------------------
# Create list of run indices
# ---------------------------------------------------------------------
range_runs = np.arange(n_runs)

# Number of runs per rank (integer division)
num_per_rank = len(range_runs) // size

# ---------------------------------------------------------------------
# Scatter equal-sized chunks across ranks
# ---------------------------------------------------------------------
data_s1 = None
if rank == 0:
    # Only rank 0 prepares the data to scatter
    data_s1 = range_runs[:num_per_rank * size]

recvbuf = np.empty(num_per_rank, dtype=np.int64)
comm.Scatter(data_s1, recvbuf, root=0)

# ---------------------------------------------------------------------
# Handle remainder runs (if n_runs not divisible by size)
# These are distributed manually from rank 0
# ---------------------------------------------------------------------
data_s2 = None
remainder = range_runs[num_per_rank * size:]

if len(remainder) != 0:
    for i in range(len(remainder)):
        if rank == 0:
            comm.send(remainder[i], dest=i)
        elif rank == i:
            data_s2 = comm.recv(source=0)

# ---------------------------------------------------------------------
# Combine scattered data and remainder (if any)
# ---------------------------------------------------------------------
if data_s2 is not None:
    data = np.insert(recvbuf, 0, data_s2)
else:
    data = recvbuf

print(f"Rank {rank} received runs: {data}")

# ---------------------------------------------------------------------
# Compile and run simulations assigned to this rank
# ---------------------------------------------------------------------
for i in data:
    exe_name = f"a_{int(i)}.exe"

    # Compile C++ source
    os.system(
        f"g++ -o {exe_name} {source_file} -lm -lgsl -lgslcblas"
    )

    # Execute compiled binary
    os.system(
        f"./{exe_name} {int(i)}"
    )

