VieClus v1.2
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-11/14-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.10+-064F8C.svg)](https://cmake.org/)
[![Build](https://github.com/KaHIP/VieClus/actions/workflows/build.yml/badge.svg)](https://github.com/KaHIP/VieClus/actions/workflows/build.yml)
[![GitHub Release](https://img.shields.io/github/v/release/KaHIP/VieClus)](https://github.com/KaHIP/VieClus/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/vieclus)](https://pypi.org/project/vieclus/)
[![Homebrew](https://img.shields.io/badge/Homebrew-available-orange)](https://github.com/KaHIP/homebrew-kahip)
[![Linux](https://img.shields.io/badge/Linux-supported-success.svg)](https://github.com/KaHIP/VieClus)
[![macOS](https://img.shields.io/badge/macOS-supported-success.svg)](https://github.com/KaHIP/VieClus)
[![GitHub Stars](https://img.shields.io/github/stars/KaHIP/VieClus)](https://github.com/KaHIP/VieClus/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/KaHIP/VieClus)](https://github.com/KaHIP/VieClus/issues)
[![Last Commit](https://img.shields.io/github/last-commit/KaHIP/VieClus)](https://github.com/KaHIP/VieClus/commits)
[![arXiv](https://img.shields.io/badge/arXiv-1802.07034-b31b1b.svg)](https://arxiv.org/abs/1802.07034)
[![Heidelberg University](https://img.shields.io/badge/Heidelberg-University-c1002a)](https://www.uni-heidelberg.de)
=====

<p align="center">
  <img src="./logo/vieclus-logo.svg" alt="VieClus Logo" width="900"/>
</p>

The graph clustering framework VieClus -- Vienna Graph Clustering. Part of the [KaHIP](https://github.com/KaHIP) organization.

Graph clustering is the problem of detecting tightly connected regions of a
graph. Depending on the task, knowledge about the structure of the graph can
reveal information such as voter behavior, the formation of new trends, existing
terrorist groups and recruitment or a natural partitioning of
data records onto pages. Further application areas
include the study of protein interaction, gene
expression networks, fraud
detection, program optimization and the spread of
epidemics---possible applications are plentiful, as
almost all systems containing interacting or coexisting entities can be modeled
as a graph. 



This is the release of our memetic algorithm, VieClus (Vienna Graph Clustering), to tackle the graph clustering problem. 
A key component of our contribution are natural recombine operators that employ ensemble clusterings as well as multi-level techniques. 
In our experimental evaluation, we show that **our algorithm successfully improves or reproduces all entries of the 10th DIMACS implementation challenge** under consideration in a small amount of time. In fact, for most of the small instances, we can improve the old benchmark result in less than a minute.
Moreover, while the previous best result for different instances has been computed by a variety of solvers, our algorithm can now be used as a single tool to compute the result. **In short our solver is the currently best modularity based clustering algorithm available.**

<p align="center">
<img src="./img/example_clustering.png"
  alt="example clustering"
  width="538" height="468">
</p>

Installation Notes
=====

### Install via Homebrew

```bash
brew install KaHIP/kahip/vieclus
```

### C++ Command Line Tool

VieClus can be compiled with or without MPI support.

#### With MPI (recommended for best solution quality)

MPI enables the parallel evolutionary algorithm which typically yields better solutions.

Prerequisites:
- OpenMPI (http://www.open-mpi.org/) -- note: due to removed progress threads in OpenMPI > 1.8, please use an OpenMPI version < 1.8 or Intel MPI to obtain a scalable parallel algorithm.

```bash
./compile_withcmake.sh
mpirun -n 2 ./deploy/vieclus examples/astro-ph.graph --time_limit=60
```

#### Without MPI (NOMPI)

If you do not have MPI installed or only need single-process execution, you can compile without MPI support. The algorithm will run on a single process using a pseudo-MPI layer.

```bash
./compile_withcmake.sh NOMPI
./deploy/vieclus examples/astro-ph.graph --time_limit=60
```

No additional dependencies are required for the NOMPI build.

For a description of the graph format please have a look into the manual.

Python Interface
=====

You can install the Python interface via pip:
``pip install vieclus``

Or build from source:
``pip install .``

### Example: Using the vieclus_graph class

```python
import vieclus

# Build a graph using the vieclus_graph helper class
g = vieclus.vieclus_graph()
g.set_num_nodes(6)

# Add edges (undirected, with weights)
g.add_undirected_edge(0, 1, 5)
g.add_undirected_edge(1, 2, 5)
g.add_undirected_edge(0, 2, 5)
g.add_undirected_edge(3, 4, 5)
g.add_undirected_edge(4, 5, 5)
g.add_undirected_edge(3, 5, 5)
g.add_undirected_edge(2, 3, 1)  # weak bridge between two communities

# Convert to CSR format and cluster
vwgt, xadj, adjcwgt, adjncy = g.get_csr_arrays()
modularity, clustering = vieclus.cluster(vwgt, xadj, adjcwgt, adjncy,
                                         time_limit=1.0)

print(f"Modularity: {modularity}")
print(f"Clustering: {clustering}")
```

### Example: Using raw CSR arrays

```python
import vieclus

# Graph in METIS CSR format (same as KaHIP)
xadj   = [0, 2, 5, 7, 9, 12]
adjncy = [1, 4, 0, 2, 4, 1, 3, 2, 4, 0, 1, 3]
vwgt   = [1, 1, 1, 1, 1]
adjcwgt = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

modularity, clustering = vieclus.cluster(vwgt, xadj, adjcwgt, adjncy,
                                         suppress_output=True,
                                         seed=0,
                                         time_limit=2.0)

print(f"Modularity: {modularity}")
print(f"Clustering: {clustering}")
```

### Parameters

The `vieclus.cluster` function takes the following arguments:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vwgt` | list | *required* | Node weights (length n) |
| `xadj` | list | *required* | CSR index array (length n+1) |
| `adjcwgt` | list | *required* | Edge weights (length m) |
| `adjncy` | list | *required* | CSR adjacency array (length m) |
| `suppress_output` | bool | `True` | Suppress console output |
| `seed` | int | `0` | Random seed |
| `time_limit` | float | `1.0` | Time limit in seconds |
| `cluster_upperbound` | int | `0` | Max cluster size (0 = no limit) |

Returns a tuple `(modularity, clustering)` where `modularity` is a float in [-1, 1] and `clustering` is a list of cluster IDs for each node.

Release Notes
=====

### v1.2
- Added Python interface (`pip install vieclus`) with pybind11 bindings
- Added `vieclus_graph` helper class for easy graph construction (same interface as KaHIP)
- Added `vieclus.cluster()` function
- Added PyPI packaging with scikit-build-core
- Added GitHub Actions CI and automated PyPI publishing
- Added NOMPI compilation support

### v1.1
- Added cmake build system
- Added option to compile without MPI support

### v1.0
- Initial release of the memetic graph clustering algorithm

Licence
=====
The program is licenced under MIT licence.
If you publish results using our algorithms, please acknowledge our work by quoting the following paper:

```
@inproceedings{BiedermannHSS18,
  author    = {Biedermann, Sonja and Henzinger, Monika and Schulz, Christian and Schuster, Bernhard},
  title     = {{Memetic Graph Clustering}},
  booktitle = {17th International Symposium on Experimental Algorithms (SEA 2018)},
  series    = {LIPIcs},
  volume    = {103},
  pages     = {3:1--3:15},
  publisher = {Schloss Dagstuhl -- Leibniz-Zentrum f{\"{u}}r Informatik},
  year      = {2018},
  doi       = {10.4230/LIPIcs.SEA.2018.3}
}
```

