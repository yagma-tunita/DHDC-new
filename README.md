# DHDC: Diagnosability, Hamiltonian Connectivity, and Throughput Analysis

Analysis and verification code for the **DHDC** (dual-hypercube-based) data center network, including its fault diagnosability, Hamiltonian connectivity, and a throughput/congestion comparison against three other server-centric networks: **HSDC**, **CSDC**, and **BCDC**.

## Key Results

### Fault Diagnosability (PMC model)

Exact closed-form bounds, each verified by exhaustive search (`proof.cpp`):

| Definition | Bound |
|---|---|
| ordinary | `t_0 = n` |
| conditional (every healthy server has at least 1 healthy neighbor) | `tc = min(2n−1, n+4)` |
| g-good-neighbor (every healthy server has at least g healthy neighbors) | `t_g = min(2n+g−2, n+3g+1)`, `1 ≤ g ≤ n−2` |

Three linear-time diagnostic algorithms are implemented: `faultisis`, `faultisisConditional`, and `faultisisG`.

### Hamiltonian Connectivity

- `DH_n` is Hamiltonian-connected for `n ≥ 4`.
- `DH_3` is bipartite, hence only Hamiltonian-laceable.
- The proof uses a column-decomposition framework (hypercube skeleton + per-column `K_{n-1} □ K_2`), a parity invariant, and a fault-tolerant Hamiltonian-path lemma for hypercubes.

### Model Comparison

Aligned by **server count N** (not by dimension n), DHDC provides:

- the **fewest ports** per server,
- the **lowest congestion** (load factor),
- the **best scale-out efficiency** (most servers per port),

at the cost of the **highest latency** (largest diameter, `4n−4`).

## Directory Layout

```
├── AlgoFault.{hpp,cpp}      # diagnostic algorithms + diagnosability formulas
├── AlgoRouting.{hpp,cpp}    # routing: DHRouting, Disjointpath, StrictDisjointpath, BFS/A*/Dijkstra
├── proof.cpp                # exhaustive verification of diagnosability
├── ham_cycle.cpp            # constructive Hamiltonian cycle
├── qpath.cpp                # hypercube Hamiltonian path
├── dhd_hamconn*.cpp         # Hamiltonian connectivity construction
├── flow_sim.{hpp,cpp}       # flow-level simulation framework
├── thp_compare.cpp          # four-model throughput curve comparison
├── measure_all.cpp          # raw metrics (per dimension n)
├── measure_scale.cpp        # metrics aligned by server count N
├── measure_nrange.cpp       # metrics over a range of N
├── throughput_compare.cpp   # throughput comparison
├── traffic_compare.cpp      # traffic / congestion comparison
├── plot_data.py             # generate comparison figures -> png/
├── plot_3d.py               # generate 3D figure -> png/
├── data/                    # CSV results
└── png/                     # generated figures
```

## Build & Run

Requires a C++17 compiler (Linux / g++) and Python 3 with `matplotlib`, `scipy`.

```bash
make            # build all executables
make csv        # run data experiments; writes CSVs into data/
python plot_data.py   # generate comparison figures into png/
python plot_3d.py     # generate the 3D figure
```
