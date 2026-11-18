# Ice Raman spectra

Raman (and IR) spectra of ice phases from CRYSTAL calculations on ensembles of proton-ordered
structures, for comparison with measured spectra.

| folder | content |
|---|---|
| `Ice_XIV/` | ice XIV (partially ordered, P2_12_12_1): ordered models derived from the experimental partial D occupancies, PBE0-D3 / POB-TZVP opt + IR/Raman frequency jobs (CRYSTAL 17, rune, 40 cores). See `Ice_XIV/README.md`. |
| `tools/` | Python tools (numpy, spglib): enumeration of H-bond topologies (`enumerate_ice.py`, `iceenum/`), ordered models from partial occupancies (`enumerate_partial_order.py`), CRYSTAL job generation (`make_crystal_jobs.py`), job progress (`job_status.py`), energies and spectra analysis (`analyze_energies.py`, `analyze_spectra.py`), identification of GenIce/P1 structures (`identify_topology.py`). Copied from the IceV_IceIh repository (`graph_enumeration/`), where the method, tests and validation are documented. |

Requirements: `pip install numpy spglib matplotlib`
