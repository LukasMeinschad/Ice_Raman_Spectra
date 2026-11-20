# Ice XII: H-bond topologies for Raman calculations (CRYSTAL 17, PBE0-D3 / POB-TZVP)

## Source structure (`ice_xii.cif`)
I-42d (No. 122), a = 8.304, c = 4.024 A, O1 4a (0,0,0), O2 8d (0.3643, 1/4, 1/8), three D sites
with occupancy 1/2 (fully H-disordered). Cell and O positions are those of Lobban, Finney & Kuhs,
Nature 391, 268 (1998) (D2O, neutron diffraction). 12 H2O per cell, 23.12 A^3 per molecule,
1.294 g/cm3 for H2O; the cell volume (277.48 A^3) equals that of the recovered ice XIV cell
(277.45 A^3, `../Ice_XIV/14_XIV_P212121_1.cif`), so both phases are calculated at the same density.

`ice_xii_genice.cif` (`genice2 --rep 1 1 1 12`): 1x1x2 cell, 7.996 x 7.996 x 7.794 A, 24 H2O.
GenIce's ice 12 lattice is uniformly contracted by ~3.5 % (a -3.7 %, c -3.2 %, volume -10 %) with
respect to the experimental cell; the fractional coordinates are the same lattice. Its H
arrangement is topology **2607** of the 1x1x2 enumeration (P2_1, multiplicity 16, net dipole
|P| = 5.96) (`tools/identify_topology.py XII --cell 1x1x2 ice_xii_genice.cif`).

## Enumeration (`tools/enumerate_ice.py XII`)
| cell | H2O | ice-rule configurations | symmetry-distinct topologies |
|---|---|---|---|
| unit (8.304 x 8.304 x 4.024) | 12 | 692 | **62** (`enumeration_unit/`) |
| 1x1x2 (GenIce cell) | 24 | 88116 | 2891 (`enumeration_1x1x2/`, structure files only for 2607) |

Space groups of the 62 unit-cell topologies: P1 33, P2_1 18, P2_12_12_1 4, Cc 3, Fdd2 2, C2 1, P2 1.

**Relation to ice XIV:** the three ordered ice XIV models of `../Ice_XIV` are ice XII topologies
28 (model 1, P2_12_12_1), 41 (model 2, P2_1) and 31 (model 3, P2_12_12_1) - ice XIV is ice XII with
the H order of topology 28 dominating (~60 %), in an orthorhombically distorted cell.

## Raman jobs (`pbe0_d3_raman/`, 40 cores each)
The 10 topologies with at least 4 symmetry operations (as for ice VI, the low-symmetry ones are
too expensive at this level of theory):

| id | space group | multiplicity | net dipole |
|---|---|---|---|
| 9, 14 | Fdd2 | 2 | yes |
| 11 | C2 | 4 | yes |
| 28, 31 | P2_12_12_1 | 4 | 0 (= ice XIV models 1 / 3, here in the tetragonal ice XII cell) |
| 29, 32 | P2_12_12_1 | 4 | 0 |
| 52, 57, 58 | Cc | 4 | yes |

Same settings as ice XIV: ATOMONLY optimisation in the experimental cell, then FREQCALC with
INTENS + INTRAMAN + INTCPHF from the optimised geometry (EXTERNAL). k-points: SHRINK 0 8 / 8 4 4
(8 along the 4.0 A axis) for the primitive P2_12_12_1 cells, isotropic SHRINK 6 12 for the
centred cells (Fdd2, C2, Cc; primitive vectors 6.2-8.3 A).

```
cd Ice_XII/pbe0_d3_raman
./submit_all.sh              # all 10 opt jobs, each freq job afterok its opt job
./submit_all.sh 28 31 9      # selected ids
```
`pbe0_d3_raman_1x1x2/config_02607/`: the GenIce H arrangement in the experimental (1x1x2) cell,
72 atoms, P2_1 - roughly 5-10x the cost of a unit-cell job; optional.
