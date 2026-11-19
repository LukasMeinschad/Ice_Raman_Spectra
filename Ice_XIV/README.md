# Ice XIV: ordered models for Raman calculations (CRYSTAL 17, PBE0-D3 / POB-TZVP)

Source structure: `14_XIV_P212121_1.cif` : D2O ice XIV,
P2_12_12_1, a = 8.3499, b = 8.1391, c = 4.0825 A, 3 O in the asymmetric unit (12 D2O / cell).

## Partial occupancies

| D site | occupancy | bond (donor -> acceptor) |
|---|---|---|
| D6, D9, D11, D15 | 1.0 | fully ordered |
| D12 / D4 | 0.593 / 0.407 | the two ends of one O1-O2 bond |
| D5 / D13 | 0.620 / 0.380 | the two ends of a second O1-O2 bond |

So 16 of the 24 bonds per cell are fixed, and 8 (2 bond types x 4) are partially ordered.
`tools/enumerate_partial_order.py` enumerates all ice-rule arrangements of the
partially ordered bonds with the fixed ones kept. In the unit cell the ice rules couple the
two partial bond types; 4 configurations remain, giving **3 symmetry-distinct ordered models**:

| model | space group | occupancy of the majority sites D12 / D5 | multiplicity | net dipole |
|---|---|---|---|---|
| 1 | P2_12_12_1 | 1.0 / 1.0 (ordered towards the experimental majority) | 1 | 0 |
| 2 | P2_1 | 0.5 / 0.5 | 2 | yes |
| 3 | P2_12_12_1 | 0.0 / 0.0 (ordered towards the minority sites) | 1 | 0 |

All models use the experimental O and D positions of the occupied sites (no idealisation).
The structures are in `ordered_models/` (CIF P1 + symmetry-adapted, CRYSTAL geometries,
summary.csv). A 1x1x2 supercell (24 H2O, 6 models, occupancies 0/0.25/0.5/0.75/1) can be
generated with `--supercell 1 1 2`.

### Weighting to the experimental partial order (~0.6)
If the partially ordered bonds are treated as independent (each with probability f = 0.6 on
the majority end), the model weights follow p_i ~ g_i * lambda^(majority ends) with
lambda = f/(1-f) = 1.5:

| model | weight |
|---|---|
| 1 | 0.36 |
| 2 | 0.48 |
| 3 | 0.16 |

which reproduces an average occupancy of 0.60. Alternatively, weights can come from the
calculated energies (Boltzmann).

## GenIce structure (`xiv_genice.cif`) and the 1x1x2 supercell

`genice2 --rep 1 1 1 14` gives a 1x1x2 cell (8.3499 x 8.1391 x 8.165 A, 24 H2O). GenIce's ice 14
lattice keeps the fully ordered bonds of ice XIV (all 32 fully occupied D sites are occupied) and
randomises only the partially ordered O1-O2 bonds; with its default depolarisation the result has
zero net dipole and occupancies D12/D4 = D5/D13 = 0.5/0.5.

The 1x1x2 supercell has 16 ice-rule configurations = 6 ordered models (`ordered_models_1x1x2/`,
`enumerate_partial_order.py --supercell 1 1 2`):

| 1x1x2 model | space group | occupancy D12/D5 | multiplicity | net dipole | note |
|---|---|---|---|---|---|
| 1 | P2_1 | 0.50 | 4 | 0 | **= the GenIce structure** (`identify_topology.py XIV --cell 1x1x2`) |
| 2 | P2_1 | 0.75 | 4 | yes | new in the supercell |
| 3 | P2_12_12_1 | 1.00 | 1 | 0 | = unit-cell model 1 |
| 4 | P2_1 | 0.25 | 4 | yes | new in the supercell |
| 5 | P2_1 | 0.50 | 2 | yes | = unit-cell model 2 |
| 6 | P2_12_12_1 | 0.00 | 1 | 0 | = unit-cell model 3 |

`pbe0_d3_raman_1x1x2/config_00001/` contains the opt + freq job for the GenIce structure (1x1x2
model 1, built from the experimental O/D positions; it has the same H-bond topology as GenIce's
coordinates, so the ATOMONLY optimisation leads to the same structure). Same settings as below,
except SHRINK 0 8 / 4 4 4 (all axes ~8 A). 72 atoms, P2_1: expect roughly 8x the cost of a
unit-cell model per SCF.

```bash
cd Ice_XIV/pbe0_d3_raman_1x1x2
./submit_all.sh
```

## Jobs (`pbe0_d3_raman/`)

* opt: OPTGEOM ATOMONLY (experimental cell kept), TOLDEG 1e-5, TOLDEX 4e-5, TOLDEE 10, MAXCYCLE 200;
  PBE0-D3, POB-TZVP, XLGRID, SHRINK 0 8 / 8 4 4 (8 k-points along the short 4.08 A axis,
  which is **a** in the standardised P2_12_12_1 setting), TOLINTEG 10 10 10 10 18, FMIXING 70, DIIS
* freq: EXTERNAL (optimised geometry incl. symmetry) + FREQCALC / INTENS / INTRAMAN / INTCPHF ... END
  (IR + Raman intensities via CPHF/CPKS; PBE0 is supported by CPHF/CPKS, unlike HSEsol in CRYSTAL 17)
* 40 cores, CPU_rune, CRYSTAL 17 (`rune/crystal/17_v1.0.2`, `runmpi17`), as for the ice Ih / VI runs

```bash
cd Ice_XIV/pbe0_d3_raman
./submit_all.sh          # 3 opt jobs + dependent freq jobs (3 x 40 cores at a time)
```

**Notes**
1. PBE0 is among the functionals supported by CPHF/CPKS in CRYSTAL 17 (HSEsol is not), so
   `INTRAMAN` + `INTCPHF` should run. Optional: `RAMANEXP` (temperature, laser wavelength) inside
   FREQCALC to match your measurement.
2. **Isotopes:** the calculation uses H masses (H2O ice). The CIF is a D2O structure. If the
   measured sample is D2O, add D masses (FREQCALC keyword `ISOTOPES`) to the frequency inputs.

## Regenerating

```bash
cd tools
python enumerate_partial_order.py ../Ice_XIV/14_XIV_P212121_1.cif --name "ice XIV" -o ../Ice_XIV/ordered_models
python make_crystal_jobs.py ../Ice_XIV/ordered_models \
    --opt-template ../Ice_XIV/templates/opt_pbe0d3.d12 \
    --freq-template ../Ice_XIV/templates/freq_pbe0d3_raman.d12 \
    --sjob-template ../Ice_XIV/templates/crystal17_rune.sjob \
    --prefix xiv --jobname XIV --ncores 40 --out ../Ice_XIV/pbe0_d3_raman
python job_status.py ../Ice_XIV/pbe0_d3_raman        # progress / time estimates
```
