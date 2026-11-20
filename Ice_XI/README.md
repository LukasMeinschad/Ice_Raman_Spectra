# Ice XI: method test (fully ordered, CRYSTAL 17, PBE0-D3 / POB-TZVP)

Ice XI (Cmc2_1) is the fully proton-ordered form of ice Ih, so a single calculation describes the
whole crystal (no disorder averaging) - a clean test of the level of theory against the Raman
spectra of ice XI (e.g. Abe & Shigenari, J. Chem. Phys. 134, 104506 / 104507 (2011)).

Structure: topology 14 of the 8-H2O orthohexagonal ice Ih cell (`enumeration_ortho/`,
`tools/enumerate_ice.py Ih --cell ortho`; 16 topologies, ice XI is the only Cmc2_1 one) =
topology 122 of the 16-H2O enumeration in the IceV_IceIh repository (HSEsol-D3 IR runs).
Cell: Rottger et al. (1994) ice Ih at 10 K in the orthohexagonal setting, a = 4.4973,
b = 7.7896, c = 7.3240 A (the orthorhombic distortion of ice XI is < 0.1 %), O at the ideal
z = 1/16, H on the ideal H-bond positions; the ATOMONLY optimisation relaxes them.

`pbe0_d3_raman/config_00014/`: same settings as ice XIV / XII (ATOMONLY opt, then FREQCALC with
INTENS + INTRAMAN + INTCPHF), SHRINK 6 12 (C-centred cell, primitive vectors 4.50 / 4.50 / 7.32 A),
40 cores.

    cd Ice_XI/pbe0_d3_raman && ./submit_all.sh
