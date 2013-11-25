#!/bin/bash


atm=3030

for resi in 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314
do
    echo "ATOM   "$atm"  N   HIS A "$resi"      19.613   8.624  -3.616  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  CA  HIS A "$resi"      19.935   7.229  -3.928  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  C   HIS A "$resi"      21.358   7.005  -4.479  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  O   HIS A "$resi"      21.994   5.997  -4.154  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  CB  HIS A "$resi"      18.925   6.644  -4.914  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  CG  HIS A "$resi"      19.256   5.247  -5.323  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  ND1 HIS A "$resi"      19.145   4.182  -4.460  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  CD2 HIS A "$resi"      19.773   4.752  -6.470  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  CE1 HIS A "$resi"      19.549   3.082  -5.069  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  NE2 HIS A "$resi"      19.931   3.401  -6.291  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  H   HIS A "$resi"      18.878   9.086  -4.132  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  HA  HIS A "$resi"      19.897   6.634  -3.016  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm" 1HB  HIS A "$resi"      17.932   6.650  -4.464  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm" 2HB  HIS A "$resi"      18.884   7.268  -5.806  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  HD2 HIS A "$resi"      20.002   5.316  -7.374  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  HE1 HIS A "$resi"      19.564   2.082  -4.636  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
    echo "ATOM   "$atm"  HE2 HIS A "$resi"      20.283   2.759  -6.986  1.00  0.00              " >> 3ANS_combo.pdb
    atm=$((atm+1))
done
