#!/bin/sh
#
# Method 1: setup for a fully dual-topology side chain residue
#

basedir=leap


tleap -f - <<_EOF
# load the AMBER force fields
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

loadamberparams $basedir/CMP.frcmod

loadoff $basedir/CMP.lib

# load the coordinates and create the systems
ligand = loadpdb $basedir/lig.pdb
m1 = loadpdb $basedir/WT.pdb
m2 = loadpdb $basedir/singleA.pdb
w = loadpdb $basedir/water_ions.pdb

protein = combine {m1 m2 w}
complex = combine {m1 m2 ligand w}

set default nocenter on

# create protein in solution
setBox protein vdw {109. 109. 109.}
savepdb protein protein.pdb
saveamberparm protein protein.parm7 protein.rst7

# create complex in solution
setBox complex vdw {109. 109. 109.}
savepdb complex complex.pdb
saveamberparm complex complex.parm7 complex.rst7

quit
_EOF
