# Geometry analyser
This is a geometry analyser script that is designed to look at the local enviornments of atoms in predominately zeolitic systems. Writes to `stdout`

## Inputs

This program takes,
1. The file name of the crystal you wish to analyse. In the form of an xyz with the 2nd line filled with the lattice parameters (a, b, c, alpha, beta, and gamma) only
2. The atom number of the atom you wish to inspect where 0 will inspect all atoms 

## Compiling

To compile with a fortran compile, use the following sequence
` types.f90 crystal_module.f90 neighbours.f90 chemistry_methods.f90 geometry_analyser.f90 -o geometry_analyser.x `
