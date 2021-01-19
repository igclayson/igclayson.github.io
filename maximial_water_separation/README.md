# Hydrate system

This is a serial F90 program that will hydrate a crystal by placing one water as close as possible to 1 atom, denoted by the atom's number within the corresponding xyz file.

## Function

This program takes a crystal structure and adds a water as far as possible to an atom of chose in a space that can safely accommodate a H2O, defined by the atom's covalent radius. 

## Input setup

The four cmd-line arguments taken are:
1. The filename, format described below.                                                                                               
2. The target atom to hydrate.
3. The fineness of the voxel grid for the unit cell to discretised to (where this will be a (N,N,N) grid with equal fractional steps).
4. The name of the output file the water is written to

This program takes a single file input, given as the first cmd-line argument, that has the format of an xyz file:
1. Number of atoms within file.
2. The lattice parameters as: `a b c alpha beta gamma`.
3. List of atoms where XX is coded to ghost atoms with a size of 0 A.

## Output

Output will consist of 1 file, an xyz of just the water positions.

## Tips and tricks

* The intel compilers are far faster than the GNU ones for intel systems (i.e., login nodes on UCL and ARCHER).
* The fineness of the grid can affected results with small-pore zeolites occasionally producing worse results if the grid is too fine due to small rings being sampled.
* Complex hydrogen bond networks are not well preserved so this code is not appropriate if the system has low stability dependent on the hydrogen-bonding network.
* The code will error if a `water.xyz` file is present as the closest reasonable water is written to this file.
* The code can be compiled with `ifort types.f90 crystal_module.f90 chemistry_methods.f90 voxel_grid.f90 maximial_water_separation.f90 -o ../bin/maximial_water_separation`

## Compiling

GNU:
`gfortran types.f90 crystal_module.f90 chemistry_methods.f90 voxel_grid.f90 maximial_water_separation.f90 -o ../bin/maximial_water_separation.x`
Intel:
`ifort types.f90 crystal_module.f90 chemistry_methods.f90 voxel_grid.f90 maximial_water_separation.f90 -o ../bin/maximial_water_separation.x`
