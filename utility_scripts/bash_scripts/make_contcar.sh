#!/bin/bash
# This is a small bash script written by Ivan which takes 2 positional variables: 1) The main file "stem" or PROJECT_NAME
# and 2) the number of atoms
# This is designed to work with CP2K with the PROJECT-pos-1.xyz and PROJECT.restart files being present in the directory
function make_contcar() {
  file="$1"
  atom_num=$( head -n 1 ${file}-pos-1.xyz )
  atoms_list=( $( tail -n $atom_num ${file}-pos-1.xyz  | sort | awk '{print $1}' | uniq) )
  echo "${file}" > ${file}.POSCAR.vasp
  echo "1.0" >> ${file}.POSCAR.vasp
  grep -A 4 "&SUBSYS" ${file}.restart | tail -n 3 | awk '{print $2 " " $3 " " $4}' >> ${file}.POSCAR.vasp
  echo ${atoms_list[@]} >> ${file}.POSCAR.vasp
  for atom in ${atoms_list[@]}; do
    number=$( tail -n $atom_num ${file}-pos-1.xyz | grep "$atom " | wc -l )
    printf "$number " >> ${file}.POSCAR.vasp
  done; printf "\n" >> ${file}.POSCAR.vasp
  echo "Cartesian" >> ${file}.POSCAR.vasp
  for atom in ${atoms_list[@]}; do
    tail -n $atom_num ${file}-pos-1.xyz | grep "$atom "
  done | awk '{print $2 " " $3 " " $4}' >> ${file}.POSCAR.vasp
}
