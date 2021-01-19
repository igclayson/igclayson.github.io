#!/bin/bash
# takes 1) name of file
function crystal_to_poscar() {
  file=$1
  atom_num=$(sed '10 !d' $file )
  echo "poscar" > $file.POSCAR.vasp
  echo "1" >> $file.POSCAR.vasp
  head -n 5 $file | sed -e '1 d' -e '$ d' >> $file.POSCAR.vasp
  atom_list=( $( tail -n $((atom_num+1)) $file | sed '$ d' | awk '{ print $1}' | uniq ) )
  echo "${atom_list[@]}" >> $file.POSCAR.vasp
  for atom in ${atom_list[@]}; do
    number=$( tail -n $((atom_num+1)) $file | sed '$ d' | awk '{ print " "$1" "}' | grep " $atom " | wc -l )
    printf "$number " >> ${file}.POSCAR.vasp
  done; printf "\n" >> ${file}.POSCAR.vasp
  echo "Cartesian" >> ${file}.POSCAR.vasp
  for  atom in ${atom_list[@]}; do
    tail -n $((atom_num+1)) $file | sed '$ d' | grep " $atom "
    done | awk '{print $2 " " $3 " " $4}' >> ${file}.POSCAR.vasp
}
