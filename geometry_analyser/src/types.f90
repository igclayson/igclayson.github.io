MODULE types
  IMPLICIT NONE
  PUBLIC :: atom
  PUBLIC :: voxel
  PUBLIC :: neighbour
  PUBLIC :: connection_mesh
! note the first value is whether it is occupied and the second is whether the
! atom in the voxel is electronegative or not (no means positve)
  TYPE atom
    CHARACTER(2) :: element
    REAL(8) :: cartesian(3)
    REAL(8) :: fractional(3)
  END TYPE
  TYPE voxel
    REAL(8) ::  center(3)
    REAL(8) ::  target_separation ! this the separation from the target atom of
!interest being hydrated
    INTEGER(2) :: target_image(3)
    INTEGER(8) ::  closest_atom ! this is the atom number which the the water
!will be closest to
    INTEGER(2) :: closest_image(3) ! this is the image which the closest atom is
!in 
    REAL(8) ::  closest_separation ! the r from the closest atom
  END TYPE
  TYPE neighbour 
    TYPE(atom)  ::  neighbour
    INTEGER(8)  ::  atom_num
    INTEGER(8)  ::  image(3)
    REAL(8) ::  separation
    REAL(8) ::  relative_vec(3)
    INTEGER(8)  ::  neighbour_type
  END TYPE
  TYPE indexArray
! note this is a pointer to the atom indices
    INTEGER(8)  ::  indices(4)
  END TYPE
  TYPE connection_mesh
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    REAL(8), ALLOCATABLE ::  bond_angle_list(:,:)
    INTEGER(8)  ::  neighbourNum
    INTEGER(8)  ::  local_environment
    INTEGER(8)  ::  metal_framework_coordation
! the above is only assigned if the metal is coordinated to the zeo
  END TYPE
END MODULE types
