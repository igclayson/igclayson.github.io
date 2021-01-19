MODULE types
  IMPLICIT NONE
  PUBLIC :: mesh_type
  PUBLIC :: atom
  PUBLIC :: voxel
  TYPE mesh_type
    REAL(8) ::  coordinates(3)
    LOGICAL :: occupation(2)
! note the first value is whether it is occupied and the second is whether the
! atom in the voxel is electronegative or not (no means positve)
  END TYPE
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
END MODULE types
