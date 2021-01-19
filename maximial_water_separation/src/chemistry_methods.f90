MODULE chemistry_methods
  USE types,  ONLY: atom, voxel
  IMPLICIT NONE
  PRIVATE
  PUBLIC  ::  overlap
  PUBLIC  ::  electronegative
CONTAINS
  LOGICAL FUNCTION overlap(Pvoxel,element)
    TYPE(voxel), INTENT(IN) ::  Pvoxel
    CHARACTER(len=2), INTENT(IN)  ::  element
    REAL(8) ::  sphere_to_sphere
    sphere_to_sphere=1.2*(determine_vdw_radii(element)+1.0691)
    IF (sphere_to_sphere >= Pvoxel%closest_separation) THEN
      overlap=.TRUE.
      !PRINT *,"OVERLAP with touching spheres: ", sphere_to_sphere, "closest separation: ", Pvoxel%closest_separation
    ELSE
      overlap=.FALSE.
    END IF
  END FUNCTION
! private
  FUNCTION determine_vdw_radii(element)
! returns the vdW radius of an input element
    CHARACTER(len=2), INTENT(IN)  ::  element
    REAL(8) ::  determine_vdw_radii
    SELECT CASE(element)
      CASE ('C')
        determine_vdw_radii=0.76
      CASE ('N')
        determine_vdw_radii=0.71
      CASE ('O')
        determine_vdw_radii=0.66
      CASE ('H')
        determine_vdw_radii=0.66!0.31
      CASE ('Co')
        determine_vdw_radii=1.26
      CASE ('Cu')
        determine_vdw_radii=1.32
      CASE ('Pd')
        determine_vdw_radii=1.39
      CASE ('XX')
        determine_vdw_radii=-1.0691
      CASE DEFAULT
        determine_vdw_radii=1.5
    END SELECT
    RETURN
  END FUNCTION
  FUNCTION separation(coord1, coord2)
    IMPLICIT NONE
    INTEGER ::  i
    REAL(8) ::  separation
    REAL(8), INTENT(IN) ::  coord1(3), coord2(3)
    separation=0
    DO i=1, 3
      separation=separation+(coord1(i)-coord2(i))**2
    END DO
    separation = SQRT(separation)
    RETURN
  END FUNCTION
  FUNCTION electronegative(element)
! returns a logical saying whether the atom is electronegative
    CHARACTER(len=2), INTENT(IN)  ::  element
    INTEGER(2) :: electronegative
    SELECT CASE(element)
      CASE ('C')
        electronegative=1
      CASE ('N')
        electronegative=1
      CASE ('O')
        electronegative=1
      CASE ('H')
        electronegative=-1
      CASE ('Co')
        electronegative=-1
      CASE ('Cu')
        electronegative=-1
      CASE ('Pd')
        electronegative=-1
      CASE DEFAULT
        electronegative=1
    END SELECT
    RETURN
  END FUNCTION
END MODULE chemistry_methods
