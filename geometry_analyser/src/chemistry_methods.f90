MODULE chemistry_methods
  USE types,  ONLY: atom, voxel, indexArray, neighbour, connection_mesh
  IMPLICIT NONE
  PRIVATE
  PUBLIC  ::  overlap
  PUBLIC  ::  electronegative
  PUBLIC  ::  determine_neighbours
  PUBLIC  ::  determine_bond_angles
CONTAINS
  FUNCTION determine_neighbours(target_atom,crystal,atomNumber)
    INTEGER(8), INTENT(IN)  ::  target_atom
    INTEGER(8), INTENT(IN)  ::  atomNumber
    TYPE(atom), INTENT(IN)  ::  crystal(3,3,3,atomNumber) 
    TYPE(atom)  ::  curr_neighbour
    TYPE(neighbour), ALLOCATABLE ::  determine_neighbours(:)
    TYPE(atom)  ::  pov_atom
    INTEGER(8)  ::  l, k, j, i, img1,img2, img3
    REAL(8) ::  r, min_coord_r
    INTEGER(8)  ::  neighbourNum, index_target(4)
    TYPE(indexArray), ALLOCATABLE ::  neighbourList(:)
    !INTEGER(8), POINTER ::  atom_indices
    ALLOCATE(neighbourList(0))
    pov_atom=crystal(2,2,2,target_atom)
    DO l= 1,3
      img3=l-2
      DO k=1,3
        img2=k-2
        DO j=1,3
          img1=j-2
          DO i=1, atomNumber
            IF(target_atom .EQ. i) THEN
              IF((l .EQ. 2) .AND. (k .EQ. 2) .AND. (j .EQ. 2)) THEN
                CYCLE
              END IF
            END IF
            r=separation(pov_atom%cartesian,crystal(j,k,l,i)%cartesian)
            IF(r .LE. &
& minimum_coordination_length(pov_atom%element, &
& crystal(l,k,j,i)%element) ) THEN
              ! so this is close enough to be counted as a neighbour <=
              !curr_neighbour crystal(l,k,j,i)
              index_target=(/ j, k, l, i /)
              CALL add_item(neighbourList,index_target)
            END IF
          END DO
        END DO
      END DO
    END DO
  neighbourNum=size(neighbourList)
  ALLOCATE(determine_neighbours(neighbourNum))
  DO i=1,neighbourNum
    !index_target=(/ ( neighbourList(i)%indices(j), j, 4 ) /)
! this is j,k,l,i or a, b, c, N
    determine_neighbours(i)%neighbour=crystal(neighbourList(i)%indices(1), &
& neighbourList(i)%indices(2), &
& neighbourList(i)%indices(3), &
& neighbourList(i)%indices(4) )
    determine_neighbours(i)%atom_num=neighbourList(i)%indices(4)
    determine_neighbours(i)%image=(/ neighbourList(i)%indices(1), &
& neighbourList(i)%indices(2), neighbourList(i)%indices(3) /)
    determine_neighbours(i)%separation=separation(pov_atom%cartesian, &
& determine_neighbours(i)%neighbour%cartesian )
    DO j=1,3
      determine_neighbours(i)%relative_vec(j)= &
& ( determine_neighbours(i)%neighbour%cartesian(j) - &
& pov_atom%cartesian(j) )
    END DO
  END DO
! so for each atom there is number of neighborus, a list of neighbours, and a
! list of bond lengths
    
  END FUNCTION determine_neighbours
  FUNCTION determine_bond_angles(connection_node)
    TYPE(connection_mesh), INTENT(IN) ::  connection_node
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8)  ::  i,j, neighbourNum
    REAL(8), ALLOCATABLE  ::  determine_bond_angles(:,:)
    REAL(8) ::  arr1(3),arr2(3),modArr1,modArr2
    neighbourNum=connection_node%neighbourNum
    neighbour_list=connection_node%neighbour_list
    ALLOCATE(determine_bond_angles(neighbourNum,neighbourNum))
    DO i=1, neighbourNum
      DO j=1, neighbourNum
        IF (i==j) CYCLE
        arr1=neighbour_list(i)%relative_vec
        arr2=neighbour_list(j)%relative_vec
        modArr1=neighbour_list(i)%separation
        modArr2=neighbour_list(j)%separation
        determine_bond_angles(i,j)=determine_angle(arr1,arr2,modArr1,modArr2)
      END DO
    END DO
  END FUNCTION determine_bond_angles
  LOGICAL FUNCTION overlap(Pvoxel,element)
    TYPE(voxel), INTENT(IN) ::  Pvoxel
    CHARACTER(len=2), INTENT(IN)  ::  element
    REAL(8) ::  sphere_to_sphere
    sphere_to_sphere=determine_vdw_radii(element)+1.0691
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
  FUNCTION minimum_coordination_length(element1,element2)
    CHARACTER(len=2), INTENT(IN)  ::  element1, element2
    INTEGER(8)  ::  ele1_indx
    REAL(8) ::  ele_row(14)
    REAL(8) ::  minimum_coordination_length
    SELECT CASE(element1)
      CASE ('C')
        ele1_indx=1
      CASE ('N')
        ele1_indx=2
      CASE ('O')
        ele1_indx=3
      CASE ('H')
        ele1_indx=4
      CASE ('Co')
        ele1_indx=5
      CASE ('Cu')
        ele1_indx=6
      CASE ('Pd')
        ele1_indx=7
      CASE ('Fe')
        ele1_indx=8
      CASE ('Ca')
        ele1_indx=9
      CASE ('Zn')
        ele1_indx=10
      CASE ('XX')
        ele1_indx=11
      CASE ('Si')
        ele1_indx=12
      CASE ('Al')
        ele1_indx=13
      CASE ('Na')
        ele1_indx=14
      CASE DEFAULT
        ele1_indx=9
    END SELECT
    SELECT CASE(element2)
      CASE ('C')
        ele_row=(/ 1.6, 1.5, 1.5, 1.12, 2.4, 2.1, 2.2, 2.4, 2.6, 2.2, &
& -1.0, 1.9, 2.1, 3.0 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('N')
        ele_row=(/ 1.5, 1.7, 1.5, 1.1, 2.1, 2.05, 2.1, 2.2, 2.5, 2.1, &
& -1.0, 1.0, 2.0, 2.6 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('O')
        ele_row=(/ 1.5, 1.5, 1.55, 1.1, 2.25, 2.4, 2.4, 2.2, 2.6, 2.6, &
& -1.0, 1.7, 2.05, 2.8 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('H')
        ele_row=(/ 1.12, 1.1, 1.1, 0.8, 1.5, 1.7, 1.7, 1.6, 2.3, 1.6, &
& -1.0, 1.5, 1.85, 1.9 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Co')
        ele_row=(/ 2.4, 2.1, 2.25, 1.5, 2.6, 2.5, 2.7, 2.6, 3.4, 2.5, &
& -1.0, 2.4, 2.5, 2.8 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Cu')
        ele_row=(/ 2.1, 2.05, 2.4, 1.7, 2.5, 2.7, 2.7, 2.6, 3.4, 2.6, &
& -1.0, 2.45, 2.6, 2.8 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Pd')
        ele_row=(/ 2.2, 2.1, 2.4, 1.7, 2.7, 2.7, 2.8, 2.8, 3.4, 2.7, &
& -1.0, 2.5, 2.7, 2.9 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Fe')
        ele_row=(/ 2.4, 2.2, 2.2, 1.6, 2.6, 2.6, 2.8, 2.8, 3.4, 2.8, &
& -1.0, 2.6, 2.7, 2.9 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Ca')
        ele_row=(/ 2.6, 2.5, 2.6, 2.3, 3.4, 3.4, 3.4, 3.4, 3.3, 3.4, &
& -1.0, 3.0, 3.2, 3.1 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Zn')
        ele_row=(/ 2.2, 2.1, 2.6, 1.6, 2.5, 2.6, 2.7, 2.8, 3.4, 2.5, &
& -1.0, 2.6, 2.8, 2.9 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('XX')
        minimum_coordination_length=-1.0
      CASE ('Si')
        ele_row=(/ 1.9, 1.0, 1.7, 1.5, 2.4, 2.45, 2.5, 2.6, 3.0, 2.6, &
& -1.0, 2.5, 2.6, 3.0 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Al')
        ele_row=(/ 2.1, 2.0, 2.05, 1.85, 2.5, 2.6, 2.7, 2.7, 3.2, 2.8, &
& -1.0, 2.6, 2.6, 3.2 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE ('Na')
        ele_row=(/ 3.0, 2.6, 2.8, 1.9, 2.8, 2.8, 2.9, 2.9, 3.1, 2.9, &
& -1.0, 3.0, 3.2, 3.1 /)
        minimum_coordination_length=ele_row(ele1_indx)
      CASE DEFAULT
        minimum_coordination_length=2.0 
    END SELECT
  END FUNCTION
  SUBROUTINE add_item(array,item)
    TYPE(indexArray), INTENT(INOUT), ALLOCATABLE  ::  array(:)
    INTEGER(8), INTENT(IN)  ::  item(4)
    TYPE(indexArray), ALLOCATABLE ::  temp(:)
    INTEGER(8)  ::  n, i
    n=size(array) 
    ALLOCATE(temp(n+1))
    DO i=1,n
      temp(i)%indices= array(i)%indices
    END DO
    n=n+1
    temp(n)%indices=item
    DEALLOCATE(array)
    ALLOCATE(array(n))
    DO i=1,n
      array(i)%indices=temp(i)%indices
    END DO
    DEALLOCATE(temp)
  END SUBROUTINE
  FUNCTION determine_angle(arr1,arr2,modArr1,modArr2)
    REAL(8), INTENT(IN) ::  arr1(3),arr2(3),modArr1,modArr2
    INTEGER(8)  ::  i
    REAL(8) ::  determine_angle
    REAL(8) ::  pi
    pi= ATAN(1.0) * 4
    determine_angle=0
    DO i=1,3
      determine_angle=determine_angle+arr1(i)*arr2(i)
    END DO 
    determine_angle=determine_angle/(modArr1*modArr2)
    determine_angle=ACOS(determine_angle)*180/pi
  END FUNCTION
END MODULE chemistry_methods

