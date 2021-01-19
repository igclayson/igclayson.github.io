! This program will take into an xyz file and the lattice parameters and
! return a,
! 1) z-matrix of the desired atoms
! 2) assignment of first-nearest neighbours
! 3) connectiity and assign geoemtry and coordination number
  PROGRAM geometry_analyser
    USE crystal_module, ONLY: generate_crystal, read_from_file
    USE chemistry_methods, ONLY: determine_neighbours, determine_bond_angles
    USE types, ONLY: atom, neighbour, connection_mesh
    USE neighbours, ONLY: catogrise
    IMPLICIT NONE
    CHARACTER(len=128) :: fileName, requestChar
    TYPE(atom), ALLOCATABLE ::  coords(:)
    INTEGER(8) :: atomNumber, atomRequest
    TYPE(atom), ALLOCATABLE ::  crystal(:,:,:,:)
    REAL(8) ::  lat_a,lat_b,lat_c
    REAL(8) ::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8)  ::  i, j, neighbourNum
    TYPE(connection_mesh), ALLOCATABLE  ::  connectivity_net(:)
    REAL(8), PARAMETER  ::  pi=2*ASIN(1.0)
    ! read from file
    CALL GET_COMMAND_ARGUMENT(1,fileName)
    CALL GET_COMMAND_ARGUMENT(2,requestChar) 
    READ(requestChar,*) atomRequest
    CALL read_from_file(fileName,coords,lat_a,lat_b,lat_c, &
    & lat_alpha,lat_beta,lat_gamma,atomNumber)
    crystal=generate_crystal(atomNumber,coords,lat_a,lat_b,lat_c, &
    & lat_alpha,lat_beta,lat_gamma)
    ALLOCATE(connectivity_net(atomNumber))
    DO i=1, atomNumber
       connectivity_net(i)%neighbour_list=determine_neighbours(i,crystal,atomNumber)
       connectivity_net(i)%neighbourNum=& 
& size(connectivity_net(i)%neighbour_list)
       connectivity_net(i)%bond_angle_list=determine_bond_angles(connectivity_net(i))
               
    END DO
    CALL catogrise(connectivity_net,crystal,atomNumber)
    IF ( atomRequest .EQ. 0 ) THEN ! i.e. print all
      DO i=1, atomNumber
        WRITE (*,*) i, crystal(2,2,2,i)%element, &
& connectivity_net(i)%local_environment
      END DO
    ELSE  
      WRITE (*,*) , crystal(2,2,2,atomRequest)%element, &
& connectivity_net(atomRequest)%local_environment, &
& connectivity_net(atomRequest)%metal_framework_coordation
      neighbourNum=connectivity_net(atomRequest)%neighbourNum
      DO i=1, neighbourNum
        WRITE (*,*) &
& connectivity_net(atomRequest)%neighbour_list(i)%atom_num
      END DO 
    END IF   
 END PROGRAM
