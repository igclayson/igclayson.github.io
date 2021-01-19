PROGRAM hydrate_system
    USE crystal_module, ONLY: generate_crystal, read_from_file, crys_separation
    USE voxel_grid, ONLY: generate_set_of_voxel, hydrate_voxel
    USE types, ONLY: atom, voxel
    IMPLICIT NONE
    CHARACTER(len=128) :: fileName, grid_char, target_char, outputFile
    TYPE(atom), ALLOCATABLE ::  coords(:)
    TYPE(atom)  ::  water(3)
    INTEGER(8) :: atomNumber
    TYPE(atom), ALLOCATABLE ::  crystal(:,:,:,:)
    REAL(8) ::  lat_a,lat_b,lat_c
    REAL(8) ::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8)  ::  grid_finess, target_atom
    TYPE(voxel), ALLOCATABLE ::  voxel_set(:,:,:)
    LOGICAL ::  file_exists
    REAL(8) ::  r_val
    ! read from file
    CALL GET_COMMAND_ARGUMENT(1,fileName)
    CALL GET_COMMAND_ARGUMENT(2,target_char)
    CALL GET_COMMAND_ARGUMENT(3,grid_char)
    CALL GET_COMMAND_ARGUMENT(4,outputFile)
    READ(target_char,*) target_atom
    READ(grid_char,*) grid_finess
    CALL read_from_file(fileName,coords,lat_a,lat_b,lat_c, &
    & lat_alpha,lat_beta,lat_gamma,atomNumber)
    crystal=generate_crystal(atomNumber,coords,lat_a,lat_b,lat_c, &
    & lat_alpha,lat_beta,lat_gamma)
    voxel_set=generate_set_of_voxel(grid_finess,lat_a,lat_b,lat_c,&
    &lat_alpha,lat_beta,lat_gamma)
    water=hydrate_voxel(voxel_set,grid_finess, &
    & crystal,target_atom,atomNumber)
    r_val=crys_separation(crystal,target_atom, water(1))
    !READ(11,*)
    !DO i=1,atomNumber
    !  READ(11,*)
    !END DO
    INQUIRE(FILE=outputFile, EXIST=file_exists)
    IF (file_exists) THEN
!      OPEN(12, FILE = 'water.xyz', STATUS="old", POSITION="append", &
!& ACTION="write")
      WRITE(*,*) "There is already an ouput file with the specified name -- ERROR"
      WRITE(*,*) outputFile
      STOP 10
    ELSE
      OPEN(12, FILE = outputFile, STATUS= 'new', ACTION="write")
      WRITE(12,*) 3
      WRITE(12,*) r_val
    END IF
    WRITE(12,*) "O", water(1)%cartesian(1), water(1)%cartesian(2), &
& water(1)%cartesian(3)
    WRITE(12,*) "H", water(2)%cartesian(1), water(2)%cartesian(2), &
& water(2)%cartesian(3)
    WRITE(12,*) "H", water(3)%cartesian(1), water(3)%cartesian(2), &
& water(3)%cartesian(3)
   CLOSE(12,STATUS='KEEP')
END PROGRAM
