MODULE neighbours
  USE types,  ONLY: atom, neighbour, connection_mesh
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: catogrise
CONTAINS
  SUBROUTINE catogrise(connectivity_net,crystal,atomNumber)
! this is going to be some decision-esaque tree to determine the enviorment of
! the neighbour (i.e., where it is a water of not
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    INTEGER(8), INTENT(IN)  ::  atomNumber
    TYPE(atom), INTENT(IN)  ::  crystal(3,3,3,atomNumber)
    TYPE(neighbour), ALLOCATABLE  ::  neighbour_list(:)
    INTEGER(8)  ::  neighbourNum, i,j 
    INTEGER(8)  ::  metal_coordination(2)
! this is a very zeolite object and can be ignored when not detaling
! with metal migations
    ! so the connectivity is a atomNumber long list of neighbour lists
    DO i=1, atomNumber ! this does the O first as that's the most important and
! only takes the single connector node
      IF (crystal(2,2,2,i)%element == 'O' ) THEN
        connectivity_net(i)%local_environment= &
& catogrise_O(connectivity_net(i))
      END IF
    END DO
    DO i=1, atomNumber
      SELECT CASE(crystal(2,2,2,i)%element)
        CASE ('O')
          CYCLE
        CASE ('H')
          connectivity_net(i)%local_environment= &
&  catogrise_H(connectivity_net,i)
        CASE ('Si')
          connectivity_net(i)%local_environment= &
& catogrise_T_atom(connectivity_net,i)
        CASE ('Al')
          connectivity_net(i)%local_environment= &
& catogrise_T_atom(connectivity_net,i)
        CASE ('P')
          connectivity_net(i)%local_environment= &
& catogrise_T_atom(connectivity_net,i)
        CASE DEFAULT
          metal_coordination=catogrise_metals(connectivity_net,i)
          connectivity_net(i)%local_environment=metal_coordination(1)
          connectivity_net(i)%metal_framework_coordation= &
& metal_coordination(2)
      END SELECT
    END DO
    !WRITE (*,*) "start"
    !DO i=1, atomNumber
    !  !neighbourNum=size(connectivity_net(i)%neighbour_list)
    !  ALLOCATE(neighbour_list(connectivity_net(i)%neighbourNum))
    !  neighbour_list=connectivity_net(i)%neighbour_list
    !  DO j=1, connectivity_net(i)%neighbourNum
    !    IF(neighbour_list(j)%neighbour%element == 'O' ) THEN
    !      !WRITE (*,*) neighbour_list(j)%image
    !      CALL catogrise_O(neighbour_list(j),connectivity_net%bond_angle_list) ! this is a  neighbour
    !    END IF
    !  END DO
    !  !neighbourNum=size(connectivity_net(atomNumber))
    !  !DO j=1,neighbourNum
    !    
    !  !END DO
    !  DEALLOCATE(neighbour_list)
    !END DO
    !WRITE (*,*) "end.."
  END SUBROUTINE
! private
  INTEGER FUNCTION catogrise_O(connection_node)
    TYPE(connection_mesh) ::  connection_node
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    REAL(8), ALLOCATABLE ::  bond_angle_list(:,:)
    INTEGER(8)  ::  neighbourNum
    REAL(8) ::  max_min_angle(2)
! so 1 = "free" bent H2O
! 2 = trigonal planar
! 9 = T-shaped
! 18 = some undefined tri-coordinate
! 3 = square planar
! 7 = deformed square planar
! 4 = H2O linear
! 5 = hydroxide
! 8 = other monobonded structre
! 6 = atomic O or oxide
! default = 15
! If there is a T shaped 
! N=2
! 10 = linear T2O; 11 = bent T2o
! 12 = silanonl; 13= other zeolite defect
! N=3
! 14 = Bronstead acid; 15= expanded coordination shell
! 16 = adsorbed water; 17 = other adsorbed species
! this is min 1st then max
    bond_angle_list=connection_node%bond_angle_list
    neighbour_list=connection_node%neighbour_list
    neighbourNum=connection_node%neighbourNum
    max_min_angle=return_max_min_angle(bond_angle_list,neighbourNum)
    SELECT CASE(neighbourNum)
       CASE (0)
        catogrise_O=6
      CASE (1)
        IF ( H_atom_count(neighbour_list,neighbourNum) .GT. 0 ) THEN
          catogrise_O=5
        ELSE
          catogrise_O=8
        END IF
      CASE (2)
        SELECT CASE (T_atom_count(neighbour_list,neighbourNum)) ! <- T atoms
          CASE (2)
            IF (max_min_angle(2) .GE. 150.0 ) THEN
              catogrise_O=10
            ELSE
              catogrise_O=11
            END IF
          CASE (1)
            IF (H_atom_count(neighbour_list,neighbourNum) .GT. 0 ) THEN
              catogrise_O=12
            ELSE
              catogrise_O=13
            END IF
          CASE (0)
            IF (max_min_angle(2) .GE. 150.0 ) THEN
              catogrise_O=4
            ELSE
              catogrise_O=1
            END IF
        END SELECT ! <- T atoms
      CASE (3)
        SELECT CASE (T_atom_count(neighbour_list,neighbourNum)) ! <- T atoms
          CASE (2)
            IF ( H_atom_count(neighbour_list,neighbourNum) .GT. 0 ) THEN
              catogrise_O=14
            ELSE
              catogrise_O=15
            END IF
          CASE (1)
            IF ( H_atom_count(neighbour_list,neighbourNum) .GT. 1 ) THEN
              catogrise_O=16
            ELSE
              catogrise_O=17
            END IF
          CASE (0)
            IF ( (max_min_angle(1) .GE. 150.0 ) .AND. &
    &  (max_min_angle(2) .LE. 95.0 ) ) THEN
              catogrise_O=9
            ELSE
              catogrise_O=2
            END IF
          CASE DEFAULT
            catogrise_O=18
        END SELECT
      CASE (4)
        IF ( (max_min_angle(1) .GE. 95.0 ) .OR. &
&  (max_min_angle(2) .LE. 85.0 ) ) THEN
          catogrise_O=7
        ELSE
          catogrise_O=3
        END IF
      CASE DEFAULT 
        catogrise_O=15
    END SELECT
  END FUNCTION catogrise_O
  INTEGER FUNCTION catogrise_H(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum
    INTEGER(8)  ::  neighboursIndx
! 1 = proton / free HO
! 2 = silanol O; 3 = bronstead H; 6 = aqueous
! 4 = hydride    
! 5 = shared proton
    
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    SELECT CASE(neighbourNum)
      CASE (0)
        catogrise_H=1
      CASE (1)
        IF ( neighbour_list(1)%neighbour%element == 'O' ) THEN
          neighboursIndx=neighbour_list(1)%atom_num
          SELECT CASE(connectivity_net(neighboursIndx)%local_environment)
            CASE (12)
              catogrise_H=2
            CASE (14)
              catogrise_H=3
            CASE DEFAULT
              catogrise_H=6
          END SELECT
        ELSE
          catogrise_H=4
        END IF
      CASE DEFAULT
        catogrise_H=5
    END SELECT
  END FUNCTION catogrise_H
  FUNCTION catogrise_metals(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum
    INTEGER(8)  ::  neighboursIndx
    INTEGER(8)  ::  catogrise_metals(2)
    REAL(8), ALLOCATABLE ::  bond_angle_list(:,:)
    REAL(8) ::  max_min_angle(2)
! note this returns the local coordination assignment and THEN how many
! framework coordination
! 1 = closely packed (i.e. more than 7 neighbours)
! 2= octahedral, 3= distorted 6-coordinate
! 14= square pyramical; 4= trigonal bipyrmiadal 
! 5= square planar; 6= tetrahedral; 7=boat or distorted tetrahedral
! 8= trigonal planar; 9= T-shaped
! 10= linear; 11= bent
! 12 = single-bound
! 13 = atomic
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    bond_angle_list=connectivity_net(pov_index)%bond_angle_list
    max_min_angle=return_max_min_angle(bond_angle_list,neighbourNum)
    SELECT CASE(neighbourNum)
      CASE (0)
        catogrise_metals(1)=13
      CASE (1)
        catogrise_metals(1)=12
      CASE (2)
        IF (max_min_angle(2) .GE. 150.0 ) THEN
          catogrise_metals(1)=10
        ELSE
          catogrise_metals(1)=11
        END IF
      CASE (3)
        IF( max_min_angle(1) .GE. 140 ) THEN
          catogrise_metals(1)=8
        ELSE
          catogrise_metals(1)=9
        END IF 
      CASE (4)
        IF ( (max_min_angle(1) .LE. 95.0 ) .OR. &
&  (max_min_angle(2) .GE. 85.0 ) ) THEN
         catogrise_metals(1)=5 
        ELSE IF ( (max_min_angle(1) .GE. 120.0) .OR. &
& (max_min_angle(2) .LE. 80.0) ) THEN
          catogrise_metals(1)=7
        ELSE
          catogrise_metals(1)=6
        END IF
      CASE (5)
        IF(rightAngle_count(bond_angle_list,neighbourNum)>15)  THEN
          catogrise_metals(1)=14
        ELSE
          catogrise_metals(1)=4
        END IF
! so see how many 81 -> 99 angles and if > 15 then is square pyramidal
!        IF( max_min_angle(1) .GE. 140 ) THEN
!          catogrise_metals(1)=14
!        ELSE
!          catogrise_metals(1)=4
!        END IF
      CASE (6)
        IF( ( max_min_angle(1) .GE. 100.0 ) .OR. &
& ( max_min_angle(2) .LE.80 ) ) THEN
          catogrise_metals(1)=3
        ELSE
          catogrise_metals(1)=2
        END IF
      CASE DEFAULT
        catogrise_metals(1)=1
    END SELECT
    catogrise_metals(2)=count_framework_O(connectivity_net,pov_index)
  END FUNCTION catogrise_metals
  FUNCTION catogrise_T_atom(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8)  ::  pov_index
    INTEGER(8)  ::  catogrise_T_atom
    INTEGER(8)  ::  neighboursIndx, neighbourNum
! consult the spreadsheet flow diagram to determine what all the numbers mean rather than the code (as it confusing!)
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    SELECT CASE(neighbourNum)
      CASE (0)
        catogrise_T_atom=0
      CASE (1)
        catogrise_T_atom=1
      CASE (2)
        catogrise_T_atom=2
      CASE (3)
        IF (T_X_bond_present(neighbour_list,neighbourNum) .EQ. .TRUE. ) THEN
          catogrise_T_atom=399
        ELSE
          IF (monobonded_O_present(connectivity_net,pov_index) .EQ. .TRUE. ) THEN
            catogrise_T_atom=398
          ELSE
            SELECT CASE(count_framework_O_T(connectivity_net,pov_index))
              CASE (3)
                catogrise_T_atom=30
              CASE (2)
                IF ( silanol_count(connectivity_net,pov_index) > 0 ) THEN
                  catogrise_T_atom=31
                ELSE
                  catogrise_T_atom=32
                END IF
              CASE (1)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (2)
                    catogrise_T_atom=33
                  CASE (1)
                    catogrise_T_atom=34
                  CASE (0)
                    catogrise_T_atom=35
                END SELECT
              CASE (0)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (3)
                    catogrise_T_atom=36
                  CASE (2)
                    catogrise_T_atom=37
                  CASE (1)
                    catogrise_T_atom=38
                  CASE (0)
                    catogrise_T_atom=39
                END SELECT
            END SELECT
          END IF
        END IF
      CASE (4)
        IF (T_X_bond_present(neighbour_list,neighbourNum) .EQ. .TRUE. ) THEN
          catogrise_T_atom=499
        ELSE
          IF (monobonded_O_present(connectivity_net,pov_index) .EQ. .TRUE. ) THEN
            catogrise_T_atom=498 
          ELSE
            SELECT CASE (count_framework_O_T(connectivity_net,pov_index))
              CASE (4)
                catogrise_T_atom=40
              CASE (3)
                IF ( silanol_count(connectivity_net,pov_index) > 0 ) THEN
                  catogrise_T_atom=41
                ELSE
                  catogrise_T_atom=42
                END IF  
              CASE (2)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (2)
                    catogrise_T_atom=43
                  CASE (1)
                    catogrise_T_atom=44
                  CASE (0)
                    catogrise_T_atom=45
                END SELECT
              CASE (1)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (3)
                    catogrise_T_atom=46
                  CASE (2)
                    catogrise_T_atom=47
                  CASE (1)
                    catogrise_T_atom=48
                  CASE (0)
                    catogrise_T_atom=49
                END SELECT
              CASE (0)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (4)
                    catogrise_T_atom=410
                  CASE (3)
                    catogrise_T_atom=411
                  CASE (2)
                    catogrise_T_atom=414
                  CASE (1)
                    catogrise_T_atom=413
                  CASE (0)
                    catogrise_T_atom=414
                END SELECT
            END SELECT
          END IF
        END IF
      CASE (5)
        IF (T_X_bond_present(neighbour_list,neighbourNum) .EQ. .TRUE. ) THEN
          catogrise_T_atom=599
        ELSE
          IF (monobonded_O_present(connectivity_net,pov_index) .EQ. .TRUE. ) THEN
            catogrise_T_atom=598
          ELSE
            SELECT CASE (count_framework_O_T(connectivity_net,pov_index))
              CASE (5)
                catogrise_T_atom=500
              CASE (4)
                IF ( silanol_count(connectivity_net,pov_index) > 0 ) THEN
                  catogrise_T_atom=501
                ELSE
                  catogrise_T_atom=502
                END IF
              CASE (3)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (2)
                    catogrise_T_atom=503
                  CASE (1)
                    catogrise_T_atom=504
                  CASE (0)
                    catogrise_T_atom=505
                END SELECT
              CASE (2)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (3)
                    catogrise_T_atom=506
                  CASE (2)
                    catogrise_T_atom=507
                  CASE (1)
                    catogrise_T_atom=508
                  CASE (0)
                    catogrise_T_atom=509
                END SELECT
              CASE (1)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (4)
                    catogrise_T_atom=510
                  CASE (3)
                    catogrise_T_atom=511
                  CASE (2)
                    catogrise_T_atom=512
                  CASE (1)
                    catogrise_T_atom=513
                  CASE (0)
                    catogrise_T_atom=514
                END SELECT
              CASE (0)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (5)
                    catogrise_T_atom=515
                  CASE (4)
                    catogrise_T_atom=516
                  CASE (3)
                    catogrise_T_atom=517
                  CASE (2)
                    catogrise_T_atom=518
                  CASE (1)
                    catogrise_T_atom=519
                  CASE (0)
                    catogrise_T_atom=520
                END SELECT
            END SELECT
          END IF
        END IF
      CASE (6)
        IF (T_X_bond_present(neighbour_list,neighbourNum) .EQ. .TRUE. ) THEN
          catogrise_T_atom=699
        ELSE
          IF (monobonded_O_present(connectivity_net,pov_index) .EQ. .TRUE. ) THEN
            catogrise_T_atom=698
          ELSE
            SELECT CASE (count_framework_O_T(connectivity_net,pov_index))
              CASE (6)
                catogrise_T_atom=600
              CASE (5)
                IF ( silanol_count(connectivity_net,pov_index) > 0 ) THEN
                  catogrise_T_atom=601
                ELSE
                  catogrise_T_atom=602
                END IF
              CASE (4)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (2)
                    catogrise_T_atom=603
                  CASE (1)
                    catogrise_T_atom=604
                  CASE (0)
                    catogrise_T_atom=605
                END SELECT
              CASE (3)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (3)
                    catogrise_T_atom=606
                  CASE (2)
                    catogrise_T_atom=607
                  CASE (1)
                    catogrise_T_atom=608
                  CASE (0)
                    catogrise_T_atom=609
                END SELECT
              CASE (2)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (4)
                    catogrise_T_atom=610
                  CASE (3)
                    catogrise_T_atom=611
                  CASE (2)
                    catogrise_T_atom=612
                  CASE (1)
                    catogrise_T_atom=613
                  CASE (0)
                    catogrise_T_atom=614
                END SELECT
              CASE (1)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (5)
                    catogrise_T_atom=615
                  CASE (4)
                    catogrise_T_atom=616
                  CASE (3)
                    catogrise_T_atom=617
                  CASE (2)
                    catogrise_T_atom=618
                  CASE (1)
                    catogrise_T_atom=619
                  CASE (0)
                    catogrise_T_atom=620
                END SELECT
              CASE (0)
                SELECT CASE(silanol_count(connectivity_net,pov_index))
                  CASE (6)
                    catogrise_T_atom=621
                  CASE (5)
                    catogrise_T_atom=622
                  CASE (4)
                    catogrise_T_atom=623
                  CASE (3)
                    catogrise_T_atom=624
                  CASE (2)
                    catogrise_T_atom=625
                  CASE (1)
                    catogrise_T_atom=626
                  CASE (0)
                    catogrise_T_atom=627
                END SELECT
            END SELECT
          END IF
        END IF
      CASE (7)
        catogrise_T_atom=700
      CASE (8)
        catogrise_T_atom=800
      CASE DEFAULT
        catogrise_T_atom=801
    END SELECT
  END FUNCTION catogrise_T_atom
  FUNCTION T_X_bond_present(neighbour_list,neighbourNum)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8)  ::  neighbourNum, i
    LOGICAL ::  T_X_bond_present
    T_X_bond_present=.FALSE.
    DO i=1,neighbourNum
      IF (neighbour_list(1)%neighbour%element .NE. 'O') THEN
        T_X_bond_present=.TRUE.
      END IF
    END DO
  END FUNCTION T_X_bond_present
  FUNCTION return_max_min_angle(bond_angle_list,neighbourNum)
    REAL(8), ALLOCATABLE, INTENT(IN) ::  bond_angle_list(:,:)
    REAL(8) ::  return_max_min_angle(2)
    INTEGER(8), INTENT(IN)  ::  neighbourNum
    INTEGER(8)  ::  i,j
    return_max_min_angle=(/ 0.0, 360.0 /)
    DO i=1,neighbourNum
      DO j=1,neighbourNum
        IF (i==j) CYCLE
        IF (bond_angle_list(i,j) .GT. return_max_min_angle(1) ) &
& return_max_min_angle(1)=bond_angle_list(i,j)
        IF (bond_angle_list(i,j) .LT.  return_max_min_angle(2) ) &
&  return_max_min_angle(2)=bond_angle_list(i,j)
      END DO
    END DO
  END FUNCTION
  FUNCTION rightAngle_count(bond_angle_list,neighbourNum)
    REAL(8), ALLOCATABLE, INTENT(IN) ::  bond_angle_list(:,:)
    INTEGER(8), INTENT(IN)  ::  neighbourNum
    INTEGER(8)  ::  i,j, rightAngle_count
    rightAngle_count=0
    DO i=1,neighbourNum
      DO j=1,neighbourNum
        IF (i==j) CYCLE
        IF ( (bond_angle_list(i,j) > 80.0 ) .AND. ( bond_angle_list(i,j) < 110.0) ) & 
& rightAngle_count=rightAngle_count+1
      END DO
    END DO
  END FUNCTION
  INTEGER FUNCTION T_atom_count(neighbour_list,neighbourNum)
    TYPE(neighbour), INTENT(IN), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  neighbourNum
    INTEGER  ::  i
    T_atom_count=0
    DO i=1, neighbourNum
      IF ( ( neighbour_list(i)%neighbour%element == 'Al' ) .OR. &
& ( neighbour_list(i)%neighbour%element == 'Si' ) ) THEN
        T_atom_count=T_atom_count+1
      END IF
    END DO
  END FUNCTION T_atom_count
  INTEGER FUNCTION H_atom_count(neighbour_list,neighbourNum)
    TYPE(neighbour), INTENT(IN), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  neighbourNum
    INTEGER  ::  i
    H_atom_count=0
    DO i=1, neighbourNum
      IF ( neighbour_list(i)%neighbour%element == 'H' ) THEN
        H_atom_count=H_atom_count+1
      END IF
    END DO
  END FUNCTION H_atom_count
  FUNCTION silanol_count(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum, i
    INTEGER(8)  ::  neighboursIndx
    INTEGER(8)  ::  silanol_count, O_environment
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    silanol_count=0
    DO i=1,neighbourNum
      neighboursIndx=neighbour_list(i)%atom_num
      O_environment=connectivity_net(neighboursIndx)%local_environment
      IF (O_environment .EQ. 12 ) THEN
        silanol_count=silanol_count+1
      END IF
    END DO
  END FUNCTION  silanol_count
  FUNCTION monobonded_O_present(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum, i
    INTEGER(8)  ::  neighboursIndx, O_environment
    LOGICAL ::  monobonded_O_present
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    monobonded_O_present=.FALSE.
    DO i=1,neighbourNum 
      neighboursIndx=neighbour_list(i)%atom_num
      O_environment=connectivity_net(neighboursIndx)%local_environment
      IF (O_environment .EQ. 8 ) THEN
        monobonded_O_present=.TRUE.
      END IF
    END DO
  END FUNCTION monobonded_O_present
  FUNCTION count_framework_O_T(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum, i
    INTEGER(8)  ::  neighboursIndx
    INTEGER(8)  ::  count_framework_O_T, O_environment
! note this is a special O_frame count for the T_atoms which exclude silanols
! and only includes 10,11 or 14
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    count_framework_O_T=0
    DO i=1,neighbourNum
      neighboursIndx=neighbour_list(i)%atom_num
      O_environment=connectivity_net(neighboursIndx)%local_environment
      SELECT CASE(O_environment)
        CASE(10)
          count_framework_O_T=count_framework_O_T+1
        CASE(11)
          count_framework_O_T=count_framework_O_T+1
        CASE(14)
          count_framework_O_T=count_framework_O_T+1
        CASE DEFAULT
          count_framework_O_T=count_framework_O_T+0
      END SELECT
   !   IF ((O_environment .GE. 10) .AND. (O_environment .LE. 15 )) THEN
   !     count_framework_O=count_framework_O+1
   !   END IF
    END DO
  END FUNCTION count_framework_O_T
  FUNCTION count_framework_O(connectivity_net,pov_index)
    TYPE(connection_mesh), ALLOCATABLE, INTENT(INOUT)  ::  connectivity_net(:)
    TYPE(neighbour), ALLOCATABLE ::  neighbour_list(:)
    INTEGER(8), INTENT(IN)  ::  pov_index
    INTEGER(8)  ::  neighbourNum, i
    INTEGER(8)  ::  neighboursIndx
    INTEGER(8)  ::  count_framework_O, O_environment
    neighbourNum=connectivity_net(pov_index)%neighbourNum
    neighbour_list=connectivity_net(pov_index)%neighbour_list
    count_framework_O=0
    DO i=1,neighbourNum
      neighboursIndx=neighbour_list(i)%atom_num
      O_environment=connectivity_net(neighboursIndx)%local_environment
      IF ((O_environment .GE. 10) .AND. (O_environment .LE. 15 )) THEN
        count_framework_O=count_framework_O+1
      END IF
    END DO
  END FUNCTION count_framework_O
END MODULE neighbours
