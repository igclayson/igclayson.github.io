MODULE crystal_module
  USE types,  ONLY: atom
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: read_from_file
  PUBLIC	:: generate_crystal
CONTAINS
  FUNCTION generate_crystal(atomNumber,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma)
    INTEGER(8), INTENT(IN)  ::  atomNumber
    TYPE(atom)  ::  coords(atomNumber)
    REAL(8), INTENT(IN)	::  lat_a,lat_b,lat_c
    REAL(8), INTENT(IN)	::  lat_alpha,lat_beta,lat_gamma
    TYPE(atom), DIMENSION(:,:,:,:), ALLOCATABLE  :: generate_crystal
    TYPE(atom)  ::  image_atom
    INTEGER(8)	::  i
    INTEGER(2)  ::  j,k,l,img1,img2,img3
    ! first wrap atoms to 000 cell
    CALL cart_to_fract(atomNumber,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma)
    DO i=1, atomNumber
     DO j=1, 3
      IF (coords(i)%fractional(j) > 1) THEN
        coords(i)%fractional(j)=coords(i)%fractional(j)-1
      END IF
      IF (coords(i)%fractional(j) < 0) THEN
        coords(i)%fractional(j) = coords(i)%fractional(j)+1
      END IF
     END DO
   END DO
   CALL fract_to_cart(atomNumber,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma)
  ALLOCATE(generate_crystal(3,3,3,atomNumber))
! remove overlapping of of ions in neighbouring cells
! do generate the 
    DO l=1,3
     img3=l-2
     DO k=1,3
     img2=k-2
      DO j=1,3
       img1=j-2
       DO i=1, atomNumber
       ! note that this is the (j k l) image for atom i with
! fractional coordinate x
        image_atom=coords(i)
        image_atom%fractional(1)=coords(i)%fractional(1)+1*img1
        image_atom%fractional(2)=coords(i)%fractional(2)+1*img2
        image_atom%fractional(3)=coords(i)%fractional(3)+1*img3
        generate_crystal(j,k,l,i)=image_atom
       END DO 
        CALL fract_to_cart(atomNumber,generate_crystal(j,k,l,1:atomNumber),lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma)
      END DO
     END DO
    END DO
  ! a function so must return generate_crystal
   ! coords have a fract and cartesian value
  END FUNCTION
  SUBROUTINE read_from_file(input,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma,atomNumber)
    CHARACTER(len=126), INTENT(IN)  ::  input
    REAL(8), INTENT(OUT) ::  lat_a,lat_b,lat_c
    REAL(8), INTENT(OUT) ::  lat_alpha,lat_beta,lat_gamma
    TYPE(atom), ALLOCATABLE, INTENT(OUT) ::  coords(:)
    INTEGER(8), INTENT(OUT) ::  atomNumber
    INTEGER(8)  ::  i
    REAL(8), PARAMETER  ::  pi=2*ASIN(1.0)
    ! read from file
    OPEN(unit=11, file=input)
    READ(11, *) atomNumber
    READ(11, *) lat_a,lat_b,lat_c,lat_alpha,lat_beta,lat_gamma
    lat_alpha=lat_alpha*pi/180
    lat_beta=lat_beta*pi/180
    lat_gamma=lat_gamma*pi/180
    ALLOCATE(coords(atomNumber))
    DO i=1, atomNumber
      READ(11, *) coords(i)%element, coords(i)%cartesian(1), &
      & coords(i)%cartesian(2), coords(i)%cartesian(3)
    END DO
    CLOSE(11,STATUS='KEEP')
  END SUBROUTINE read_from_file

! private
  SUBROUTINE cart_to_fract(atomNumber,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma) 
    TYPE(atom)  ::  coords(:)
    INTEGER(8), INTENT(IN)	::  atomNumber
    REAL(8), INTENT(IN)	::  lat_a,lat_b,lat_c
    REAL(8), INTENT(IN)	::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8)  ::  i
!   Note give in radians
    REAL(8)	::  const_u_1,const_u_2, const_v_1, const_v_2, const_w_1
    REAL(8) ::  const_vol_1,const_vol_2, vol
!   These constants are set to make the fractional coords
    const_vol_1=(1-COS(lat_alpha)**2-COS(lat_beta)**2-COS(lat_gamma)**2)
    const_vol_2=2*COS(lat_alpha)*COS(lat_beta)*COS(lat_gamma)
    vol=lat_a*lat_b*lat_c*SQRT(const_vol_1+const_vol_2)
    const_u_1= -COS(lat_gamma)/(lat_a*SIN(lat_gamma))
    const_u_2=lat_b*lat_c*(COS(lat_alpha)*COS(lat_gamma)-COS(lat_beta))/(vol*SIN(lat_gamma))
    const_v_1=1/(lat_b*SIN(lat_gamma))
    const_v_2=lat_a*lat_c*(COS(lat_beta)*COS(lat_gamma)-COS(lat_alpha))/(vol*SIN(lat_gamma))
    const_w_1=lat_a*lat_b*SIN(lat_gamma)/vol
    DO i=1, atomNumber
!     u
      coords(i)%fractional(1) = coords(i)%cartesian(1)/lat_a &
      & + coords(i)%cartesian(2)*const_u_1 &
      & + coords(i)%cartesian(3)*const_u_2
!add1+add2+add3
!     v
      coords(i)%fractional(2) = coords(i)%cartesian(2)*const_v_1 &
      & + coords(i)%cartesian(3)*const_v_2
!     w
      coords(i)%fractional(3) = coords(i)%cartesian(3)*const_w_1
    END DO
  END SUBROUTINE cart_to_fract

  SUBROUTINE fract_to_cart(atomNumber,coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma)
    TYPE(atom)  ::  coords(:)
    INTEGER(8), INTENT(IN)	::  atomNumber
    REAL(8), INTENT(IN)	::  lat_a,lat_b,lat_c
    REAL(8), INTENT(IN)	::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8)  ::  i
!   Note give in radians
    REAL(8)	::  const_x_1,const_x_2,const_y_1,const_y_2,const_z_1
    REAL(8) ::  const_vol_1,const_vol_2, vol
!   These constants are set to make the fractional coords
    const_vol_1=(1-COS(lat_alpha)**2-COS(lat_beta)**2-COS(lat_gamma)**2)
    const_vol_2=2*COS(lat_alpha)*COS(lat_beta)*COS(lat_gamma)
    vol=lat_a*lat_b*lat_c*SQRT(const_vol_1+const_vol_2)
    const_x_1=lat_b*COS(lat_gamma)
    const_x_2=lat_c*COS(lat_beta)
    const_y_1=lat_b*SIN(lat_gamma)
    const_y_2=lat_c*(COS(lat_alpha)-COS(lat_beta)*COS(lat_gamma))/SIN(lat_gamma)
    const_z_1=vol/(lat_a*lat_b*SIN(lat_gamma))
    DO i=1, atomNumber
      coords(i)%cartesian(1) = coords(i)%fractional(1)*lat_a &
      & + coords(i)%fractional(2)*const_x_1 &
      & + coords(i)%fractional(3)*const_x_2
      coords(i)%cartesian(2) = &
      & coords(i)%fractional(2)*const_y_1 &
      & + coords(i)%fractional(3)*const_y_2
      coords(i)%cartesian(3) = coords(i)%fractional(3)*const_z_1
    END DO
  END SUBROUTINE fract_to_cart
!  FUNCTION separation(
  
END MODULE crystal_module
