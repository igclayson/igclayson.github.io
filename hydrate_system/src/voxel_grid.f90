MODULE voxel_grid
  USE types,  ONLY: atom, voxel 
  USE chemistry_methods, ONLY: electronegative, overlap
  IMPLICIT NONE
  PRIVATE
  PUBLIC  ::  generate_set_of_voxel
  PUBLIC  ::  hydrate_voxel
  CONTAINS
  FUNCTION generate_set_of_voxel(grid_finess,lat_a,lat_b,lat_c,&
&lat_alpha,lat_beta,lat_gamma)
    REAL(8), INTENT(IN) ::  lat_a,lat_b,lat_c
    REAL(8), INTENT(IN) ::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8), INTENT(IN)  ::  grid_finess ! this how much to divide my grid 
    TYPE(voxel), ALLOCATABLE ::  generate_set_of_voxel(:,:,:)
    INTEGER(2)  :: i, j, k
    REAL(8) ::  center_point(3)
    ALLOCATE(generate_set_of_voxel(grid_finess,grid_finess,grid_finess))
    DO i=1,grid_finess
      center_point(1)=i
      center_point(1)=(2*center_point(1)-1)/(2*grid_finess)
      DO j=1,grid_finess
        center_point(2)=j
        center_point(2)=(2*center_point(2)-1)/(2*grid_finess)
        DO k=1,grid_finess
          center_point(3)=k
          center_point(3)=(2*center_point(3)-1)/(2*grid_finess)
          generate_set_of_voxel(i,j,k)%center=fract_to_cart(center_point, &
        & lat_a,lat_b,lat_c,lat_alpha,lat_beta,lat_gamma)
        END DO
      END DO
    END DO
  END FUNCTION
  FUNCTION hydrate_voxel(voxel_grid,grid_finess,crystal,target_atom,atomNumber)
! this finds for all voxel centers the target separation, closest atom
! separation and the closet atom enumeration and the enumerate atoms image

    INTEGER(8), INTENT(IN)  ::  target_atom,atomNumber,grid_finess
    TYPE(atom), INTENT(IN)  ::  crystal(3,3,3,atomNumber)
    TYPE(voxel), INTENT(IN) ::  voxel_grid(grid_finess,grid_finess,grid_finess)
    INTEGER(4)  ::  i,j,k, vPoint(3),target_imageV(3),direction, closest_imageV(3)
    REAL(8) ::  r_hold, normV, alpha
    REAL(8) ::  target_vector(3), per_vector(3), centerV(3)
    LOGICAL ::  occupied
    TYPE(atom)  ::  hydrate_voxel(3) ! 1=O,2=H,3=H
    CHARACTER(len=2) :: element
    r_hold=1000
    DO i=1,grid_finess
      DO j=1,grid_finess
        DO k=1,grid_finess
          CALL examine_voxel(voxel_grid(i,j,k),crystal,target_atom,atomNumber)
          element=crystal(2,2,2,voxel_grid(i,j,k)%closest_atom)%element
          occupied=overlap(voxel_grid(i,j,k),element)
!          PRINT *, "occupied: ", occupied, "separation: ", &
!& voxel_grid(i,j,k)%target_separation
          IF ( (.NOT. occupied)&
& .AND. voxel_grid(i,j,k)%target_separation < r_hold) THEN
            r_hold=voxel_grid(i,j,k)%target_separation
            PRINT *, "r_hold: ", r_hold
            vPoint=(/i,j,k/)
            PRINT *, "vPoint: ", vPoint
            PRINT *, "voxel center: ", voxel_grid(i,j,k)%center
            PRINT *,""
          END IF
        END DO
      END DO
    END DO
    closest_imageV=voxel_grid(vPoint(1),vPoint(2),vPoint(3))%closest_image
    PRINT *, "closest image : ", closest_imageV
    PRINT *, "closeset atom ", voxel_grid(vPoint(1),vPoint(2),vPoint(3))%closest_atom
    target_imageV=voxel_grid(vPoint(1),vPoint(2),vPoint(3))%target_image
    element=crystal(closest_imageV(1)+2,closest_imageV(2)+2,closest_imageV(3)+2, &
    & voxel_grid(vPoint(1),vPoint(2),vPoint(3))%closest_atom)%element
    direction=electronegative(element)
    centerV=voxel_grid(vPoint(1),vPoint(2),vPoint(3))%center
    !PRINT *,"closest voxel: ",centerV
    normV=0
    !PRINT *, crystal(1,1,1,1)%cartesian(1)
!    PRINT *, &
!&  crystal(target_imageV(1),target_imageV(2),target_imageV(3),target_atom)
    DO i=1,3
      target_vector(i)=crystal(target_imageV(1)+2,target_imageV(2)+2,target_imageV(3)+2,target_atom)%cartesian(i)- &
      & centerV(i)
  !    PRINT *, target_vector(i)
      normV=normV+target_vector(i)**2
    END DO
    normV=SQRT(normV)
    DO i=1,3
      target_vector(i)=target_vector(i)/normV
    END DO
    alpha=(target_vector(1)+target_vector(2))/target_vector(3)
    alpha=SQRT(1/(1+1+alpha**2))
    per_vector(1)=alpha
    per_vector(2)=alpha
    per_vector(3)=-1*alpha*((target_vector(1)+target_vector(2))/target_vector(3))
    hydrate_voxel(1)%element="O"
    hydrate_voxel(2)%element="H"
    hydrate_voxel(3)%element="H"
    DO i=1,3
      hydrate_voxel(1)%cartesian(i)= centerV(i)
      hydrate_voxel(2)%cartesian(i)= centerV(i) -&
      & direction*target_vector(i)*0.4939 - direction*per_vector(i)*0.7591
      hydrate_voxel(3)%cartesian(i)= centerV(i) -&
      & direction*target_vector(i)*0.4939+ direction*per_vector(i)*0.7591
    END DO
   ! PRINT *,hydrate_voxel(1)%cartesian
   ! PRINT *,hydrate_voxel(2)%cartesian
   ! PRINT *,hydrate_voxel(3)%cartesian
  END FUNCTION hydrate_voxel
! private
  SUBROUTINE examine_voxel(Pvoxel,crystal,target_atom,atomNumber)
! note Pvoxel just means "this particular voxel"
    TYPE(voxel) ::  Pvoxel
    TYPE(atom), INTENT(IN)  ::  crystal(:,:,:,:)
    INTEGER(8), INTENT(IN)  ::  target_atom,atomNumber
    REAL(8) ::  r_min, r_target, r_hold
    INTEGER(8)  ::  i,j,k,a, n_min, imageC(3), imageT(3)
    r_target=1000
    r_min=1000
    DO i=1,3
      DO j=1,3
        DO k=1,3
          DO a=1,atomNumber
            r_hold=r(crystal(i,j,k,a)%cartesian,Pvoxel%center)
            IF (r_hold < r_min) THEN
                r_min=r_hold
                n_min=a
                imageC=(/i-2,j-2,k-2/)
            END IF
            IF (a==target_atom) THEN
              IF(r_hold < r_target) THEN
                r_target=r_hold
                imageT=(/i-2,j-2,k-2/)
              END IF
            END IF 
          END DO
        END DO
      END DO
    END DO
    Pvoxel%target_separation=r_target
    Pvoxel%target_image=imageT
    Pvoxel%closest_atom=n_min
    Pvoxel%closest_image=imageC
    Pvoxel%closest_separation=r_min
  END SUBROUTINE examine_voxel
  FUNCTION fract_to_cart(coords,lat_a,lat_b,lat_c, &
  & lat_alpha,lat_beta,lat_gamma) ! this is a lower memory version
    REAL(8) ::  fract_to_cart(3)
    REAL(8), INTENT(IN) ::  coords(3)
    REAL(8), INTENT(IN) ::  lat_a,lat_b,lat_c
    REAL(8), INTENT(IN) ::  lat_alpha,lat_beta,lat_gamma
    INTEGER(8)  ::  i
!   Note give in radians
    REAL(8) ::  const_x_1,const_x_2,const_y_1,const_y_2,const_z_1
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
    DO i=1, 3
      fract_to_cart(1)= coords(1)*lat_a &
      & + coords(2)*const_x_1 &
      & + coords(3)*const_x_2
      fract_to_cart(2) = &
      & coords(2)*const_y_1 &
      & + coords(3)*const_y_2
      fract_to_cart(3) = coords(3)*const_z_1
    END DO
  END FUNCTION
  REAL  FUNCTION r(arr1,arr2)
    REAL(8) :: arr1(3),arr2(3)
    INTEGER(2)  ::  i
    r=0
    DO i=1,3
      r=r+(arr1(i)-arr2(i))**2
    END DO
    r=SQRT(r)
  END FUNCTION 
END MODULE voxel_grid
