MODULE PDAFomi_dim_obs_l

  USE PDAFomi_obs_f, ONLY: obs_f, r_earth, pi, debug, n_obstypes, error
  USE PDAFomi_obs_l, ONLY: obs_l, obs_l_all, firstobs, offset_obs_l
  USE PDAF_mod_filtermpi, ONLY: mype, npes_filter

  USE omp_lib

  IMPLICIT NONE
  SAVE

  INTERFACE PDAFomi_init_dim_obs_l
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_iso
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso_locweights
  END INTERFACE

CONTAINS

  SUBROUTINE PDAFomi_init_dim_obs_l_iso(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius              !< Localization cut-off radius
    REAL, INTENT(in) :: sradius              !< Support radius of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of current observation vector

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_init_dim_obs_l_iso

  SUBROUTINE PDAFomi_check_dist2_loop(thisobs_l, thisobs, coordsA, cnt_obs, mode)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(in) :: thisobs       !< Data type with full observation
    REAL, INTENT(in) :: coordsA(:)           !< Coordinates of current analysis domain (ncoord)
    INTEGER, INTENT(inout) :: cnt_obs        !< Count number of local observations
    INTEGER, INTENT(in) :: mode              !< 1: count local observations
                                             !< 2: initialize local arrays
    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_check_dist2_loop

  SUBROUTINE PDAFomi_init_dim_obs_l_noniso(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius(:)           !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)           !< Vector of support radii of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of current observation vector

! *** Local variables ***
    REAL :: maxcoords_l, mincoords_l         ! Min/Max domain coordinates to check geographic coords
    REAL :: maxocoords_l, minocoords_l       ! Min/Max observation coordinates to check geographic coords
    INTEGER :: cnt_obs                       ! Counter for valid local observations

    ! WRITE(0,*) "RSE: ENTER PDAFomi_init_dim_obs_l_noniso"

    doassim: IF (thisobs%doassim == 1) THEN

! ***********************************************
! *** Check offset in full observation vector ***
! ***********************************************
       
       ! Check consistency of dimensions
       IF (SIZE(cradius) /= thisobs%ncoord) THEN          
          error = 12
       END IF
       IF (SIZE(sradius) /= thisobs%ncoord) THEN          
          error = 13
       END IF
       IF (thisobs%ncoord/=3 .AND. thisobs%disttype>=10) THEN          
          error = 14
       END IF


! **************************************
! *** Store localization information ***
! **************************************

       thisobs_l%locweight = locweight

       ! Allocate vectors for localization radii and store their values
       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(thisobs%ncoord))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(thisobs%ncoord))

       thisobs_l%nradii = thisobs%ncoord
       thisobs_l%cradius(:) = cradius(:)
       thisobs_l%sradius(:) = sradius(:)


! **************************************
! *** Count valid local observations ***
! **************************************       

       cnt_obs = 0
       IF (thisobs_l%nradii==1) THEN
          ! 1D but with radius specified as array
          STOP "RSE: INTERRUPT #1"
          CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 1)
       ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN
          ! Nonisotropic in 2 or 3 dimensions
          CALL PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coords_l, cnt_obs, 1)
       ELSE          
          error = 10
          STOP "RSE: INTERRUPT #2"
       END IF


! ************************************************
! *** Initialize local observation for PDAFomi ***
! ************************************************

       CALL PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs)


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************       

       ! Count local observations and initialize index and distance arrays
       IF (thisobs_l%dim_obs_l>0) THEN

          cnt_obs = 0
          IF (thisobs_l%nradii==1) THEN
             ! 1D but with radius specified as array
             STOP "RSE: INTERRUPT #3"
             CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 2)
          ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN
             ! Nonisotropic in 2 or 3 dimensions
             CALL PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coords_l, cnt_obs, 2)
          ELSE             
             error = 11
             STOP "RSE: INTERRUPT #4"
          END IF
       END IF              

    END IF doassim

    ! WRITE(0,*) "RSE: EXIT  PDAFomi_init_dim_obs_l_noniso"

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso

  SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights(thisobs_l, thisobs, coords_l, locweights, cradius, &
       sradius, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweights(:)     !< Types of localization function
    REAL, INTENT(in) :: cradius(:)           !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)           !< Vector of support radii of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of current observation vector

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights

  SUBROUTINE PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coordsA, cnt_obs, mode)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coordsA(:)        !< Coordinates of current analysis domain (ncoord)
    INTEGER, INTENT(inout) :: cnt_obs     !< Count number of local observations
    INTEGER, INTENT(in) :: mode              !< 1: count local observations
                                             !< 2: initialize local arrays

! *** Local variables ***
    INTEGER :: i, k                 ! Counters
    ! INTEGER :: verbose              ! verbosity flag
    INTEGER :: domsize              ! Flag whether domainsize is set
    LOGICAL :: distflag             ! Flag whether distance in a coordinate direction is within cradius
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    REAL :: distance2               ! square distance
    REAL :: cradius2                ! cut-off radius on ellipse or ellipsoid
    REAL :: phi, theta              ! Angles in ellipse or ellipsoid
    REAL :: dist_xy                 ! Distance in xy-plan in 3D case

    ! REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    ! REAL :: coordsB(thisobs%ncoord) ! Array for coordinates of a single observation
    REAL :: dists(16)   ! Distance vector between analysis point and observation
    REAL :: coordsB(16) ! Array for coordinates of a single observation

    REAL :: cradius                 ! Directional cut-off radius
    REAL :: sradius                 ! Directional support radius
    LOGICAL :: checkdist            ! Flag whether distance is within cut-off radius

! **********************
! *** RSE CODE BLOCK ***
! **********************
    INTEGER :: local_cnt_obs
    INTEGER :: global_cnt_obs
    INTEGER :: local_mode
    LOGICAL :: allocFlags(16)    
    INTEGER :: accStreamId
    
    domsize = 1
    local_mode = mode
    global_cnt_obs = cnt_obs    
    allocFlags = .FALSE.

    accStreamId = omp_get_thread_num() + 1    

! **********************
! *** Initialization ***
! **********************

    ! WRITE(0,'(1x, a, 1x, I6, 1x, I6, 1x, I6)') "RSE: ENTER PDAFomi_check_dist2_noniso_loop", &
    !                                             thisobs%dim_obs_f, mode, cnt_obs
    
    IF( (local_mode == 2) .AND. (.NOT.ALLOCATED(thisobs_l%id_obs_l)) ) THEN
        WRITE(0,*) "RSE: WARNING ##1"
    END IF

#if 0
    !$ACC ENTER DATA CREATE(thisobs, thisobs_l, global_cnt_obs, coordsA(:)) ASYNC(accStreamId)
    !$ACC UPDATE     DEVICE(thisobs, thisobs_l, global_cnt_obs, coordsA(:)) ASYNC(accStreamId)
    !$ACC ENTER DATA CREATE(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) ASYNC(accStreamId)
    !$ACC UPDATE     DEVICE(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) ASYNC(accStreamId)
#else
    !$ACC ENTER DATA CREATE(thisobs_l, global_cnt_obs, coordsA(:)) ASYNC(accStreamId)
    !$ACC UPDATE     DEVICE(thisobs_l, global_cnt_obs, coordsA(:)) ASYNC(accStreamId)
    !$ACC ENTER DATA CREATE(thisobs_l%cradius(:), thisobs_l%sradius(:)) ASYNC(accStreamId)
    !$ACC UPDATE     DEVICE(thisobs_l%cradius(:), thisobs_l%sradius(:)) ASYNC(accStreamId)
#endif
    IF(ALLOCATED(thisobs_l%id_obs_l)) THEN
        !$ACC ENTER DATA CREATE(thisobs_l%id_obs_l(:), thisobs_l%distance_l(:), thisobs_l%cradius_l(:)) &
        !$ACC            CREATE(thisobs_l%sradius_l(:), thisobs_l%dist_l_v(:)) ASYNC(accStreamId)
        !$ACC UPDATE     DEVICE(thisobs_l%id_obs_l(:), thisobs_l%distance_l(:), thisobs_l%cradius_l(:)) &
        !$ACC            DEVICE(thisobs_l%sradius_l(:), thisobs_l%dist_l_v(:)) ASYNC(accStreamId)
    ENDIF

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
    !$ACC PRESENT(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) &
    !$ACC PRESENT(thisobs_l%id_obs_l(:), thisobs_l%distance_l(:), thisobs_l%cradius_l(:)) &
    !$ACC PRESENT(thisobs_l%sradius_l(:), thisobs_l%dist_l_v(:)) &
    !$ACC PRESENT(thisobs_l, thisobs, global_cnt_obs, coordsA(:)) &
    !$ACC PRIVATE(k, distflag, slon, slat, distance2, cradius2, phi, theta, dist_xy, dists(:)) &
    !$ACC PRIVATE(coordsB(:), cradius, sradius, checkdist, local_cnt_obs) &
    !$ACC FIRSTPRIVATE(local_mode,  domsize) &
    !$ACC ASYNC(accStreamId)
    scancount: DO i = 1, thisobs%dim_obs_f       

       ! Initialize distance flag
       checkdist = .FALSE.    ! Whether an observation lies within the local box
       distflag = .TRUE.      ! Whether an observation lies within the local radius (ellipse, ellipsoid)       

       ! Observation coordinates
       coordsB = thisobs%ocoord_f(1:thisobs%ncoord, i)


! ************************
! *** Compute distance ***
! ************************

       norm: IF ((thisobs%disttype==0 .OR. thisobs%disttype==10) .OR. &
            ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==0)) THEN
          WRITE(0,*) "DEBUG RSE: ##0"
          ! *** Compute Cartesian distance ***          

       ELSEIF ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==1) THEN norm

          ! WRITE(0,*) "DEBUG RSE: ###5"

          ! *** Compute periodic Cartesian distance ***          

          IF (thisobs%ncoord==3) THEN
             ! WRITE(0,*) "DEBUG RSE: ###6"
             
             IF (thisobs%domainsize(3)<=0.0) THEN 
                dists(3) = ABS(coordsA(3) - coordsB(3))
             ELSE
                dists(3) = MIN(ABS(coordsA(3) - coordsB(3)), &
                     ABS(ABS(coordsA(3) - coordsB(3))-thisobs%domainsize(3)))
             END IF

             IF (dists(3)>thisobs_l%cradius(3)) THEN
                distflag = .FALSE.
             ELSE
                IF (thisobs%domainsize(2)<=0.0) THEN 
                   dists(2) = ABS(coordsA(2) - coordsB(2))
                ELSE
                   dists(2) = MIN(ABS(coordsA(2) - coordsB(2)), &
                        ABS(ABS(coordsA(2) - coordsB(2))-thisobs%domainsize(2)))
                END IF
                
                IF (dists(2)>thisobs_l%cradius(2)) THEN
                   distflag = .FALSE.
                ELSE
                   IF (thisobs%domainsize(1)<=0.0) THEN 
                      dists(1) = ABS(coordsA(1) - coordsB(1))
                   ELSE
                      dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                           ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
                   END IF
                   IF (dists(1)>thisobs_l%cradius(1)) THEN
                      distflag = .FALSE.
                   ELSE
                      ! full squared distance
                      distance2 = 0.0
                      IF (thisobs%disttype<10) THEN
                         ! full 3D localization
                         ! !$ACC LOOP SEQ
                         DO k = 1, thisobs%ncoord
                            distance2 = distance2 + dists(k)*dists(k)
                         END DO
                         ! !$ACC END LOOP
                      ELSE
                         ! factorized 2+1D localization
                         !$ACC LOOP SEQ
                         DO k = 1, thisobs%ncoord-1
                            distance2 = distance2 + dists(k)*dists(k)
                         END DO
                         !$ACC END LOOP
                      END IF
                   END IF
                END IF
             END IF
          ELSEIF (thisobs%ncoord==2) THEN
             WRITE(0,*) "DEBUG RSE: ###7"
             
          ELSEIF (thisobs%ncoord==1) THEN
             WRITE(0,*) "DEBUG RSE: ###8"
             
          END IF

       ELSEIF (thisobs%disttype==2 .OR. thisobs%disttype==12) THEN norm
          WRITE(0,*) "DEBUG RSE: ###9"          

       ELSEIF (thisobs%disttype==3 .OR. thisobs%disttype==13) THEN norm
          WRITE(0,*) "DEBUG RSE: ###12"

       END IF norm


! ***************************************************************************
! *** Compute directional cut-off and support radii and set distance flag ***
! ***************************************************************************       

       dflag: IF (distflag) THEN        

          nrad: IF (thisobs_l%nradii == 2 .OR. (thisobs_l%nradii == 3 .AND. thisobs%disttype >= 10)) THEN
             
             ! WRITE(0,*) "DEBUG RSE: ###20"

             IF ((thisobs_l%cradius(1) == thisobs_l%cradius(2)) .OR. &
                  (thisobs_l%sradius(1) == thisobs_l%sradius(2))) THEN
                ! 2D isotropic case

                ! WRITE(0,*) "DEBUG RSE: ###30"

                cradius2 = thisobs_l%cradius(1) * thisobs_l%cradius(1)

                IF (distance2 <= cradius2) THEN
                   
                   ! Set flag for valid observation
                   checkdist = .TRUE.
                   ! checkCounterRed = checkCounterRed + 1

                   !$ACC ATOMIC CAPTURE
                   global_cnt_obs = global_cnt_obs + 1
                   local_cnt_obs = global_cnt_obs
                   !$ACC END ATOMIC

                   cradius = thisobs_l%cradius(1)
                   sradius = thisobs_l%sradius(1)

                END IF
             ELSE

                WRITE(0,*) "DEBUG RSE: ###31"                

             END IF

          ELSE IF (thisobs_l%nradii == 3  .AND. thisobs%disttype < 10) THEN nrad
             
             WRITE(0,*) "DEBUG RSE: ###21"

          ELSEIF (thisobs_l%nradii == 1) THEN nrad

             WRITE(0,*) "DEBUG RSE: ###22"

          END IF nrad

          IF (local_mode==2 .AND. checkdist) THEN

             ! For internal storage (use in prodRinvA_l and likelihood_l)

             thisobs_l%id_obs_l(local_cnt_obs) = i                       ! node index
             thisobs_l%distance_l(local_cnt_obs) = SQRT(distance2)       ! distance
             thisobs_l%cradius_l(local_cnt_obs) = cradius                ! directional cut-off radius
             thisobs_l%sradius_l(local_cnt_obs) = sradius                ! directional support radius
             IF (thisobs_l%locweight_v>0 .AND. thisobs_l%nradii==3) THEN                
                thisobs_l%dist_l_v(local_cnt_obs) = dists(3)             ! distance in vertical direction
             END if
          END IF
    END IF dflag

 END DO scancount
 ! !$ACC END PARALLEL LOOP

 IF(ALLOCATED(thisobs_l%id_obs_l)) THEN
    !$ACC UPDATE      HOST(thisobs_l%id_obs_l(:), thisobs_l%distance_l(:), thisobs_l%cradius_l(:)) &
    !$ACC             HOST(thisobs_l%sradius_l(:), thisobs_l%dist_l_v(:)) ASYNC(accStreamId)
    !$ACC EXIT DATA DELETE(thisobs_l%id_obs_l(:), thisobs_l%distance_l(:), thisobs_l%cradius_l(:)) &
    !$ACC           DELETE(thisobs_l%sradius_l(:), thisobs_l%dist_l_v(:)) ASYNC(accStreamId)
 ENDIF

#if 0
 !$ACC UPDATE      HOST(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) ASYNC(accStreamId)
 !$ACC EXIT DATA DELETE(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) ASYNC(accStreamId)
 !$ACC UPDATE      HOST(thisobs, thisobs_l, global_cnt_obs, coordsA(:)) ASYNC(accStreamId)
 !$ACC EXIT DATA DELETE(thisobs, thisobs_l, global_cnt_obs, coordsA) ASYNC(accStreamId)
#else
#if 0
 !$ACC EXIT DATA DELETE(thisobs_l%cradius(:), thisobs_l%sradius(:), thisobs%domainsize(:), thisobs%ocoord_f(:,:)) ASYNC(accStreamId)
 !$ACC UPDATE      HOST(thisobs_l, global_cnt_obs) ASYNC(accStreamId)
 !$ACC EXIT DATA DELETE(thisobs, thisobs_l, global_cnt_obs, coordsA) ASYNC(accStreamId)
#else
 !$ACC EXIT DATA DELETE(thisobs_l%cradius(:), thisobs_l%sradius(:)) ASYNC(accStreamId)
 !$ACC UPDATE      HOST(thisobs_l, global_cnt_obs) ASYNC(accStreamId)
 !$ACC EXIT DATA DELETE(thisobs_l, global_cnt_obs, coordsA) ASYNC(accStreamId)
#endif
#endif
 !$ACC WAIT(accStreamId)

 cnt_obs = global_cnt_obs  
 
#if 0
 IF(ALLOCATED(thisobs_l%id_obs_l  )) allocFlags( 1) = .TRUE.
 IF(ALLOCATED(thisobs_l%distance_l)) allocFlags( 2) = .TRUE.
 IF(ALLOCATED(thisobs_l%cradius_l )) allocFlags( 3) = .TRUE.
 IF(ALLOCATED(thisobs_l%sradius_l )) allocFlags( 4) = .TRUE.
 IF(ALLOCATED(thisobs_l%dist_l_v  )) allocFlags( 5) = .TRUE.
 IF(ALLOCATED(thisobs_l%cradius   )) allocFlags( 6) = .TRUE.
 IF(ALLOCATED(thisobs_l%sradius   )) allocFlags( 7) = .TRUE.
 IF(ALLOCATED(thisobs%domainsize  )) allocFlags( 8) = .TRUE.
 IF(ALLOCATED(thisobs%ocoord_f    )) allocFlags( 9) = .TRUE.
 ! IF(ALLOCATED(coordsA           )) allocFlags(10) = .TRUE. 
 WRITE(0,*) "ALLOCATE: ", allocFlags(1:9)
 IF(ALLOCATED(thisobs_l%id_obs_l  )) WRITE(0,*) "thisobs_l%id_obs_l   ", SIZE(thisobs_l%id_obs_l  )
 IF(ALLOCATED(thisobs_l%distance_l)) WRITE(0,*) "thisobs_l%distance_l ", SIZE(thisobs_l%distance_l)
 IF(ALLOCATED(thisobs_l%cradius_l )) WRITE(0,*) "thisobs_l%cradius_l  ", SIZE(thisobs_l%cradius_l )
 IF(ALLOCATED(thisobs_l%sradius_l )) WRITE(0,*) "thisobs_l%sradius_l  ", SIZE(thisobs_l%sradius_l )
 IF(ALLOCATED(thisobs_l%dist_l_v  )) WRITE(0,*) "thisobs_l%dist_l_v   ", SIZE(thisobs_l%dist_l_v  )
 IF(ALLOCATED(thisobs_l%cradius   )) WRITE(0,*) "thisobs_l%cradius    ", SIZE(thisobs_l%cradius   )
 IF(ALLOCATED(thisobs_l%sradius   )) WRITE(0,*) "thisobs_l%sradius    ", SIZE(thisobs_l%sradius   )
 IF(ALLOCATED(thisobs%domainsize  )) WRITE(0,*) "thisobs%domainsize   ", SIZE(thisobs%domainsize  )
 IF(ALLOCATED(thisobs%ocoord_f    )) WRITE(0,*) "thisobs%ocoord_f     ", SIZE(thisobs%ocoord_f,  1), SIZE(thisobs%ocoord_f,   2)
 IF(.TRUE.                         ) WRITE(0,*) "coordsA              ", SIZE(coordsA             ) 
 WRITE(0,'(1x, a, 1x, I3, 1x, I6)') "RSE: PDAFomi_check_dist2_noniso_loop ", mode, cnt_obs
#endif 

  END SUBROUTINE PDAFomi_check_dist2_noniso_loop

  SUBROUTINE PDAFomi_set_localization(thisobs_l, cradius, sradius, locweight)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: cradius         !< Localization cut-off radius
    REAL, INTENT(in) :: sradius         !< Support radius of localization function
    INTEGER, INTENT(in) :: locweight    !< Type of localization function

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_set_localization

  SUBROUTINE PDAFomi_set_localization_noniso(thisobs_l, nradii, cradius, sradius, locweight, locweight_v)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: nradii            !< Number of radii to consider for localization
    REAL, INTENT(in) :: cradius(nradii)      !< Localization cut-off radius
    REAL, INTENT(in) :: sradius(nradii)      !< Support radius of localization function
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    INTEGER, INTENT(in) :: locweight_v       !< Type of localization function in vertical direction (only for nradii=3)

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_set_localization_noniso

  SUBROUTINE PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of observation vector over all obs. types
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of single observation type vector

    ! WRITE(0,'(1x, a, 1x, I6, 1x, I6)') "RSE: ENTER PDAFomi_set_dim_obs_l", cnt_obs_l_all, cnt_obs_l

    ! Store ID of first observation type that calls the routine
    ! This is reset in PDAFomi_deallocate_obs
    IF (firstobs == 0) THEN
       firstobs = thisobs%obsid
    END IF
 
    ! Reset offset of currrent observation in overall local obs. vector
    IF (thisobs%obsid == firstobs) THEN
       offset_obs_l = 0
       cnt_obs_l_all = 0
    END IF

    ! Store offset
    thisobs_l%off_obs_l = offset_obs_l

    ! Initialize pointer array
    IF (thisobs%obsid == firstobs) THEN
       IF (ALLOCATED(obs_l_all)) DEALLOCATE(obs_l_all)
       ALLOCATE(obs_l_all(n_obstypes))
    END IF

    ! Set pointer to current observation
    obs_l_all(thisobs%obsid)%ptr => thisobs_l

    ! Store local observation dimension and increment offset
    thisobs_l%dim_obs_l = cnt_obs_l
    offset_obs_l = offset_obs_l + cnt_obs_l
    cnt_obs_l_all = cnt_obs_l_all + cnt_obs_l

    ! Allocate arrays to store information on local observations
    IF (ALLOCATED(thisobs_l%id_obs_l)) DEALLOCATE(thisobs_l%id_obs_l)
    IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
    IF (ALLOCATED(thisobs_l%cradius_l)) DEALLOCATE(thisobs_l%cradius_l)
    IF (ALLOCATED(thisobs_l%sradius_l)) DEALLOCATE(thisobs_l%sradius_l)

    haveobs: IF (cnt_obs_l>0) THEN
       ALLOCATE(thisobs_l%id_obs_l(cnt_obs_l))
       ALLOCATE(thisobs_l%distance_l(cnt_obs_l))
       ALLOCATE(thisobs_l%cradius_l(cnt_obs_l))
       ALLOCATE(thisobs_l%sradius_l(cnt_obs_l))
       IF (thisobs_l%locweight_v>0) THEN
          IF (ALLOCATED(thisobs_l%dist_l_v)) DEALLOCATE(thisobs_l%dist_l_v)
          ALLOCATE(thisobs_l%dist_l_v(cnt_obs_l))
       END IF

    ELSE
       ALLOCATE(thisobs_l%id_obs_l(1))
       ALLOCATE(thisobs_l%distance_l(1))
       ALLOCATE(thisobs_l%cradius_l(1))
       ALLOCATE(thisobs_l%sradius_l(1))
       IF (ALLOCATED(thisobs_l%dist_l_v)) DEALLOCATE(thisobs_l%dist_l_v)
       ALLOCATE(thisobs_l%dist_l_v(1))
    END IF haveobs

    ! WRITE(0,'(1x, a, 1x, I6, 1x, I6)') "RSE: EXIT  PDAFomi_set_dim_obs_l", cnt_obs_l_all, cnt_obs_l

  END SUBROUTINE PDAFomi_set_dim_obs_l

  SUBROUTINE PDAFomi_store_obs_l_index(thisobs_l, idx, id_obs_l, distance, &
       cradius_l, sradius_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: idx       !< Element of local observation array to be filled
    INTEGER, INTENT(in) :: id_obs_l  !< Index of local observation in full observation array
    REAL, INTENT(in) :: distance     !< Distance between local analysis domain and observation
    REAL, INTENT(in) :: cradius_l    !< cut-off radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: sradius_l    !< support radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_store_obs_l_index

  SUBROUTINE PDAFomi_store_obs_l_index_vdist(thisobs_l, idx, id_obs_l, distance, &
       cradius_l, sradius_l, vdist)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: idx       !< Element of local observation array to be filled
    INTEGER, INTENT(in) :: id_obs_l  !< Index of local observation in full observation array
    REAL, INTENT(in) :: distance     !< Distance between local analysis domain and observation
    REAL, INTENT(in) :: cradius_l    !< cut-off radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: sradius_l    !< support radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: vdist        !< support radius in vertical direction for 2+1D factorized localization

    STOP "RSE: INTERRRUPT"

  END SUBROUTINE PDAFomi_store_obs_l_index_vdist

END MODULE PDAFomi_dim_obs_l
