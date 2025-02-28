!> PDAF-OMI observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! The subroutines in this module are for the particular handling of
!! a single observation type:
!!     EN4 salinity
!! The routines are called by the different call-back routines of PDAF
!! usually by callback_obs_pdafomi.F90
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data types (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! **Using this template:**
!!   To be able to distinguish the observation type and the routines in this module,
!!   we recommend to rename the module according to the observation module-type.
!!   Further,we recommend to replace 'OBSTYPE' in the routine names according to the
!!   type of the observation so that they can be identified when calling them from 
!!   the call-back routines.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_OBSTYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_OBSTYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routines, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_OBSTYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_OBSTYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2024-03 - Lars Nerger - Initial code for ICON
!! * Later revisions - see repository log
!!
MODULE obs_EN4_sao_pdafomi

  USE mod_parallel, &
       ONLY: mype_filter, npes_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_EN4_sao = .FALSE. !< Whether to assimilate this data type
  REAL    :: rms_obs_EN4_sao         !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


! *********************************************************
! *** Data type obs_f defines the full observations by  ***
! *** internally shared variables of the module         ***
! *********************************************************

! Relevant variables that can be modified by the user:
!   TYPE obs_f
!      ---- Mandatory variables to be set in INIT_DIM_OBS ----
!      INTEGER :: doassim                    ! Whether to assimilate this observation type
!      INTEGER :: disttype                   ! Type of distance computation to use for localization
!                                            ! (0) Cartesian, (1) Cartesian periodic
!                                            ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                     ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!
!      ---- Optional variables - they can be set in INIT_DIM_OBS ----
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!      ---- Variables with predefined values - they can be changed in INIT_DIM_OBS  ----
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!      REAL :: inno_omit=0.0                ! Omit obs. if squared innovation larger this factor times
!                                           !     observation variance
!      REAL :: inno_omit_ivar=1.0e-12       ! Value of inverse variance to omit observation
!   END TYPE obs_f

! Data type obs_l defines the local observations by internally shared variables of the module

! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could rename the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - array with indices of module-type observation in process-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!! * thisobs\%inno_omit   - Omit obs. if squared innovation larger this factor times observation variance
!!                          (default: 0.0, omission is active if >0) 
!! * thisobs\%inno_omit_ivar - Value of inverse variance to omit observation
!!                          (default: 1.0e-12, change this if this value is not small compared to actual obs. error)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE init_dim_obs_EN4_sao(step, dim_obs)

    USE mod_memcount, &     ! Counting allocated memory
         ONLY: memcount
    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assimilation, &
         ONLY: screen, cradius, rmsfileflag, rms_min, debug, &
         disttype, scale_vert
    USE mod_model, &
         ONLY: dim_xyz_p, dimxy_p, dimz, &
         idx_wet_in_all, dims_2d_p, offs_2d_p, pi, rad2deg
    USE mod_statevector_pdaf, &
         ONLY: sfields, id
    USE parser, &
         ONLY: handle, parse
    USE mod_io_pdaf, &
         ONLY: replace_slash, nf_check

    IMPLICIT NONE

    INCLUDE 'netcdf.inc' 

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: ixy, iz, rec              ! Counters
    INTEGER :: irms, obsrec              ! Counters
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    REAL, ALLOCATABLE :: obs_xyz(:,:)    ! Observations read from file
    REAL, ALLOCATABLE :: rms_xyz(:,:)    ! Observation errors read from file or initialize directly
    REAL, ALLOCATABLE :: lon_1d(:), lat_1d(:)  ! Longitude and latitude of observations

    INTEGER :: dimt
    INTEGER :: countv(3), startv(3)
    INTEGER :: ncid_in, id_dim, id_state, id_rms, id_lon, id_lat
    CHARACTER(len=120) :: obspath
    CHARACTER(len=120) :: obsdir
    CHARACTER(len=120) :: obsfile
    CHARACTER(len=150) :: ncfile_in


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - EN4_sao'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_EN4_sao) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = disttype   ! 11: 2+1D factorized Cartesian periodic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 3

    ! Initialize flag for type of full observations
    IF (npes_filter>1) thisobs%use_global_obs = 0

    IF (disttype==1 .OR. disttype==11) THEN
       ALLOCATE(thisobs%domainsize(thisobs%ncoord))
       thisobs%domainsize(1) = 360.0
       thisobs%domainsize(2) = 0.0
       thisobs%domainsize(3) = 0.0
    END IF


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Path to and name of file holding observations (DEFAULT Setting)
    obspath = ''
    obsdir = ''
    obsfile = ''

    ! Parse obs data input dir and file name 
    handle = 'obsdir' 
    CALL parse(handle, obsdir)
    handle = 'obsfile' 
    CALL parse(handle, obsfile)

    CALL replace_slash(obsdir, obspath)

    ncfile_in = TRIM(obspath)//'/'//TRIM(obsfile)
    IF (mype_filter==0) &
         WRITE (*,*) 'Read observations from file: ',TRIM(ncfile_in)


! *************************
! *** Read observations ***
! *************************

    CALL nf_check(NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in))

    CALL nf_check(NF_INQ_DIMID(ncid_in, 'time', id_dim))
    CALL nf_check(NF_INQ_DIMLEN(ncid_in, id_dim, dimt))
    CALL nf_check(NF_INQ_VARID(ncid_in, sfields(id%sao)%cvar, id_state))
    IF (rmsfileflag == 1) THEN
       CALL nf_check(NF_INQ_VARID(ncid_in, sfields(id%sao)%crms, id_rms))
    END IF
    CALL nf_check(NF_INQ_VARID(ncid_in, 'clon', id_lon))
    CALL nf_check(NF_INQ_VARID(ncid_in, 'clat', id_lat))


    ! *** Read observation values ***

    IF (debug) WRITE (*,'(/1x,a)') '------- Read observations -------------'

    ALLOCATE(obs_xyz(dimxy_p,dimz))
    CALL memcount(4, 'd', dimxy_p*dimz)

    ! Read observation array
    startv(3) = 1
    countv(3) = dimt 
    startv(2) = 1
    countv(2) = dimz
    startv(1) = offs_2d_p(mype_filter+1)+1
    countv(1) = dims_2d_p(mype_filter+1)

    CALL nf_check(NF_GET_VARA_DOUBLE(ncid_in, id_state, startv, countv, obs_xyz))

    !monitor
    IF (debug) WRITE(*,*) 'Read complete for '//TRIM(sfields(id%sao)%cvar)// &
         ': min(xyz):', MINVAL(obs_xyz) &
         ,', max(xyz):',MAXVAL(obs_xyz) 


    ! *** Read coordinates ***

    ALLOCATE(lon_1d(dimxy_p))
    ALLOCATE(lat_1d(dimxy_p))
    CALL memcount(4, 'd', 2*dimxy_p)

    startv(1) = offs_2d_p(mype_filter+1)+1
    countv(1) = dims_2d_p(mype_filter+1)

    CALL nf_check(NF_GET_VARA_DOUBLE(ncid_in, id_lon, startv(1), countv(1), lon_1d))
    CALL nf_check(NF_GET_VARA_DOUBLE(ncid_in, id_lat, startv(1), countv(1), lat_1d))

    ! shift lontitude to range 0 to 360 degrees
    IF (disttype==0 .OR. disttype==1 .OR. disttype==10 .OR. disttype==11) THEN
       lon_1d = lon_1d * rad2deg + 180.0
       lat_1d = lat_1d * rad2deg
    ELSE
       lon_1d = lon_1d + pi
       lat_1d = lat_1d
    END IF


    ! *** Read observation errors ***

    ALLOCATE(rms_xyz(dimxy_p,dimz))
    CALL memcount(4, 'd', dimxy_p*dimz)

    IF (rmsfileflag == 1) THEN                    ! read from file if flag is set

       IF (mype_filter==0) &
            WRITE(*,*) 'Read rms_obs from file.'

       ! Read state array
       startv(3) = 1
       countv(3) = dimt 
       startv(2) = 1
       countv(2) = dimz
       startv(1) = offs_2d_p(mype_filter+1)+1
       countv(1) = dims_2d_p(mype_filter+1)

       CALL nf_check(NF_GET_VARA_DOUBLE(ncid_in, id_rms, startv, countv, rms_xyz))

       !monitor
       WRITE(*,*) 'Read complete for '//TRIM(sfields(id%sao)%crms)// &
            ': min(xyz):', MINVAL(rms_xyz) &
            ,', max(xyz):',MAXVAL(rms_xyz)

    ELSE        

       ! *** take rms_obs from command line if flag is not set ***
       IF (mype_filter==0) &
            WRITE(*,'(8x,a,es10.3)') 'rms_obs is fixed with rms_obs_EN4_sao:', rms_obs_EN4_sao

       rms_xyz(:,:) = rms_obs_EN4_sao

    END IF

    ! *** Close input files
    CALL nf_check(nf_close(ncid_in))


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    IF (mype_filter==0) &
         WRITE(*,'(8x,a)') 'Scanning for valid observations'

    obsrec = 0
    DO iz = 1, dimz
       DO ixy = 1, dimxy_p
          IF (obs_xyz(ixy,iz) > -1e33) THEN
             obsrec = obsrec + 1
          END IF
       END DO
    END DO

    dim_obs_p = obsrec

    IF (npes_filter==1) THEN
       WRITE (*,'(8x, a, i7)') '--- number of observations from EN4_sao: ', dim_obs_p
    ELSE
       IF (screen>2) THEN
          WRITE (*,'(8x, a, i4, 2x, a, i7)') 'PE', mype_filter, &
               '--- number of observations from EN4_sao: ', dim_obs_p
       END IF
    END IF
    

    ! *** Initialize vector of observations on the process sub-domain ***

    IF (dim_obs_p>0) THEN

       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(ocoord_p(thisobs%ncoord, dim_obs_p))
       ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))
       ALLOCATE(ivar_obs_p(dim_obs_p))
       CALL memcount(4, 'd', 2*dim_obs_p + thisobs%ncoord*dim_obs_p)
       CALL memcount(4, 'i', dim_obs_p)

       irms = 0
       obsrec = 0
       rec = 0
       DO iz = 1, dimz
          DO ixy = 1, dimxy_p

             rec = rec + 1

             IF (obs_xyz(ixy, iz) > -1e33) THEN

                obsrec = obsrec + 1

                obs_p(obsrec) = obs_xyz(ixy,iz)
                thisobs%id_obs_p(1, obsrec) = idx_wet_in_all(rec) + sfields(id%sao)%off
                ocoord_p(1, obsrec) = lon_1d(ixy)
                ocoord_p(2, obsrec) = lat_1d(ixy)
                ocoord_p(3, obsrec) = scale_vert*REAL(iz)

                IF (rms_xyz(ixy, iz) >= rms_min) THEN
                   ivar_obs_p(obsrec) = rms_xyz(ixy, iz)
                ELSE
                   ivar_obs_p(obsrec) = rms_min
                   irms = irms + 1
                END IF
             END IF
          END DO
       END DO

       IF (irms>0) THEN
          WRITE(*,*) TRIM(sfields(id%sao)%crms),' < rms_min (',rms_min,'):', irms
       END IF

    ELSE
       ! No observation in process sub-domain
       ALLOCATE(obs_p(1))
       ALLOCATE(ocoord_p(thisobs%ncoord, 1))
       ALLOCATE(thisobs%id_obs_p(1, 1))
       ALLOCATE(ivar_obs_p(1))
    END IF


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_xyz, rms_xyz)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

  END SUBROUTINE init_dim_obs_EN4_sao



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_EN4_sao(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    ! Example: Observation operator for observed grid point values
    CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

  END SUBROUTINE obs_op_EN4_sao



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_EN4_sao(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    ! one can also set observation-specific values for the localization.
    USE mod_assimilation, &   
         ONLY: coords_l, locweight, cradius, cradius_z, iso_loc

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

! *** Local variables ***
    REAL :: lradius_3d(3), sradius_3d(3) ! non-isotropic localization radii

    WRITE(0,*) "RSE: ENTER init_dim_obs_l_EN4_sao"

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! For cradius and sradius:
    ! If these are defined as scalar values, isotropic localization is used.
    ! If these are vectors, nonisotropic localization is used
    !   (their length has to be equal to thisobs%ncoord)

    lradius_3d(1) = cradius
    lradius_3d(2) = cradius
    lradius_3d(3) = cradius_z

    sradius_3d = lradius_3d

    IF (.NOT.iso_loc) THEN
       ! Use non-isotropic localization
       WRITE(0,*) "RSE: ENTER init_dim_obs_l_EN4_sao noniso branch"
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, lradius_3d, sradius_3d, dim_obs_l)
       WRITE(0,*) "RSE: EXIT  init_dim_obs_l_EN4_sao noniso branch"
    ELSE
       ! Use isoptropic localization (usually needs scaling of vertical coordinates)
       WRITE(0,*) "RSE: ENTER init_dim_obs_l_EN4_sao iso branch"
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, cradius, cradius, dim_obs_l)
       WRITE(0,*) "RSE: EXIT  init_dim_obs_l_EN4_sao iso branch"
    END IF

    WRITE(0,*) "RSE: EXIT  init_dim_obs_l_EN4_sao"

  END SUBROUTINE init_dim_obs_l_EN4_sao



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE localize_covar_EN4_sao(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: cradius, locweight, sradius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_EN4_sao_pdafomi_TEMPLATE.F90: Apply covariance localization'

! *************************************
! *** Apply covariance localization ***
! *************************************

    ! Here one has to specify the three localization variables
    ! which can be different for each observation type.

    ! For cradius and sradius:
    ! If these are defined as scalar values, isotropic localization is used.
    ! If these are vectors, nonisotropic localization is used
    !   (their length has to be equal to thisobs%ncoord)

    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius, sradius, &
         coords_p, HP_p, HPH)

  END SUBROUTINE localize_covar_EN4_sao

END MODULE obs_EN4_sao_pdafomi
