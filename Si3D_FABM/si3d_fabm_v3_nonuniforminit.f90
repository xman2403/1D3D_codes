!************************************************************************
                           MODULE si3d_fabm
!************************************************************************
!
!  Purpose: Procedures that implement routines in FABM
!
!-------------------------------------------------------------------------
   USE si3d_boundaryconditions
   USE si3d_types
   USE fabm
   USE fabm_config
   USE fabm_types

   IMPLICIT NONE
   SAVE

   TYPE (type_model) :: model
   TYPE (type_bulk_variable_id), ALLOCATABLE, DIMENSION(:) :: id_FABMtr_out ! index for output
   REAL(selected_real_kind(13)), ALLOCATABLE, DIMENSION(:,:), TARGET :: horizontal_wind_speed, SW2d, PAR2d
   REAL(selected_real_kind(13)), ALLOCATABLE, DIMENSION(: , : , : , :) , TARGET:: tracer3d, source3d
   REAL(selected_real_kind(13)), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: tempp, salinity1, depth3d, QswFr3d, PARflux3d, layerh3d, ssm3d
   REAL, ALLOCATABLE, DIMENSION(:,:), TARGET :: si3d_mask2d
   REAL, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: si3d_mask3d

   CHARACTER(LEN=256), ALLOCATABLE, DIMENSION(:) :: get_fabm_output_vars, get_fabm_init_vars
  
CONTAINS

!************************************************************************
SUBROUTINE input_fabm
!************************************************************************
!
!  Purpose: To read model details from file
!
!------------------------------------------------------------------------

    CHARACTER(LEN=64) :: var_input_file1 = "si3d_fabm_vars.txt"
    CHARACTER(LEN=64) :: var_input_file2 = "si3d_fabm_init_tracers.txt"
    INTEGER :: istat

    IF (idbg == 1) PRINT *, "entering reading nml"
    model = fabm_create_model_from_file(100);
    IF (idbg == 1) PRINT *, "exiting reading nml"

    ! Initialize (reads FABM configuration from fabm.yaml)
    ! After this the number of set of biogeochemical variables is fixed.
    ! (access variable metadata in model%state_variables, model%diagnostic_variables)
    ! ntr<=SIZE(model%state_variables)
    nfabmVar=SIZE(model%state_variables);

    ! output variables must be smaller than interior variables
    IF (ntr>nfabmVar) THEN
        PRINT *, "error, ntr>nfabmvar", ntr, ">", nfabmVar
        STOP
    ENDIF

    ! ... allocate variable name storage
    ALLOCATE(get_fabm_output_vars(ntr), get_fabm_init_vars(rdtr),STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 351 )

    ! ... read the list of output variables from file
    CALL ReadVarList(var_input_file1,ntr,get_fabm_output_vars)

    ! ... read the list of input tracers in si3d_init.txt
    CALL ReadVarList(var_input_file2,rdtr,get_fabm_init_vars)

    ! total number of tracers in FABM to be simulated
    nFABMtr_out = ntr;

    ! redefine ntr to the max number of tracers in FABM
    ntr = nfabmVar;


END SUBROUTINE input_fabm

!************************************************************************
SUBROUTINE init_fabm
!************************************************************************
!
!  Purpose: To initialize fabm
!
!------------------------------------------------------------------------

    INTEGER :: ivar
    INTEGER :: istat, j, k, i, nn

    ! ... allocate memory to a bulk matrix for tracers
    ALLOCATE (tracer3d(im1,jm1,km1,nfabmVar), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 301 )

    ALLOCATE (si3d_mask3d(im1,jm1,km1), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 302 )

    ALLOCATE (si3d_mask2d(im1,jm1), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 303 )

    ! ... allocate memory to a bulk matrix for tracers
    ALLOCATE (source3d(im1,jm1,km1,nfabmVar), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 304 )

    ! ... allocate memory to temperature matrix
    ALLOCATE (tempp(im1,jm1,km1), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 305 )

    ! ... allocate memory to temperature matrix
    ALLOCATE (salinity1(im1,jm1,km1), depth3d  (im1,jm1,km1), &
            & QswFr3d  (im1,jm1,km1), PARflux3d(im1,jm1,km1), &
            & layerh3d (im1,jm1,km1), ssm3d    (im1,jm1,km1), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 306 )

    ! ... allocate memory to settling velocity for tracers
    ALLOCATE (wVELtr(im1,jm1,km1,nfabmVar), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 307 )

    ! ... allocate memory to tracer IDs which are linked with the Si3D tracer (itr)
    ! clear id_tracer_out that defined earlier in the si3d_procedure
    DEALLOCATE (id_tracer_out)
    DEALLOCATE (id_init_var)
    ALLOCATE (id_FABMtr_out(nFABMtr_out), id_tracer_out(nFABMtr_out), &
            & id_init_var(rdtr), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 308 )

    ! ... allocate memory to variable name of FABM
    ALLOCATE (varnames_with_fabm(nfabmVar), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 309 )

    ! ... allocate memory to tracer IDs which are linked with the Si3D tracer (itr)
    ALLOCATE (horizontal_wind_speed(im1,jm1), SW2d(im1,jm1), &
            & PAR2d (im1,jm1), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 310 )

    IF (idbg == 1) PRINT *, "entering set domain"
    ! Provide extents of the spatial domain.
    CALL fabm_set_domain(model,im1,jm1,km1)
    IF (idbg == 1) PRINT *, "exiting set domain"

    IF (idbg == 1) PRINT *, "entering set indexes"
    ! Specify vertical index of surface and bottom
    CALL model%set_surface_index(1)
    CALL model%set_bottom_index(km1)
    IF (idbg == 1) PRINT *, "exiting set indexes"

    IF (idbg == 1) PRINT *, "entering set mask"
    ! send mask to FABM
    CALL define_3Dmask_for_fabm(si3d_mask3d,si3d_mask2d)
    CALL fabm_set_mask(model,si3d_mask3d,si3d_mask2d)
    IF (idbg == 1) PRINT *, "exiting set mask"

    ! allocate scalars to bulk matrix
    tracer3d=0;
    id_tracer_out = 0;
    id_init_var = 0;
    IF (idbg == 1) PRINT *, "entering link vars"
    ! Send pointers to state variable data in FABM
    DO ivar=1,ntr

        CALL fabm_link_interior_state_data(model,ivar,tracer3d(:,:,:,ivar))

        ! ... get all allocated variable names
        varnames_with_fabm(ivar)=model%state_variables(ivar)%name;
        PRINT *, "State variable name", ivar, varnames_with_fabm(ivar)

        ! assign the output vars' id in the solver matrix to id_tracer_out
        DO i = 1, nFABMtr_out
            ! ... Assign indexes to output variables
            IF (TRIM(varnames_with_fabm(ivar)) == get_fabm_output_vars(i)) THEN
                id_tracer_out (i) = ivar;
            ENDIF
            ! ... Assign indexes to output variables (should be always 1)
            IF (TRIM(varnames_with_fabm(ivar)) == 'aed_oxygen_oxy') THEN
                LDO = ivar;
            ENDIF
        ENDDO

        ! assign the init tracers' id in the solver matrix to id_init_var
        DO j = 1, rdtr
            ! ... Assign indexes to init tracers
            IF (TRIM(varnames_with_fabm(ivar)) == get_fabm_init_vars(j)) THEN
                id_init_var (j) = ivar;
            ENDIF
        ENDDO

    ENDDO

    ! ... check if all indexes are assigned
    DO ivar = 1, nFABMtr_out
        IF (id_tracer_out(ivar) == 0) THEN
            PRINT *, "Output ID unassigned: ", ivar
            STOP
        ENDIF
    ENDDO
    DO ivar = 1, rdtr
        IF (id_init_var(ivar) == 0) THEN
            PRINT *, "Tracer input ID unassigned: ", ivar
            STOP
        ENDIF
    ENDDO


    ! convert temperature and depth matrix
    IF (idbg == 1) PRINT *, "entering temperature conversion"
    CALL kl2xyz(sal,   tempp)
    CALL kl2xyz(zfromt,depth3d)
    CALL kl2xyz(QswFr, QswFr3d)
    CALL kl2xyz(hp,layerh3d)

    ! ... start with linking bulk and surface data

    ! Transfer pointer to environmental data
    ! Array temp with extents nx,ny,nz is assumed to be declared and accessible.
    ! Do this for all variables on FABM's standard variable list that the model can provide.
    ! For this list, visit http://fabm.net/standard_variables
    IF (idbg == 1) PRINT *, "entering link temperature"
    CALL model%link_bulk_data(standard_variables%temperature,tempp)
    ! assume no salinity for fresh waters (3D)
    salinity1 = 0;
    CALL model%link_bulk_data(standard_variables%practical_salinity,salinity1)

    ! Link depth 3D to FABM
    CALL model%link_bulk_data(standard_variables%depth,depth3d)

    ! Link eta 3D to FABM
    CALL model%link_bulk_data(standard_variables%attenuation_coefficient_of_shortwave_flux,QswFr3d)

    ! Link PAR light 3D to FABM
    PARflux3d = 0; ! w/m2
    CALL model%link_bulk_data(standard_variables%downwelling_photosynthetic_radiative_flux,PARflux3d)

    ! Link layer thickness 3D to FABM
    CALL model%link_bulk_data(standard_variables%cell_thickness,layerh3d)

    ssm3d = 0; ! mg/L
    CALL model%link_bulk_data(standard_variables%mass_concentration_of_suspended_matter,ssm3d)

    ! ... start with surface data linking
    ! Initilize wind speed - no wind
    horizontal_wind_speed = 0;
    ! Link wind speed 2D to FABM
    CALL model%link_horizontal_data(standard_variables%wind_speed,horizontal_wind_speed)

    ! Link short wave radiation (assume uniform)
    SW2d = 0;
    CALL model%link_horizontal_data(standard_variables%surface_downwelling_shortwave_flux,SW2d)

    ! Link PAR light radiation to FABM
    PAR2d = 0;
    CALL model%link_horizontal_data(standard_variables%surface_downwelling_photosynthetic_radiative_flux,PAR2d)

    IF (idbg == 1) PRINT *, "entering check model"
    ! Check whether FABM has all dependencies fulfilled
    ! (i.e., whether all required calls for fabm_link_*_data have been made)
    CALL fabm_check_ready(model)

    IF (idbg == 1) PRINT *, "entering initialising model"
    ! Initialise tracers
    DO k = 1, km1; DO j = 1, jm1
        call fabm_initialize_state(model,1,im1,j,k)
    ENDDO; ENDDO
    IF (idbg == 1) PRINT *, "exiting initialising model"

    IF (idbg == 1) PRINT *, "entering tracer assigning from FABM"
    DO ivar = 1, ntr
        ! DO in FABM will be overwriten | modified to reassign below
        !IF (ivar == LDO) CYCLE
        ! FABM to tracers
        CALL xyz2kl(tracer3d(:,:,:,ivar),tracer(:,:,ivar)) ! O2 mmol/m**3
    ENDDO

    IF (idbg == 1) PRINT *, "entering tracer overwriten by file"
    ! over write some part in SUBROUTINE InitializeScalarFields
    DO nn = 1, rdtr
        DO k = 1, km1
          tracer(k,:,id_init_var(nn)) = ScalarProfile(k,nn+1)
        ENDDO
    END DO ! ... End loop over tracers

    tracerpp = tracer;

    IF (idbg == 1) PRINT *, "end init_fabm"


END SUBROUTINE init_fabm

!************************************************************************
SUBROUTINE define_3Dmask_for_fabm(si3d_mask3d,si3d_mask2d)
!************************************************************************
!
!  Purpose: To pass the 3D mask into fabm
!
!------------------------------------------------------------------------

    REAL, DIMENSION(im1,jm1), INTENT(INOUT) :: si3d_mask2d
    REAL, DIMENSION(im1,jm1,km1), INTENT(INOUT) :: si3d_mask3d
    INTEGER :: i, j, kms, k1s

    ! ... set all to zeros (assume all cells are dry)
    si3d_mask2d = 0;
    si3d_mask3d = 0;

    ! ... Loop over the domain and check wet cells
    DO i = 1 , im1; DO j = 1 , jm1
        IF (mask2d(i,j)) THEN
            kms=kmz(i,j);
            k1s=k1z(i,j);
            ! set wet celll for the wet column
            si3d_mask2d(i,j        ) = 1;
            si3d_mask3d(i,j,k1s:kms) = 1;
        ENDIF
    ENDDO; ENDDO

END SUBROUTINE define_3Dmask_for_fabm


!************************************************************************
SUBROUTINE do_fabm_aed
!************************************************************************
!
!  Purpose: To implement fabm
!
!------------------------------------------------------------------------

    INTEGER :: itr, j, k
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,nfabmVar) :: flux_sf,flux_bt
    REAL(selected_real_kind(13)) :: sms_sf(im1,jm1,SIZE(model%surface_state_variables))
    REAL(selected_real_kind(13)) :: sms_bt(im1,jm1,SIZE(model%bottom_state_variables))
    REAL :: rs

    ! Update tracer3d with the current states
    DO itr = 1, ntr
        CALL kl2xyz(tracer(:,:,itr),tracer3d(:,:,:,itr)) ! O2 mmol/m**3
    ENDDO

    ! update temperature
    CALL kl2xyz(sal,tempp)

    ! update bulk data
    CALL kl2xyz(hp,layerh3d)
    CALL kl2xyz(zfromt,depth3d)
    CALL kl2xyz(QswFr, QswFr3d)

    ! update wind speed
    horizontal_wind_speed = SQRT(uair**2 + vair**2);

    ! update radiation
    rs = parab(0.,thrs,surfbc1(2,:),dtSurfbc/3600);
    ! short wave
    SW2d = rs;
    ! PAR lights
    PARflux3d = rs * swR_to_PAR;    PARflux3d = PARflux3d*QswFr3d;
    PAR2d = rs*swR_to_PAR;

    ! update domain mask
    CALL define_3Dmask_for_fabm(si3d_mask3d,si3d_mask2d)

    IF (idbg == 1) PRINT *, "update source terms"
    ! fabm update pelagic sources
    source3d = 0;
    wVELtr = 0;
    DO k = 1, km1; DO j = 1, jm1
        CALL fabm_do(model,1,im1,j,k,source3d(:,j,k,:))
        ! get vertical movement and use for advection equation of tracers later - added on 22/03/2016 Chris
        CALL fabm_get_vertical_movement(model,1,im1,j,k,wVELtr(:,j,k,:))
    ENDDO; ENDDO

    ! fabm update surface and bottom fluxes - added on 22/03/2016 Chris
    flux_sf = 0; flux_bt = 0; sms_sf = 0; sms_bt = 0;
    DO j = 1, jm1
        CALL fabm_do_surface(model,1,im1,j,flux_sf(:,j,:),sms_sf(:,j,:))
        CALL fabm_do_bottom(model,1,im1,j,flux_bt(:,j,:),sms_bt(:,j,:))
    ENDDO

    ! assign surface source term to source3d
    CALL layer2sourceterm(flux_sf,k1z,source3d)
    ! assign bottom source term to source3d
    CALL layer2sourceterm(flux_bt,kmz,source3d)

    IF (idbg == 1) PRINT *, "finished update"

    IF (idbg == 1) PRINT *, "entering tracer source term allocation"
    ! assign values to source and sink storage, await for imtracer
    DO itr = 1, ntr
        CALL xyz2kl(source3d(:,:,:,itr),sourcesink(:,:,itr));
    ENDDO



END SUBROUTINE do_fabm_aed

!************************************************************************
SUBROUTINE layer2sourceterm(flux,kxz,source3d)
!************************************************************************
!
!  Purpose: To convert surface flux into source term
!
!------------------------------------------------------------------------
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,nfabmVar), INTENT(IN) :: flux
    INTEGER, DIMENSION(im1,jm1), INTENT(IN) :: kxz
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,km1,nfabmVar), INTENT(INOUT) :: source3d
    INTEGER :: l, i, j, kss

    DO l= 1, lm
        ! ... 3D-(i,j) indexes for l
        i = l2i(l); j = l2j(l);
        ! ... get surface/bottom index
        kss = kxz(i,j);

        ! assign values
        source3d(i,j,kss,:)=source3d(i,j,kss,:) + flux(i,j,:)/hp(kss,l)
    ENDDO

END SUBROUTINE layer2sourceterm

!************************************************************************
SUBROUTINE kk2xyz(var1d,var3d)
!************************************************************************
!
!  Purpose: To convert matrices from 1d to 3d
!
!------------------------------------------------------------------------

    REAL, DIMENSION(km1), INTENT(IN) :: var1d
    INTEGER :: l
    INTEGER :: i, j
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,km1), INTENT(OUT) :: var3d

    var3d = 0.0;

    DO l= 1, lm

        ! ... 3D-(i,j) indexes for l
        i = l2i(l); j = l2j(l);
        ! assign variable for the water column
        var3d(i,j,:) = var1d(:);

    ENDDO

END SUBROUTINE kk2xyz

!************************************************************************
SUBROUTINE kl2xyz(var2d,var3d)
!************************************************************************
!
!  Purpose: To initialize fabm
!
!------------------------------------------------------------------------

    REAL, DIMENSION(km1,lm1), INTENT(IN) :: var2d
    INTEGER :: l
    INTEGER :: i, j
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,km1), INTENT(OUT) :: var3d

    var3d = 0.0;

    DO l= 1, lm

        ! ... 3D-(i,j) indexes for l
        i = l2i(l); j = l2j(l);
        ! assign variable for the water column
        var3d(i,j,:) = var2d(:,l);

    ENDDO

END SUBROUTINE kl2xyz

!************************************************************************
SUBROUTINE xyz2kl(var3d,var2d)
!************************************************************************
!
!  Purpose: To initialize fabm
!
!------------------------------------------------------------------------
    REAL(selected_real_kind(13)), DIMENSION(im1,jm1,km1), INTENT(IN) :: var3d
    INTEGER :: l, kk
    INTEGER :: i, j, kms, k1s
    REAL, DIMENSION(km1,lm1), INTENT(OUT) :: var2d

    var2d = 0.0;

    DO i= 1, im1 ; DO j= 1, jm1

        ! ... 3D indexes to l
        l=ij2l(i,j);
        ! assign variable for the water column
        var2d(:,l) = var3d(i,j,:);

    ENDDO; ENDDO

END SUBROUTINE xyz2kl

!************************************************************************
SUBROUTINE ReadVarList(var_input_file,n_out,fabmvar_out)
!************************************************************************
!
!  Purpose: To initialize fabm
!
!------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: n_out
    CHARACTER(LEN=64), INTENT(IN) :: var_input_file
    CHARACTER(LEN=256),DIMENSION(n_out), INTENT(OUT) :: fabmvar_out
    CHARACTER(LEN=256) :: fline, title
    INTEGER :: i100=100, i, ios, istat, nn
    LOGICAL :: done

    ! ... File open
    OPEN (UNIT=i100, FILE=var_input_file, STATUS="OLD", IOSTAT=ios)
    IF(ios /= 0) CALL open_error ( "Error opening "//var_input_file, ios )

    ! ... Read headers
    READ (UNIT=i100, FMT='(//(A)/)',IOSTAT=ios) title
    IF (ios /= 0) CALL input_error ( ios, 350 )

    ! ... initialize
    done=.FALSE.
    i=0;

    IF (idbg == 1) PRINT *, "entering variables loop when opening si3d_fabm_vars"
    DO WHILE (.NOT. done)
        READ (UNIT=i100, FMT='(A)', IOSTAT=ios) fline
        IF (ios /= 0) CALL input_error ( ios, 352 )

        IF (TRIM(fline)=='' .OR. i==n_out) THEN
            CLOSE (UNIT=i5)
            EXIT
        ENDIF

        IF (idbg == 1) PRINT *, "appending values to fabmvar_out"
        ! assign number and variables
        i=i+1
        fabmvar_out(i)=TRIM(fline);
    ENDDO



END SUBROUTINE ReadVarList


END MODULE si3d_fabm

