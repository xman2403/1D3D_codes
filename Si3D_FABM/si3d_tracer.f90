!************************************************************************
                           MODULE si3d_tracer
!************************************************************************
!
!  Purpose: Inject tracer in the domain
!
!-------------------------------------------------------------------------

   USE si3d_types
   USE si3d_boundaryconditions


   IMPLICIT NONE
   SAVE

  
CONTAINS

!************************************************************************
SUBROUTINE read_tracer_from_txt_file()
!************************************************************************
!
!  Purpose: To initialize customized tracers
!
!------------------------------------------------------------------------

    INTEGER :: i100=100, i, j, ios, istat

    ! ... File open
    OPEN (UNIT=i100, FILE="si3d_inittracer.txt", STATUS="OLD", IOSTAT=ios)
    IF(ios /= 0) CALL open_error ( "Error opening si3d_inittracer.txt", ios )



    READ (UNIT=i100, FMT='(/(A))', IOSTAT=ios) title
    ! ... Read no. of columns, each with tracer injected
    READ (UNIT=i100, FMT='(///( 14X,I20))', IOSTAT=ios) ninj
    IF (ios  /= 0) CALL input_error ( ios, 702 )
    ! ... time instant (s) for tracer injected
    READ (UNIT=i100, FMT='(14X,G20.2)', IOSTAT=ios) tinj
    IF (ios  /= 0) CALL input_error ( ios, 703 )

    ! ... allocate space for settling velocity
    ALLOCATE (wsettle(ntr), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 704 )
    wsettle = 0.0;
    ! ... settling velocity for tracer injected
    READ (UNIT=i100, FMT='(14X,G20.2)', IOSTAT=ios) wsettle(iinj)
    IF (ios  /= 0) CALL input_error ( ios, 704 )

    IF ( iinj > ninj .OR. ninj < 1) THEN
        PRINT *, '*** ERROR ***'
        PRINT *, 'No. of groups cannot be > '
        PRINT *, 'No. of water columns with tracer injected'
        STOP
    ENDIF

    ! ... Allocate space for arrays holding location
    !     of tracer coordinates, concentration, and group id
    ALLOCATE ( iintr  (ninj), jintr   (ninj), &
                kintr1 (ninj), kintrm (ninj), cintr (ninj), ginj(ninj), STAT=istat)
    IF (istat /= 0) CALL allocate_error ( istat, 705 )

    ! ... Read in locations & characteristics of diffusers
    !     At this point, they are pressumed constants in time
    READ (UNIT=i100, FMT='(A)', IOSTAT=ios) commentline

    i = 0;
    DO j = 1, ninj
        IF (ios/= 0) CALL input_error ( ios, 706 )
        READ (UNIT=i100, FMT='(9X,5I4,G11.2)', IOSTAT=ios) &
                iintr(j), jintr(j), kintr1(j), kintrm(j), ginj(j), cintr(j)
        IF (ios/= 0) CALL input_error ( ios, 707 )
        IF (ginj(j).NE.i) i = ginj(j)
    ENDDO

    ! ... may need to be removed! Chris
    IF (i .NE. iinj) THEN
        PRINT *, '*** ERROR ***'
        PRINT *, 'No. of tracers causing groups'
        STOP
    ENDIF
    ! ...

    !.....Close tracer file.....
    CLOSE (UNIT=i100)

    !.....tracer activated
    isinjected = 0;

END SUBROUTINE read_tracer_from_txt_file

!************************************************************************
SUBROUTINE do_tracer_injection()
!************************************************************************
!
!  Purpose: To initialize customized tracers
!
!------------------------------------------------------------------------

    INTEGER :: inn, kinjtop, kinjbottom, kms, i, j, l


    DO inn=1, ninj

        ! ... Define i,j,l indexes
        i = iintr(inn);
        j = jintr(inn);
        l = ij2l(i,j);

        ! ... Define bottom cells & no. of layers
        kms = kmz(i,j) ;

        ! may cause errors with allocation -  REMOVED
!        IF (kintr1 (inn) < kms) THEN
!            kinjtop = kintr1 (inn);
!        ELSE
!            kinjtop = 2;
!        ENDIF
!
!        IF (kintrm(inn) <= kms) THEN
!            kinjbottom = kintrm(inn);
!        ELSE
!            kinjbottom = kms;
!        ENDIF

        ! ... assign bottom and top cells
        kinjtop = kintr1 (inn);
        kinjbottom = kintrm(inn);

        ! assign tracer
        tracer(kinjtop:kinjbottom,l,id_tracer_out(ginj(inn))) = cintr(inn);

        IF (idbg == 1) PRINT *,kinjtop, kinjbottom, ginj(inn), cintr(inn)

    ENDDO

    ! assign to the p-point
    tracerpp = tracer;

    ! mark as injected
    isinjected = 1;

END SUBROUTINE do_tracer_injection

!************************************************************************
SUBROUTINE do_tracer_vertical_movement()
!************************************************************************
!
!  Purpose: calculate settling rate by Chris
!   - under construction NOT USED
!
!------------------------------------------------------------------------

    !INTEGER :: kms, i, j, l, nt, k

    !sourcesink(k,l,nt)=-trwsettle*(tracerpp(k+1,l,nt)-tracerpp(k-1,l,nt))/2./ddz;

END SUBROUTINE do_tracer_vertical_movement

END MODULE si3d_tracer

