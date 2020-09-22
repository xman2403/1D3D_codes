!************************************************************************
                          MODULE si3d_procedures
!************************************************************************
!
!  Purpose: Procedures for the semi-implicit 3-D (si3d) hydrodynamic
!           model.
!
!-------------------------------------------------------------------------

   USE si3d_types
   USE si3d_boundaryconditions
   USE si3d_mixing
   USE si3d_ecomod
   USE si3d_fabm
   USE si3d_tracer

   IMPLICIT NONE
   SAVE
  
CONTAINS

!************************************************************************
SUBROUTINE input
!************************************************************************
!
!  Purpose: To read si3d model input parameters.
!
!------------------------------------------------------------------------

   !.....Local variables.................................................
   CHARACTER(LEN=12) :: input_file = "si3d_inp.txt"
   INTEGER :: ios, nn, i, j, istat

   !.....Open input parameter file.....
   OPEN (UNIT=i5, FILE=input_file, STATUS="OLD", IOSTAT=ios)
   IF (ios /= 0) CALL open_error ( "Error opening "//input_file, ios )


   !.....Read header record containing comments about run................
   READ (UNIT=i5, FMT='(/(A))', IOSTAT=ios) title
   IF (ios /= 0) CALL input_error ( ios, 1)

   !.....Read & define start date of the run.............................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) iyr0,imon0,iday0,ihr0 
   IF (ios /= 0) CALL input_error ( ios, 2 )
   CALL compute_date (0.0) ! idt real

   !.....Read space-time domains, cell size & time step .................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) xl,yl,zl,tl,idx,idy, &
   & idz,dzmin, datadj, zetainit,idt, ibathyf, xslope, yslope
   IF (ios /= 0) CALL input_error ( ios, 3 )

   ! ... Read parameters controlling solution algorithm .................
   READ (UNIT=i5,FMT='(///(14X,G20.2))',IOSTAT=ios) itrap,niter,ismooth, &
       & beta, iturb, Av0, Dv0, iadv, itrmom, ihd, Ax0, Ay0, f, theta,   &
       & ibc,isal, itrsca, cd, ifsurfbc, dtsurfbc, cw, wa, phi, radischm,&
       & rdtr, isseddep, idbg
   IF (ios /= 0) CALL input_error ( ios, 4 )

   !.....Read node numbers for time series output .......................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) ipt
   IF (ios /= 0) CALL input_error ( ios, 5 )
   IF (ipt > 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) nnodes
     IF (ios /= 0) CALL input_error ( ios, 5 )
     READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios) (inode(nn), nn = 1, nnodes)
     IF (ios /= 0) CALL input_error ( ios, 5 )
     READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios) (jnode(nn), nn = 1, nnodes)
     IF (ios /= 0) CALL input_error ( ios, 5 )
   ENDIF

   ! ... Read nodes for horizontal plane output ...........................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) iop
   IF (ios /= 0) CALL input_error ( ios, 6 )
   IF (iop /= 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) n_planes
     IF (ios /= 0) CALL input_error ( ios, 6 )
     IF (n_planes > max_planes) THEN
       PRINT *,'ERROR: # of planes requested > maximum allowed' 
       STOP
     END IF    
     DO j = 1, n_planes
       ! ... Read number of cells in X-section j
       READ (UNIT=i5, FMT='(14X,I20)' , IOSTAT=ios) p_out(j)
     ENDDO
   ENDIF

   ! ... Read nodes for vertical plane output ...........................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) iox
   IF (ios /= 0) CALL input_error ( ios, 7 )
   IF (iox /= 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) n_sections
     IF (ios /= 0) CALL input_error ( ios, 7 )
     IF (n_sections > max_sections) THEN
       PRINT *,'ERROR: # of sections requested > maximum allowed' 
       STOP
     END IF    
     DO j = 1, n_sections
       ! ... Read number of cells in X-section j
       READ (UNIT=i5, FMT='(/14X,I20)' , IOSTAT=ios) n_section_cells(j)
       IF (ios /= 0) CALL input_error ( ios, 7 )
       ! ... Read i coordinates for cells in X-section j
       READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios)                       &
       &    (xinode(j, nn), nn = 1, n_section_cells(j) )
       IF (ios /= 0) CALL input_error ( ios, 7 )
       ! ... Read j coordinates for cells in X-section j
       READ (UNIT=i5, FMT='(14X,10I5)', IOSTAT=ios)                       &
       &    (xjnode(j, nn),nn=1,n_section_cells(j))
       IF (ios /= 0) CALL input_error ( ios, 7 )
     END DO
   ENDIF

   !! ... Read toggles for 3D output  .....................................
   READ (UNIT=i5, FMT='(///(14X,I20))', IOSTAT=ios) ipxml
   IF (ios /= 0) CALL input_error ( ios, 7 )
   READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) itspf
   IF (ios /= 0) CALL input_error ( ios, 7 )

   !.....Read info on open boundaries.....................................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) nopen
   IF (ios /= 0) CALL input_error ( ios, 8 )
   IF (nopen > 0) THEN
      READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) dtsecopenbc
      IF (ios /= 0) CALL input_error ( ios, 8 )
      DO nn = 1, nopen
         READ (UNIT=i5, FMT='(/(14X,I20))', IOSTAT=ios) iside(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) itype(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) isbc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jsbc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iebc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jebc(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
      END DO
   END IF
 
   !.....Read info to output nested grid boundaries ......................
   IF (ioNBTOGGLE > 0 ) THEN
     READ (UNIT=i5,  FMT='(///14X,I20)', IOSTAT=ios) nxNBO
     IF (ios /= 0) CALL input_error ( ios, 8 )
     IF (nxNBO > 0) THEN
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) ioNBO
       IF (ios /= 0) CALL input_error ( ios, 8 )
       READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) xxNBO
       IF (ios /= 0) CALL input_error ( ios, 8 )
       DO nn = 1, nxNBO
         READ (UNIT=i5, FMT='(/(14X,I20))', IOSTAT=ios) isdNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) isbcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jsbcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iebcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
         READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) jebcNBO(nn)
         IF (ios /= 0) CALL input_error ( ios, 8 )
       END DO
     END IF
   END IF

   !.... Read info for tracers ...........................................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) ntr
   IF (ios  /= 0) CALL input_error ( ios, 9 )
   iotr = 0; ! Default value
   IF (ntr > 0) THEN
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) ecomod
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iotr
     IF (ios  /= 0) CALL input_error ( ios, 9 )
     !.... Read info for the tracer group of interest (match ginj()) ......
     READ (UNIT=i5, FMT='(14X,I20)', IOSTAT=ios) iinj
     IF (ios  /= 0) CALL input_error ( ios, 701 )
     ! read tracer profiles
     IF (iinj > 0) THEN
         CALL read_tracer_from_txt_file
     ENDIF
   ENDIF

   !.... Read info for plume models & oxygenation systems ................
   READ (UNIT=i5, FMT='(///14X,I20)', IOSTAT=ios) iopss
   IF (ios  /= 0) CALL input_error ( ios, 10 )
   IF (iopss ==0) k4sod=0.0
   IF (iopss > 0) THEN
     ! ... Read no. of devices, each with its own inflow/outflow rate
     READ (UNIT=i5, FMT='( 14X,I20)', IOSTAT=ios) npssdev
     IF (ios  /= 0) CALL input_error ( ios, 10 )
     IF ( npssdev > iopss .OR. npssdev < 1) THEN
       PRINT *, '*** ERROR ***'
       PRINT *, 'No. of devices creating point sources & sinks cannot be > '
       PRINT *, 'No. of water columns with point sources & sinks'
       STOP
     ENDIF
     ! Read time in seconds between consecutive records from time files
     READ (UNIT=i5, FMT='(14X,G20.2)', IOSTAT=ios) dtsecpss

     ! ... Allocate space for arrays holding location 
     !     of point sources and sinks and devices
     ALLOCATE ( ipss  (iopss), jpss   (iopss), &
                iodev (iopss), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 100 )

     ! ... Read in locations & characteristics of diffusers 
     !     At this point, they are pressumed constants in time
     READ (UNIT=i5, FMT='(A)', IOSTAT=ios) commentline
     i = 0;
     DO j = 1, iopss
       IF (ios/= 0) CALL input_error ( ios, 10)
       READ (UNIT=i5, FMT='(9X,3I4)', IOSTAT=ios)  &
            ipss(j), jpss(j), iodev(j)
       IF (ios/= 0) CALL input_error ( ios, 10 )
       IF (iodev(j).NE.i) i = iodev(j)
     END DO
     IF (i .NE. npssdev) THEN
       PRINT *, '*** ERROR ***'
       PRINT *, 'No. of devices causing point sources'
       STOP
     ENDIF
   ENDIF

   ! ... Input instructions & parameters controlling the solution of 
   !     the tracer transport equations. 
   IF (ntr > 0) THEN

     nFABMtr_out=ntr;
     ! ... temperally allocate memoery to output storage (can be overwritten by the Si3D-FABM module)
     ALLOCATE (id_tracer_out(nFABMtr_out),id_init_var(nFABMtr_out), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 11 )
     ! ... Assign default values
     DO i = 1, nFABMtr_out
         id_tracer_out(i) = i;
         id_init_var(i)=i;
     ENDDO


     ! initialise DO location in tracers (to be rectify later)
     LDO = 1;
     ! .... Define other input variables
     SELECT CASE (ecomod)
     CASE (-1) ! Tracer Cloud Releases 
        PRINT *, 'Tracer Cloud Modelling activated'
        CALL trcinput 
     CASE (1) ! Water quality routines
        PRINT *, 'Water Quality Model activated'
        CALL WQinput
     CASE (2) ! Size Structure distribution
        PRINT *, 'Size Structure Model activated'
        CALL szinput
     CASE (3) ! Sediment transport routines
        PRINT *, 'Sediment transport activated'
        CALL sdinput
     CASE (4) ! FABM - AED
        PRINT *, 'AED activated'
        CALL input_fabm
        ! ntr has been redefined as the max number of tracers in FABM
     END SELECT


     ! ... Allocate space for some arrays - they are initialized
	 !     only if ecomod < 0, but used allways in determining
	 !     if the tracer transport equation is used or not in subroutine fd.
     ALLOCATE ( trct0(ntr), trcpk(ntr), trctn(ntr), &
                trcx0(ntr), trcy0(ntr), trcz0(ntr), &
                trcsx(ntr), trcsy(ntr), trcsz(ntr), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 121 )

     ! .... Initialize trct0 and trctn to default values
	 trct0 = 1E7;
	 trctn =   0;

     ! ... Read sources & sinks specifications -
     CALL PointSourceSinkInput

   ENDIF


   !.....Close input file.....
   CLOSE (UNIT=i5)

   !.....Define frequently used numerical constants and coefficients.....
   dt=idt; dx=idx; dy=idy; ddz=idz; twodt=2.*dt; dtdx=dt/dx; dtdy=dt/dy
   gdtdx=g*dtdx; gdtdy=g*dtdy; gdt2dx2=gdtdx*dtdx; gdt2dy2=gdtdy*dtdy
   dxdy = dx * dy; ! Smagorinsky
   !gthx=gdtdx*theta; gthy=gdtdy*theta; gth1x=gdtdx*2.*(1.-theta)
   !gth1y=gdtdy*2.*(1.-theta); 
   cwind=2.*dt*cw*rhoair*wa*wa; alp4=(1.-alp)/4.
   twodx=2.*dx; twody=2.*dy; fourdx=2.*twodx; fourdy=2.*twody
   dxdx=dx*dx; dydy=dy*dy; twodxdx=2.*dxdx; twodydy=2.*dydy; beta2=beta/2.
   chi1=1.-chi; twochi1=2.*chi1
   im=nint(xl/dx)+1; im1=im+1; i1=2; ndx=im1-i1
   jm=nint(yl/dy)+1; jm1=jm+1; j1=2; ndy=jm1-j1
   nts=tl/dt+.5; apxml = ABS(ipxml)
   isec0=REAL(ihr0)/100.*3600.; ! dt_min=idt/60. ! idt real

   ! ... Generate grid dimensions in Z-direction 
   CALL ZGridDimensions

END SUBROUTINE input

!************************************************************************
SUBROUTINE ZGridDimensions
!************************************************************************
!
!  Purpose: To read in depths of levels between consecutive layers
!
!------------------------------------------------------------------------

   !.....Local variables.................................................
   CHARACTER(LEN=14) :: input_file = "si3d_layer.txt"
   INTEGER :: ios, k, istat

   ! ... Variable layer thickness
   IF (ibathyf < 0 ) THEN

     !.....Open input file.....
     OPEN (UNIT=i5, FILE=input_file, STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening "//input_file, ios )

     !.....Read information ...
     ! Skip over first line (header)  
     READ (UNIT=i5, FMT='(//)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 40 )
   
     ! Read number of layers from second header record
     READ (UNIT=i5, FMT='(10X,I11)', IOSTAT=ios) km1
     IF (ios /= 0) CALL input_error ( ios, 41 )
  
     !..... Allocate space for zlevel array
     ALLOCATE (zlevel(1:km1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 0 )

     ! .... Read array with levels to layer interfaces 
     DO k = 1, km1
       READ (UNIT=i5, FMT='(10X,G11.2)', IOSTAT=ios) zlevel(k)
       IF (ios /= 0) CALL input_error ( ios, 42 )
     END DO

     ! ... Generate grid dimensions in z-direction
     km  = km1 - 1  ;
     k1  = 2        ;
     ndz = km1 - k1 ;

   ! ... Constant layer thickness (original si3d) 
   ELSE

     ! ... Generate grid dimensions in z-direction
     !km=CEILING((zl-dzmin)/ddz)+1
     km=CEILING((zl)/ddz)+1
     km1=km+1; k1=2; ndz=km1-k1

     !..... Allocate space for zlevel array
     ALLOCATE (zlevel(1:km1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 0 )

     !.....Initialize arrays with levels to layer interfaces 
     zlevel(k1) = 0.0; 
     DO k = k1+1, km1
       zlevel(k)=zlevel(k-1)+ddz
     END DO
     zlevel(k1) = -100. 
     zlevel(1 ) = -100.

   ENDIF


END SUBROUTINE ZGridDimensions

!************************************************************************
SUBROUTINE AllocateSpace
!************************************************************************
!
!  Purpose: To allocate space for model arrays at a size determined during 
!           execution. These arrays are 'permanently'
!           allocated for the entire duration of a model run.
!
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: istat, one=1

   !.....Allocate logical mask arrays.....
   ALLOCATE ( mask2d(im1,jm1    ), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 1 )

   !..... Allocate geometry arrays
   ALLOCATE (kmz   (im1,jm1), k1z(im1,jm1), &
             k1u   (im1,jm1), k1v(im1,jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 2 )

   !.....Allocate matrix solution arrays.....
   ALLOCATE ( ubar(1), rparm(16), iparm(25), ip(1), jp(1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 3 )

   !.....Allocate other miscellaneous arrays.....Map 2D-l into 3D-(i,j) indexes
   ALLOCATE ( ds(km), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 4 )

   !.....Allocate arrays with variables at zeta-pts in 3D space &
   !     arrays used in solution procedures....
   ALLOCATE ( s   (im1,jm1), sp  (im1,jm1), &
            & spp (im1,jm1),                &
            & sx  (im1,jm1), sy  (im1,jm1), &
            & dd  (im1,jm1), qq  (im1,jm1), &
            & eagx(im1,jm1), earx(im1,jm1), &
            & eagy(im1,jm1), eary(im1,jm1), &
            & rr  (im1,jm1), hhs (im1,jm1), &
            & hhu (im1,jm1), hhv (im1,jm1), & 
            & uair(im1,jm1), vair(im1,jm1), &
            & cdw (im1,jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 5 )
 
   ! .... Allocate arrays used in model output
   ALLOCATE(  uout (km1) , vout (km1) , wout(km1),  &
            & Avout(km1) , Dvout(km1) , sal1(ndz),  &
            & uhout(km1) , scout(km1),  trout(km1,ntrmax), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 7 )

   ! ...  Allocate space for output routines
   ALLOCATE ( interior_plane_points   ( n_planes  ),  &
            & interior_section_points ( n_sections), STAT = istat  )
   IF (istat /= 0) CALL allocate_error ( istat, 8 )

   ! .... Allocate arrays used to store velocity boundary conditions
   IF (nopen > 0) THEN 
     ALLOCATE(  uhWB  (km1,jm1), uhEB  (km1,jm1),   &
                huWB  (km1,jm1), huEB  (km1,jm1),   &
                vhSB  (km1,im1), vhNB  (km1,im1),   &
                hvSB  (km1,im1), hvNB  (km1,im1),   &
                uhWBpp(km1,jm1), uhEBpp(km1,jm1),   &
                huWBpp(km1,jm1), huEBpp(km1,jm1),   &
                vhSBpp(km1,im1), vhNBpp(km1,im1),   &
                hvSBpp(km1,im1), hvNBpp(km1,im1),   &
                uhWBp (km1,jm1), uhEBp (km1,jm1),   &
                huWBp (km1,jm1), huEBp (km1,jm1),   &
                vhSBp (km1,im1), vhNBp (km1,im1),   &
                hvSBp (km1,im1), hvNBp (km1,im1), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 9 )
   ENDIF

END SUBROUTINE AllocateSpace

!************************************************************************ 
SUBROUTINE AllocateSpace2
!************************************************************************
!
!  Purpose: To allocate space for model arrays at a size determined during 
!           execution. These arrays are 'permanently'
!           allocated for the entire duration of a model run.
!
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   ! ... Local variables
   INTEGER:: istat

   ALLOCATE ( ij2l(im1,jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 10 )
 
   ALLOCATE ( l2i(lm1), l2j(lm1), & 
   &          lEC(lm1), lWC(lm1), &
   &          lNC(lm1), lSC(lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 11 )

   !.....Allocate arrays at u-pts.....
   ALLOCATE ( uh(km1,lm1), uhp(km1,lm1), uhpp(km1,lm1),      & 
            & u (km1,lm1), up (km1,lm1), upp (km1,lm1),      &
            & ex(km1,lm1), agx(km1,lm1), arx (km1,lm1),      &
            & hu(km1,lm1), hup(km1,lm1), hupp(km1,lm1),STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 12 )

   !.....Allocate arrays at v-pts.....
   ALLOCATE ( vh(km1,lm1), vhp(km1,lm1), vhpp(km1,lm1),      & 
            & v (km1,lm1), vp (km1,lm1), vpp (km1,lm1),      &
            &              agy(km1,lm1), ary (km1,lm1),      &
            & hv(km1,lm1), hvp(km1,lm1), hvpp(km1,lm1),STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 13 )

   !.....Allocate arrays at pressure-pts.....
   ALLOCATE ( h   (km1,lm1), hp  (km1,lm1), hpp   (km1,lm1), & 
            & sal (km1,lm1), salp(km1,lm1), salpp (km1,lm1), &
            & rhop(km1,lm1), zfromt(km1,lm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 14 )

   !.....Allocate arrays at vertical interfaces between pressure-pts.....
   ALLOCATE ( Av  (km1,lm1), Dv  (km1,lm1),            &
              wp  (km1,lm1), STAT=istat    )
   IF (istat /= 0) CALL allocate_error ( istat, 15 )

   ! ....Allocate arrays at pressure points  
   ALLOCATE ( QswFr     (km1,lm1), &
           &  HeatSource(km1,lm1), STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 16 )
   
   !.....Allocate horizontal diffusion, th, and th1 arrays.....
   ALLOCATE (haypp(km1,lm1),haxpp(km1,lm1), &
             hdypp(km1,lm1),hdxpp(km1,lm1), &
             th   (km1,lm1),th1  (km1,lm1),STAT=istat)
   IF (istat /=0) CALL allocate_error ( istat, 17 )

   !.... Allocate space for Higher-Order turbulence models......
   IF (iturb > 0) THEN   
     ALLOCATE ( q2 (km1,lm1), q2p (km1,lm1 ), q2pp (km1,lm1),   &
                q2l(km1,lm1), q2lp(km1,lm1 ), q2lpp(km1,lm1),   &
                dsT(km1)    , aaT (3, km1+1), sal1T(ndz+1), STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 18 )
   ELSE IF (iturb < 0) THEN                	 	! si3dgotm
     ALLOCATE ( si3dtke (km1,lm1), &	    		! si3dgotm
                si3deps (km1,lm1), &        		! si3dgotm
                si3dlen (km1,lm1), STAT = istat )	! si3dgotm
     IF (istat /= 0) CALL allocate_error ( istat, 18 )  ! si3dgotm
   ENDIF

   ! .... Allocate arrays used in advection for tracers & plumes & 
   !      when scalar transport is done with flux-limiters ..............
   IF (iopss > 0 .OR. ntr > 0 .OR. itrsca >= 4 ) THEN
     ALLOCATE(  fluxX (km1,lm1),fluxY(km1,lm1),fluxZ(km1,lm1),  &
                STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 19 )
   ENDIF

   ! .... Allocate arrays for tracers ...................................
   IF (ntr > 0 ) THEN
     ALLOCATE(  tracer     (km1,lm1,ntr),  &
                tracerpp   (km1,lm1,ntr),  &
                sourcesink (km1,lm1,ntr),  &
                STAT = istat )
     IF (istat /= 0) CALL allocate_error ( istat, 20 )
   ENDIF 

   ! .... Allocate arrays used in oxygenation simulations - these are allocated
   !      here and not in pssbc0 or bcpss0, since when I call to those other 
   !      routines (at the time of reading IN, km1 has not been initialized yet 
   !      (at that time one still needs to load the bathymetry information) .....
   IF (iopss > 0) THEN

     ALLOCATE(Qpss(km1,iopss), & 
	          Tpss(km1,iopss), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 21 ) 	 
     IF (ntr > 0) THEN
       ALLOCATE(Rpss(km1,iopss,ntr), STAT=istat)
       IF (istat /= 0) CALL allocate_error ( istat, 22 )
       ALLOCATE(oteffnn(npssdev,iopss), STAT=istat)
       IF (istat /= 0) CALL allocate_error ( istat, 22 )
     ENDIF

   ENDIF 

   ! .... Allocate kt array
   ALLOCATE(kt(im1,jm1), STAT=istat)
   IF (istat /= 0 ) CALL allocate_error ( istat, 23)

   ! .... Allocate epskt array
   ALLOCATE(epskt(im1,jm1), STAT=istat)
   IF (istat /= 0 ) CALL allocate_error ( istat, 24)

   ! .... Allocate k4sodout array
   ALLOCATE(k4sodout(im1,jm1), STAT=istat)
   IF (istat /= 0 ) CALL allocate_error ( istat, 25)

   ! .... Allocate tracerbot array
   ALLOCATE(tracerbot(im1,jm1), STAT=istat)
   IF (istat /= 0 ) CALL allocate_error ( istat, 26)

   ! .... Allocate Hnet array
   ALLOCATE(Hnet(im1,jm1), STAT=istat)
   IF (istat /= 0 ) CALL allocate_error ( istat, 27)

END SUBROUTINE AllocateSpace2

!************************************************************************
SUBROUTINE bathy
!************************************************************************
!
!  Purpose: To read the file with bathymetry for the basin. The depths
!           are assumed to be defined at the corners of each computational
!           cell and stored in the array 'h4'. The average depths along
!           each cell face are computed and used to define the 3-D layer
!           thickness array hhs of bottom depths from datum at zeta points. 
!           2-D logical mask arrays are defined with .TRUE. values for all 
!           wet cells and .FALSE. values for all dry cells. The 2-D arrays, 
!           hhu & hhv, of depths at u- & v-points are defined from hhs. 
!           The depths are in meters below a datum.
!  13/11/08 Input depths at zeta point
!
!-------------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER(LEN=50) :: bathymetry_file = "h"
   INTEGER :: i, j, k, l, ios, istat, kb, is, ie, js, je
   INTEGER :: imm, jmm, ncols, ncols1, nc, nn, ia, ib
   REAL :: hs1, udepth, vdepth, ztop
   CHARACTER(LEN=13) :: fmt
   REAL, DIMENSION(:,:), POINTER :: h4

   !.....Open bathymetry file.....
   OPEN (UNIT=i5, FILE=bathymetry_file, STATUS="OLD", IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//bathymetry_file, ios )

   !.....Read header information.....
   READ (UNIT=i5, FMT='(37X,I5,6X,I5,8X,I4//)', IOSTAT=ios) imm, jmm, ncols
   IF (ios /= 0) CALL input_error ( ios, 11 )

   !.....Check grid dimensions against input parameters.....
   IF ((im /= imm+1) .OR. (jm /= jmm+1)) THEN
      PRINT *, " ****ERROR -- Grid size computed from input file does not"
      PRINT *, "              agree with the header in the bathymetry file"
      PRINT '(4(A,I5))', " im=", im, " imm+1=", imm+1, " jm=", jm, " jmm+1=", jmm+1
      PRINT '(A/)', " "
      PRINT *, " ****STOPPING si3d in SUBROUTINE bathy"
      STOP
   END IF

   !.....Write data format to an internal file.....
   IF (ibathyf == 1) THEN ! Stockton
     WRITE (UNIT=fmt, FMT='("(5X,", I4, "G4.0)")') ncols ! idt real
   ELSE                   !General case (Deeper lakes > 100 m)
     WRITE (UNIT=fmt, FMT='("(5X,", I4, "G5.0)")') ncols ! idt real
   ENDIF

   !.....Allocate space for bathymetry array.....
   ALLOCATE ( h4(1:im1, 1:jm1), STAT=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 15 )

   !.....Read bathymetry.....
   ncols1 = ncols     - 1; 
   nc     = imm/ncols + 1; 
   ia     = -ncols1   + 1; 
   DO nn = 1, nc
      ia = ia + ncols
      ib = ia + ncols1
      IF ( ia > im ) EXIT
      IF ( ib > im ) ib = im
      DO j = jm, j1, -1 ! change 1 to j1
         READ (UNIT=i5, FMT=fmt, IOSTAT=ios) (h4(i,j), i = ia, ib)
         IF (ios /= 0) CALL input_error ( ios, 12 )
      END DO
   END DO

   !.....Close bathymetry file.....
   CLOSE (UNIT=i5)

   !.....Convert data (in dm to meters) & 
   !     add adjustment to bathymetry datum ..... 
   h4 = h4 * 0.1 + datadj ! idt real
 
   !.....Add fictitious row/column of depths around grid.....
   h4(1,j1:jm  ) = h4(2,j1:jm )         ! west side
   h4(i1:im,jm1) = h4(i1:im,jm)         ! north side
   h4(im1,j1:jm) = h4(im,j1:jm)         ! east side
   h4(i1:im,1  ) = h4(i1:im,2 )         ! south side
   ! Take care of corners
   h4(  1,  1) = h4( 2, 2); h4(  1,jm1) = h4( 2,jm) 
   h4(im1,jm1) = h4(im,jm); h4(im1,  1) = h4(im, 2)

   !.....Define mask2d array.....
   DO j = 1, jm1; DO i = 1, im1
     IF ( h4(i,j) > 0.0 ) THEN 
        mask2d(i,j) = .TRUE.
     ELSE
        mask2d(i,j) = .FALSE.
     END IF
   END DO; END DO

   !.....Be sure mask2d is false in the fictitious
   !     row/column around the outer edge of the grid.....
   mask2d(  1,1:jm1) = .FALSE.; mask2d(1:im1,jm1) = .FALSE.
   mask2d(im1,1:jm1) = .FALSE.; mask2d(1:im1,  1) = .FALSE.

   !.....Define layer No. for bottom cell (kmz) & 
   !     bottom depth from datum at zeta-points (hhs).....
   hhs = ZERO; 
   DO j = 1,jm1; DO i = 1,im1

     SELECT CASE (mask2d(i,j))

     CASE (.FALSE.)      ! Cells with all dry layers

       kmz(i,j) = km1

     CASE (.TRUE.)      ! Cells with wet layers

       ! ... Take depth at zeta-point as given in h file
       hs1 = h4(i,j)

       ! ... Compute No. of wet layers & depths at zeta-points
       ! --- A. constant layer thickness (ibathyf >= 0) 
       IF (ibathyf >= 0) THEN
         kmz(i,j) = FLOOR(hs1/ddz)
         IF( (hs1-kmz(i,j)*ddz) > dzmin ) THEN
           kmz(i,j) = kmz(i,j) + 1;
           hhs(i,j) = hs1
         ELSE
           hhs(i,j) = kmz(i,j)*ddz
         ENDIF
         ! Add the fictitious first layer above the water surface
         kmz(i,j) = kmz(i,j) + 1

       ! --- B. Variable layer thickness (ibathyf < 0)
       ELSE 
         DO k = k1, km         
           IF (zlevel(k+1)>=hs1) THEN
             ! ... Option 1 - it works
             hhs(i,j) = zlevel(k)+MAX(dzmin,(hs1-zlevel(k)));
             kmz(i,j) = k 
             EXIT
             !! ... Option 2 - preferable from a theoretical stand point
             !IF ( hs1-zlevel(k) > dzmin) THEN
             !  hhs(i,j) = hs1
             !  kmz(i,j) = k 
             !ELSE
             !  hhs(i,j) = zlevel(k)
             !  kmz(i,j) = k - 1
             !  IF ( kmz(i,j) == 1 ) THEN
             !     hhs(i,j) = dzmin
             !     kmz(i,j) = k1
             !  ENDIF
             !ENDIF
             !EXIT
           ENDIF
         ENDDO
       ENDIF

     END SELECT

   END DO; END DO

   !.....Define bottom depths from datum at u- and v- points
   hhu = ZERO; 
   hhv = ZERO;
   DO j = 1,jm; DO i = 1,im;   
      IF(mask2d(i+1,j)) THEN         
        hhu(i,j) = MIN(hhs(i+1,j), hhs(i,j))
      ENDIF 
      IF(mask2d(i,j+1)) THEN         
        hhv(i,j) = MIN(hhs(i,j+1), hhs(i,j))
      ENDIF 
   END DO; END DO

   !.....Process and output the bathymetry needed for
   !     graphics and particle tracking if ioutg=1.....
   IF ( ipxml > 0 ) CALL outg ( h4 )

   !.....Compute the first and last column and row of grid 
   !     with wet points. Variables used in subr. solver
   ifirst = im; ilast = i1; jfirst = jm; jlast = j1
   DO j = j1, jm; DO i = i1, im
     IF ( .NOT. mask2d(i,j) ) CYCLE   ! Ignore dry points
     IF ( i < ifirst ) ifirst = i
     IF ( i > ilast  ) ilast  = i
     IF ( j < jfirst ) jfirst = j
     IF ( j > jlast  ) jlast  = j
   END DO; END DO

   ! ... Find dimensions for 2D-lk arrays 
   l = 0
   DO i = i1, im; DO j = j1, jm; 
     IF(.NOT. mask2d(i,j)) CYCLE
     l = l + 1
   END DO; END DO;
   lm = l ; lm1 = lm + 1;

   ! ... Allocate space for arrays in 2D-lk coordinates
   CALL AllocateSpace2

   ! ... Mapping functions from/to 3D-(i,j)- to/from 2D-l-indexes
   l = 0
   ij2l = lm1; ! For dry(i,j) the map function ij2l will yield lm1
   DO i = i1, im; DO j = j1, jm; 
     IF(.NOT. mask2d(i,j)) CYCLE
     l = l + 1
     l2i (l  ) = i ; ! Goes from ipl to index i in the i,j plane
     l2j (l  ) = j ; ! Goes from ipl to index j in the i,j plane
     ij2l(i,j) = l ; ! Goes from i,j to the ipl index in the ipl line
   END DO; END DO;

   ! ... Assign E,W,N,S colums for each l-column in the 2D-l space
   l = 0
   DO i = i1, im; DO j = j1, jm; 
     IF(.NOT. mask2d(i,j)) CYCLE
     l = l + 1
     lEC(l) = ij2l(i+1,j); ! Defines water column East  of l
     lWC(l) = ij2l(i-1,j); ! Defines water column West  of l
     lNC(l) = ij2l(i,j+1); ! Defines water column North of l 
     lSC(l) = ij2l(i,j-1); ! Defines water column South of l
   END DO; END DO;

   !.....Deallocate h4 pointer array.....
   DEALLOCATE ( h4 )

END SUBROUTINE bathy

!***********************************************************************
SUBROUTINE fd
!***********************************************************************
!
!  Purpose: fd  is the supervisory subroutine for the finite
!           difference method. It advances the solution one time
!           step. Vertical diffusion is treated implicitly.
!
!-----------------------------------------------------------------------

   ! ... Local variables
   INTEGER :: itr,itr2

   ! ... Define tz used in all following routines ...............
   tz = 1.0/istep


   IF (idbg == 1) PRINT *, " readbcNGB"

   ! ... Read in nested boundary conditions, if specified .......
   CALL readbcNGB

   IF (idbg == 1) PRINT *, " openbcUVH"

   !.....Assign new values of s or u/v along open boundaries.....
   CALL openbcUVH

   IF (idbg == 1)  PRINT *, " PointSourceSinkSolver"

   !.... Solve Plume Model (takes rhop & estimates Qpss & Tpss) - iopss
   CALL PointSourceSinkSolver

   IF (idbg == 1) PRINT *, " surfbc"

   ! ... Assign boundary conditions at free surface .............
   CALL surfbc 

   IF (idbg == 1) PRINT *, " exmom(1)"

   !.....Evaluate explicit terms in x-momentum eq................
   CALL exmom(1)

   IF (idbg == 1) PRINT *, " matmom(1)"

   !.....Solve a tridiagonal system at each horizontal node
   !     to obtain the matrices for the x-momentum equations.....
   CALL matmom(1)

   IF (idbg == 1) PRINT *, " exmom(2)"

   !.....Evaluate explicit terms in y-momentum eq................
   CALL exmom(2)

   IF (idbg == 1) PRINT *, " matmom(2)"

   !.....Solve a tridiagonal system at each horizontal node
   !     to obtain the matrices for the y-momentum equations.....
   CALL matmom(2)

   IF (idbg == 1) PRINT *, " matcon"

   !.....Calculate matrix coefficients for implicit free surface
   CALL matcon

   IF (idbg == 1) PRINT *, "solverSPARSE"

   !.....Solve implicit system of equations for zeta............
   CALL solverSPARSE

   IF (idbg == 1) PRINT *, " openbcUVH"

   !.....Reassign new values of s or u/v along open boundaries.....
   CALL openbcUVH

   IF (idbg == 1) PRINT *, " vel"

   !.....Solve for velocities explicitly. If DRYING occurs 
   !     at this point, dry cells are removed and k1z/k1u/k1v
   !     recalculated (wetting is done after finishing the 
   !     calculations for a given time step or iteration) ......
   CALL vel

   IF (idbg == 1) PRINT *, " Scalar Transport"

   ! ... Solve for active scalar transport
   IF (isal /= 0) THEN
     CALL exsal 
     CALL imsal
     CALL openbcSCA
   END IF

   IF (idbg == 1) PRINT *, " Tracers, ntr, niter, lastiter ",ntr,niter,lastiter

   !.....Solve for non-active scalar transport 
   IF (ntr      >  0 .AND. & 
       niter    >  0 .AND. &
       lastiter == 1      ) THEN

       IF      (ecomod <  0) THEN
         CALL srcsnk00         
       ELSE IF (ecomod == 0) THEN 
         CALL srcsnk00 
       ELSE IF (ecomod == 1) THEN
         CALL srcsnkWQ
       ELSE IF (ecomod == 2) THEN
         CALL srcsnkSZ
       ELSE IF (ecomod == 3) THEN
         CALL srcsnkSD
       ELSE IF (ecomod == 4) THEN
         CALL do_fabm_aed
       ENDIF

       DO itr = 1, ntr
         IF (ecomod < 0 .AND. ( trct0(itr) > n .OR. trctn(itr) < n ) ) CYCLE
         CALL exTracer (itr)
         CALL ImTracer (itr)
       ENDDO

!       ! ... openbc tracers are updated when they are specified. - Chris DISABLED
!       IF (idbg == 1) PRINT *, " entering openbcTracer assignment rdtr = ", rdtr
!       DO itr2= 1, rdtr
!         CALL openbctracer (itr2)
!       ENDDO

        IF (idbg == 1) PRINT *, " determining tracer injctn tinj & its ", tinj,its
        IF (iinj > 0 .AND. isinjected == 0) THEN
          IF (tinj <= its) THEN
            CALL do_tracer_injection()
          ENDIF
        ENDIF

   ENDIF

   IF (idbg == 1) PRINT *, " vel2"

   ! ... Recalculate near surface velocity to account 
   !     for WETTING & recalculate k1z, k1u & k1v ..............
   CALL vel2

   IF (idbg == 1) PRINT *, " smooth"

   !.....Smooth solution on leapfrog step if ismooth>=1.........
   IF (ismooth >= 1 .AND. istep == 1) CALL smooth

   IF (idbg == 1) PRINT *, " Mixing Coefficients"

   !.....Assing eddy viscosity and diffusivity at n+1 ...........
   CALL UpdateMixingCoefficients

END SUBROUTINE fd

!***********************************************************************
SUBROUTINE exmom ( ieq )
!***********************************************************************
!
!  Purpose: To evaluate the explicit terms (advection, Coriolis, and
!           horizontal diffusion) in the momentum equations. The sum
!           of these terms are saved in the 3-D array  ex(i,j,k)
!           which is the primary output from this subroutine. The
!           horizontal advection terms can be evaluated by either upwind
!           or centered differencing.
!           (note: ex(i,j,k)  is a temporary workspace array that is
!           used for storing the explicit terms from both the x- and y-
!           momentum eqs and is also used later for storing the explicit
!           terms from the salinity eq in SUBR. exsal.)
!
!  Dummy argument:
!  ieq    = Parameter indicating whether the explicit terms in the
!           X  or  Y  momentum equation are to be evaluated.
!           (1=x-momentum, 2=y-momentum)
!
!-----------------------------------------------------------------------
   !.....Argument.....
   INTEGER, INTENT(IN) :: ieq

   !.....Local variables.....
   REAL :: twodt1
   REAL :: corx, cory, advx, advy, hdx, hdy, uE, uW, vN, vS, wU, wD, &
           scW, scE, scN, scS, scU, scD, gxh, gyh
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne
   INTEGER :: kmne, nwlayers, nn, is, ie, js, je

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)
   n_exmom = n_exmom + 1

   !.....Constant.....
   twodt1 = twodt*tz

   ex = 0.0 ! FJRGOTM

   !.....Choose to evaluate explicit terms for either the x- or y-mom eq.....
   SELECT CASE (ieq)

   ! -----X-momentum equation-----
   CASE (1)
        
      !......Calulate the explicit terms by sweeping over interior u-pts &
      !      store results in ex- array                 

      !.....Calculate the coefficient array haypp for use
      !     in horizontal diffusion term & th1,th for use
      !     in vertical advection term in the x-momentum
      haxpp(:,lm1) = 0.0; 
      haypp(:,lm1) = 0.0;
      DO l = 1, lm

        ! ... Map 3D-(i,j) from 2D-l indexes
        i = l2i(l); j = l2j(l);

        ! ... Cycle if E-column is dry
        IF ( .NOT. mask2d(i+1,j) ) CYCLE

        !.....Compute layer number for the bottom wet u-pt.......
        kmx = MIN(kmz(i+1,j), kmz(i,j))
        k1x =                 k1u(i,j)

        ! ... Horizontal diffusion ...........................
        DO k = k1x,kmx
          haxpp(k,l)=Ax0*MIN(hupp(k,lEC(l)),hupp(k,l)) 
          haypp(k,l)=Ay0*MIN(hupp(k,lNC(l)),hupp(k,l))
        ENDDO

        SELECT CASE (itrmom)

        CASE(1)

          !.....Calculate weighting arrays for vertical advection 
          DO k = k1x, kmx
            th (k,l) = hup(k-1,l)/(hup(k-1,l)+hup(k,l))
            th1(k,l) = 1.-th(k,l)
          ENDDO

          !.....Set th=th1 at the free surface & bottom 
          th (k1x  ,l) = 0.0; 
          th1(k1x  ,l) = 0.0;
          th (kmx+1,l) = 0.5; 
          th1(kmx+1,l) = 0.5;

        CASE DEFAULT

        END SELECT

      END DO


      !......Calulate the explicit terms by sweeping over interior u-pts.....
      DO l = 1, lm
 
         ! ... Map 2D-l into 3D-(i,j) indexes
         i = l2i(l); j = l2j(l);

         ! ... Cycle if E-column is dry
         IF ( .NOT. mask2d(i+1,j) ) CYCLE

         ! Compute the layer number for the bottom wet u-pt
         kmx = MIN(kmz(i+1,j), kmz(i,j))
         k1x =                 k1u(i,j)

         ! Compute explicit term	     
         DO k = k1x,kmx
                                       
            ! ... For u-layers connecting wett & dry cells neglect
            !     contribution from advective, coriolis & diffusion
            IF ( hp(k,l) <= ZERO .OR. hp(k,lEC(l)) <= ZERO) THEN 
              ex(k,l) = uhpp(k,l)            
              CYCLE
            ENDIF

            !.....Coriolis.....
            corx = 0.25 * f * (vhp(k,     lEC(l) ) + vhp(k,    l )       &
                  &           +vhp(k, lSC(lEC(l))) + vhp(k,lSC(l)))
                      
            !.....Advection   
            uE = uhp(k,    lEC(l) ) +uhp(k,    l )
            uW = uhp(k,        l  ) +uhp(k,lWC(l))
            vN = vhp(k,    lEC(l) ) +vhp(k,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k,    l ); IF ( k == k1x ) wU = 0.0; 
            wD = wp (k+1,  lEC(l) ) +wp (k+1  ,l ); IF ( k == kmx ) wD = 0.0;  

            SELECT CASE (itrmom)

            CASE (1)   ! Centered differences with th & th1 

              scE = up(k,lEC(l))+up(k,    l ) 
              scW = up(k,lWC(l))+up(k,    l ) 
              scN = up(k,lNC(l))+up(k,    l ) 
              scS = up(k,lSC(l))+up(k,    l )
              advx = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy 
              advx=advx+(wU*(th (k  ,l)* up(k  ,l)  +          &
                             th1(k  ,l)* up(k-1,l)) -          &
                         wD*(th (k+1,l)* up(k+1,l)  +          &
                             th1(k+1,l)* up(k  ,l)) ) / 2.

            CASE (2) ! Upwinding all

              advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                       (uE-ABS(uE))* upp(k,lEC(l)) -           &
                       (uW+ABS(uW))* upp(k,lWC(l)) -           &
                       (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                    +( (vN+ABS(vN))* upp(k,    l ) +           &
                       (vN-ABS(vN))* upp(k,lNC(l)) -           &
                       (vS+ABS(vS))* upp(k,lSC(l)) -           &
                       (vS-ABS(vS))* upp(k,    l ) ) / fourdy  &  
                    +( (wU+ABS(wU)) * upp(k  ,l) +             &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.         & 
                     -((wD+ABS(wD)) * upp(k+1,l) +             &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.

            CASE (3)  ! Centered differences - avoid computation of th1 and th

              scE = up(k,lEC(l))+up(k,    l ) 
              scW = up(k,lWC(l))+up(k,    l ) 
              scN = up(k,lNC(l))+up(k,    l ) 
              scS = up(k,lSC(l))+up(k,    l )
              scU = (up(k  ,l)*hup(k  ,l)+        &
                     up(k-1,l)*hup(k-1,l))/       &
                    (hup(k ,l)+hup(k-1,l)) 
              scD = (up(k  ,l)*hup(k  ,l)+        &
                     up(k+1,l)*hup(k+1,l))/       &
                    (hup(k ,l)+hup(k+1,l)) 
              advx = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy +          &
                     (wU * scU - wD * scD ) / 2.

            CASE (4) ! Upwinding for horizontal & centered for vertical  

              advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                       (uE-ABS(uE))* upp(k,lEC(l)) -           &
                       (uW+ABS(uW))* upp(k,lWC(l)) -           &
                       (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                    +( (vN+ABS(vN))* upp(k,    l ) +           &
                       (vN-ABS(vN))* upp(k,lNC(l)) -           &
                       (vS+ABS(vS))* upp(k,lSC(l)) -           &
                       (vS-ABS(vS))* upp(k,    l ) ) / fourdy    
              advx=advx+(wU*(th (k  ,l)* up(k  ,l)  +          &
                             th1(k  ,l)* up(k-1,l)) -          &
                         wD*(th (k+1,l)* up(k+1,l)  +          &
                             th1(k+1,l)* up(k  ,l)) ) / 2.

            END SELECT
	      		                  
            !.....Horizontal diffusion.....
            IF ( ihd == 0) THEN
              hdx = 0.0E0
            ELSEIF ( ihd == 1) THEN
               hdx  =   (haypp(k,    l )*(upp(k,lNC(l))-upp(k,    l ))       &
                    & -  haypp(k,lSC(l))*(upp(k,    l )-upp(k,lSC(l))))/dydy &
                    & + (haxpp(k,    l )*(upp(k,lEC(l))-upp(k,    l ))       &
                    & -  haxpp(k,lWC(l))*(upp(k,    l )-upp(k,lWC(l))))/dxdx
            ELSEIF ( ihd >  1) THEN
              hdx  = 2.* (haxpp(k,   (l))*(upp(k,lEC(    l ))-upp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(upp(k,        l  )-upp(k,lWC(l))))/dxdx & 
                 &   +   (haypp(k,    l )*(upp(k,lNC(    l ))-upp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(upp(k,        l  )-upp(k,lSC(l))))/dydy & 
                 &   +   (haypp(k,    l )*(vpp(k,lEC(    l ))-vpp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(vpp(k,lEC(lSC(l)))-vpp(k,lSC(l))))/dxdy
            ENDIF

            !! ... Adjust terms to account for boundary conditions (North Delta SALMON)
            !IF((                           j>=  jm1 - 20          )   .OR. &
            !   (i<=591               .AND. j<=  60                )   .OR. &
            !   (i<=20                                             )   .OR. &
            !   (i>=990 .AND. i<=1010 .AND. j>=  74 .AND. j <=  81 )   .OR. &
            !   (                           j<=  20                )   .OR. &
            !   (i>=678 .AND. i<= 682 .AND. j>= 100 .AND. j <= 115 )) THEN
            !   corx = 0.0
            !   advx = 0.0
            !   hdx  = 2.*hdx
            !ENDIF
            !IF (j < 20 .OR. j > jm1-20) THEN; hdx = 0.0; corx = 0.0; end IF; 

            ! ... Gravity along x-direction (open channel flows)
            gxh = g * xslope * hup(k,l)

            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv-corx-hdx-gxh)
                
         END DO
      END DO

      ! ... Recalculate ex for near bdry. cells (Comment out for North Delta Salmon)
      IF (nopen > 0) CALL MODexmom4openBCX        
      
   ! -----Y-momentum equation-----      
   CASE (2)

      !.....Calculate the explicit terms by sweeping over interior v-pts &
      !     store results in ex-array

      !.....Calculate the coefficient array haxpp for use
      !     in horizontal diffusion term & th1,th for use
      !     in vertical advection term in the y-momentum 
      haxpp(:,lm1) = 0.0; 
      haypp(:,lm1) = 0.0;
      DO l = 1, lm

         ! ... Map 3D-(i,j) from 2D-l indexes
         i = l2i(l); j = l2j(l);

         ! ... Cycle if N-column is dry
         IF ( .NOT. mask2d(i,j+1) ) CYCLE

         !.....Compute layer number for top & bottom wet v-pt.....
         kmy = MIN(kmz(i,j+1), kmz(i,j))
         k1y =                 k1v(i,j)

         ! ... Horizontal diffusion
         DO k = k1y, kmy
           haypp(k,l)=Ay0*MIN(hvpp(k,lNC(l)),hvpp(k,l)) 
           haxpp(k,l)=Ax0*MIN(hvpp(k,lEC(l)),hvpp(k,l))
         ENDDO

         SELECT CASE (itrmom)

         CASE (1)

           !.....Calculate weighting arrays for vertical advection 
           DO k = k1y, kmy
              th (k,l) = hvp(k-1,l)/(hvp(k-1,l)+hvp(k,l))
              th1(k,l) = 1.-th(k,l)
           ENDDO
             
           !.....Set th=th1 at the free surface & bottom 
           th (k1y  ,l) = 0.0; 
           th1(k1y  ,l) = 0.0;
           th (kmy+1,l) = 0.5; 
           th1(kmy+1,l) = 0.5;              
         
         CASE DEFAULT

         END SELECT

      END DO
      !......Calulate the explicit terms by sweeping over interior u-pts.....
      DO l = 1, lm

         ! ... Map 2D-l into 3D-(i,j) indexes
         i = l2i(l); j = l2j(l);
         
         ! ... Cycle if N-column is dry
         IF ( .NOT. mask2d(i,j+1) ) CYCLE

         ! Compute the layer number for top & bottom wet v-pt
         kmy = MIN(kmz(i,j+1), kmz(i,j))
         k1y =                 k1v(i,j)

         ! Compute explicit term	     
         DO k = k1y,kmy
                                    
            ! ... For v-layers connecting wett & dry cells neglect
            !     contribution from advective, coriolis & diffusion
            IF ( hp(k,l) <= ZERO .OR. hp(k,lNC(l)) <= ZERO) THEN 
               ex(k,l) = vhpp(k,l)
               CYCLE
            ENDIF
         
            !.....Coriolis.....
            cory = 0.25 * f * (uhp(k,     lNC(l) ) + uhp(k,    l )       &
                  &           +uhp(k, lWC(lNC(l))) + uhp(k,lWC(l)))
                      
            !.....Advection  
            uE = uhp(k,    lNC(l) ) +uhp(k,    l );
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l));
            vN = vhp(k,    lNC(l) ) +vhp(k,    l );
            vS = vhp(k,        l  ) +vhp(k,lSC(l));
            wU = wp (k,    lNC(l) ) +wp (k,    l ); IF ( k == k1y ) wU = 0.0;
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l ); IF ( k == kmy ) wD = 0.0;

            SELECT CASE ( itrmom)

            CASE (1)  ! Centered differences using th & th1 factors

              scE = vp(k,lEC(l)) + vp(k,    l )
              scW = vp(k,lWC(l)) + vp(k,    l ) 
              scN = vp(k,lNC(l)) + vp(k,    l ) 
              scS = vp(k,lSC(l)) + vp(k,    l ) 
              advy = (uE * scE - uW * scW ) / fourdx +       &
                     (vN * scN - vS * scS ) / fourdy   
              advy = advy  +                                 &
                     (wU*(th (k  ,l)* vp(k  ,l)  +           &
                          th1(k  ,l)* vp(k-1,l)) -           &
                      wD*(th (k+1,l)* vp(k+1,l)  +           &
                          th1(k+1,l)* vp(k  ,l)) ) / 2.

            CASE(2) ! Upwinding
			      
              advy = ( (uE+ABS(uE)) * vpp(k,    l ) +          &
                  &    (uE-ABS(uE)) * vpp(k,lEC(l)) -          &
                  &    (uW+ABS(uW)) * vpp(k,lWC(l)) -          &
                  &    (uW-ABS(uW)) * vpp(k,    l ) ) / fourdx &
                  & +( (vN+ABS(vN)) * vpp(k,    l ) +          &
                  &    (vN-ABS(vN)) * vpp(k,lNC(l)) -          &
                  &    (vS+ABS(vS)) * vpp(k,lSC(l)) -          &
                  &    (vS-ABS(vS)) * vpp(k,    l ) ) / fourdy &
                  & +( (wU+ABS(wU)) * vpp(k  ,l) +             &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.         &
                    -( (wD+ABS(wD)) * vpp(k+1,l) +             &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.

            CASE (3)  ! Centered differences - avoid computation of th1 and th

              scE = vp(k,lEC(l)) + vp(k,    l )
              scW = vp(k,lWC(l)) + vp(k,    l ) 
              scN = vp(k,lNC(l)) + vp(k,    l ) 
              scS = vp(k,lSC(l)) + vp(k,    l ) 
              scU = (vp(k  ,l)*hvp(k  ,l)+                      &
                     vp(k-1,l)*hvp(k-1,l))/                     &
                    (hvp(k ,l)+hvp(k-1,l)) 
              scD = (vp(k  ,l)*hvp(k  ,l)+                      &
                     vp(k+1,l)*hvp(k+1,l))/                     &
                    (hvp(k ,l)+hvp(k+1,l)) 
              advy = (uE * scE - uW * scW ) / fourdx +          &
                     (vN * scN - vS * scS ) / fourdy +          &
                     (wU * scU - wD * scD ) / 2.

            CASE(4) ! Upwinding only for horizontal advection
			      
              advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                  &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                  &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                  &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                  & +( (vN+ABS(vN))* vpp(k,    l ) +           &
                  &    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                  &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                  &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  
              advy = advy +                                    &
                     (wU*(th (k  ,l)* vp(k  ,l)  +             &
                          th1(k  ,l)* vp(k-1,l)) -             &
                      wD*(th (k+1,l)* vp(k+1,l)  +             &
                          th1(k+1,l)* vp(k  ,l)) ) / 2.
		
            END SELECT

            !.....Horizontal diffusion.....
            IF ( ihd == 0) THEN
              hdy = 0.0E0
            ELSEIF ( ihd == 1) THEN  		                  
              hdy  =   (haypp(k,    l )*(vpp(k,lNC(l))-vpp(k,    l ))              &
                   & -  haypp(k,lSC(l))*(vpp(k,    l )-vpp(k,lSC(l))))/dydy        &
                   & + (haxpp(k,    l )*(vpp(k,lEC(l))-vpp(k,    l ))              &
                   & -  haxpp(k,lWC(l))*(vpp(k,    l )-vpp(k,lWC(l))))/dxdx
            ELSEIF ( ihd >  1) THEN
              hdy  = 2.* (haypp(k,   (l))*(vpp(k,lNC(    l ))-vpp(k,    l ))       &
                 &     -  haypp(k,lSC(l))*(vpp(k,        l  )-vpp(k,lSC(l))))/dydy & 
                 &   +   (haxpp(k,    l )*(vpp(k,lEC(l     ))-vpp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(vpp(k,        l  )-vpp(k,lWC(l))))/dxdx & 
                 &   +   (haxpp(k,    l )*(upp(k,lNC(    l ))-upp(k,    l ))       &
                 &     -  haxpp(k,lWC(l))*(upp(k,lNC(lWC(l)))-upp(k,lWC(l))))/dxdy
            ENDIF

            ! ... Adjust terms to account for boundary conditions (North Delta SALMON)
            !IF((                           j>=  jm1 - 20          )   .OR. &
            !   (i<=591               .AND. j<=  60                )   .OR. &
            !   (i<=20                                             )   .OR. &
            !   (i>=990 .AND. i<=1010 .AND. j>=  74 .AND. j <=  81 )   .OR. &
            !   (                           j<=  20                )   .OR. &
            !   (i>=678 .AND. i<= 682 .AND. j>= 100 .AND. j <= 115 )) THEN
            !   cory = 0.0
            !   advy = 0.0
            !   hdy  = 2.*hdy
            !ENDIF
            !IF (j < 20 .OR. j > jm1-20) THEN; hdy = 0.; cory = 0.0; end IF
               
            ! ... Gravity along y-direction (open channel flows)
            gyh = g * yslope * hvp(k,l)
               
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv+cory-hdy-gyh)
               
         END DO
      END DO

      ! ... Recalculate ex for near bdry. cells (Comment out for North Delta Salmon)
      IF (nopen > 0) CALL MODexmom4openBCY    
    
   END SELECT

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exmom = t_exmom + (etime - btime)

END SUBROUTINE exmom

!***********************************************************************
SUBROUTINE matmom ( ieq )
    !***********************************************************************
    !
    !  Purpose:   To define the matrices for the momentum equations.
    !
    !  Algorithm: The x-momentum equations at each vertical array of
    !             u-pts are first expressed in the compact matrix form
    !
    !               [aa] [uh] = [gg] - g*dt/dx*rho*(s(i+1,j)-s(i,j))*[hh]
    !
    !             by defining the three matrices  [hh],  [gg], and  [aa].
    !             (Because the  [aa]  matrix is tridiagonal, only the
    !             diagonals are stored.) The above system of equations is
    !             then rearranged into the form
    !
    !               [uh] = [ag] - g*dt/dx*rho*(s(i+1,j)-s(i,j))*[ar]
    !
    !             by indirect solution using the tridiagonal solver  trid.
    !             The matrices  [ag]  and  [ar]  are the output from the  trid
    !             subroutine along with their summation over the depth at
    !             each horizontal node point,  [eag]  and [ear]. The matrices
    !             for the x-momentum eq are stored in the fortran arrays
    !             agx,  arx,  eagx, and  earx. Everything is similar for
    !             the y-momentum eq.
    !
    !  Dummy argument:
    !  ieq    = Parameter indicating whether the matrices for the
    !           x  or  y  momentum equation are to be evaluated.
    !           (1=x-momentum, 2=y-momentum)
    !
    !  23/04/2008   F.J. Rueda    Recompute baroclinic term at bottom
    !  23/04/2008   F.J. Rueda    Do not set to zero baroclinic term at
    !                             wett/dry cells near the surface
    !  26/05/2009   F.J. Rueda    Include ptype option 3 (circular plume)
    !-----------------------------------------------------------------------

    !.....Argument.....
    INTEGER, INTENT(IN) :: ieq

    !.....Local variables.....
    INTEGER :: i, j, k, l, kmx, kmy, k1x, k1y, nwlayers, inn
    REAL    :: twodt1, wsx0, wsy0, tausx, taubx, tausy, tauby,      &
        hpmin, usurf, vsurf,Uhvalue, Usource, Vsource,       &
        Vhvalue, cwx, cwy
    REAL, DIMENSION(km1) :: bclncx, bclncy, vdiffx, vdiffy,         &
        & hupdrho, hvpdrho, deltaz, Avxdudz, Avydvdz, Avx, Avy,   &
        & rhopx, rhopy, gg, hh, ar, ag
    REAL, DIMENSION(3,km1) :: aa
    REAL, PARAMETER :: rhop0 = 1000.

    !.....Timing.....
    REAL, EXTERNAL :: TIMER
    REAL :: btime, etime
    btime = TIMER(0.0)

    !.....Constants.....
    twodt1 = twodt*tz;

    SELECT CASE ( ieq )
   
        !                -----X-momentum equation-----
        CASE (1)
            !eagx = 0.0; earx = 0.0;
            !.....Loop over interior u-pts .....
            DO l = 1, lm
            
                ! ... Map 2D-l into 3D-(i,j) indexes
                i = l2i(l); j = l2j(l);

                ! ... Skip for dry u-pts
                IF (.NOT.mask2d(i+1,j)) CYCLE

                ! ... Compute the layer number for top & bottom wet u-pt
                kmx = MIN(kmz(i+1,j), kmz(i,j))
                k1x =                 k1u(i,j)
                nwlayers = (kmx-k1x) + 1
                IF(nwlayers < 1) CYCLE

                ! ... Compute eddy viscosity at interfaces vertically between u-pts
                Avx = 0.0
                DO k = k1x+1,kmx
                    Avx(k) = 0.5 * (Av(k,lEC(l))+Av(k,l))
                ENDDO

                ! ... Define average layer density at u-pt (in kg/m**3) ...........
                rhopx(k1x:kmx) = 1000. ! Neglect vertical density variations

                ! ... Compute explicit portion of water surface slope term ........
                wsx0 = rhopx(k1x) * gdtdx * (spp(i+1,j) - spp(i,j))

                ! ... Compute baroclinic term .....................................
                SELECT CASE (ibc)
                    CASE (0)     ! No baroclinic term
                        bclncx(k1x:kmx) = 0.0
                    CASE (1:)
                        DO k = k1x, kmx
                            hupdrho(k) = gdtdx*hup(k,l)*(rhop(k,lEC(l))-rhop(k,l))
                          !IF(hp(k,l) <= ZERO .OR. hp(k,lEC(l)) <= ZERO) hupdrho(k) = 0.0
                        ENDDO
                        bclncx(k1x) = hupdrho(k1x)
                        IF (kmx > k1x) THEN   ! Two or more wet layers
                            ! Recompute bottom layer baroclinic term along horizontal plane
                            CALL bclnc_km (l, kmx, 1, hupdrho(kmx) )
                            DO k = k1x+1, kmx
                                bclncx(k) = bclncx(k-1) + hupdrho(k-1) + hupdrho(k)
                            END DO
                        END IF
                END SELECT
                ! ... Compute explicit portion of vertical diffusion term ........
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        vdiffx(k1x) = 0.0
                    CASE (2:)   ! Two or more wet layers (hupp->hup)
                        DO k = k1x+1,kmx
                            deltaz(k) = hup(k-1,l) + hup(k,l)
                            Avxdudz(k)= Avx(k)*(upp(k-1,l)-upp(k,l))/deltaz(k)
                        ENDDO
                        Avxdudz(k1x)      = 0.0  ! Set value at free surface to zero
                        Avxdudz(kmx+1)    = 0.0  ! Set value at bottom boundary to zero
                        vdiffx(k1x:kmx) = twodt*(Avxdudz(k1x:kmx) - Avxdudz(k1x+1:kmx+1)) &
                            *2.*(1.-theta) ! This factor accounts for semi-implicitness
                END SELECT

                !.....Form  [hh]  matrix...........................................
                hh   (k1x:kmx) = hup(k1x:kmx,l)/rhopx(k1x:kmx)
 
                !.....Form [gg]  matrix............................................
                gg(k1x:kmx) = ex(k1x:kmx,l)                                       &
                    -hh(k1x:kmx  ) * (bclncx(k1x:kmx)+wsx0) * tz         &
                    +vdiffx(k1x:kmx)                        * tz
 
                !.....Form [aa]  matrix............................................
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        aa(2,k1x) = 1.0
                    CASE (2:)   ! Two or more wet layers (hu->hup)
                        ! Define upper diagonal terms
                        aa(3,k1x:kmx-1)= -twodt1*Avx(k1x+1:kmx)/(hup(k1x+1:kmx,l)*      &
                            (hup(k1x:kmx-1,l)+hup(k1x+1:kmx,l)))*2.*theta
                        aa(3,kmx)      =  0.0
                        ! Define lower diagonal terms
                        aa(1,k1x+1:kmx)= -twodt1*Avx(k1x+1:kmx)/(hup(k1x:kmx-1,l)*      &
                            (hup(k1x:kmx-1,l)+hup(k1x+1:kmx,l)))*2.*theta
                        aa(1,k1x)      =  0.0
                        ! Define center diagonal terms
                        aa(2,k1x:kmx)  =  1.0                                           &
                            -(hup(k1x-1:kmx-1,l)/hup(k1x:kmx,l))*aa(1,k1x:kmx)    &
                            -(hup(k1x+1:kmx+1,l)/hup(k1x:kmx,l))*aa(3,k1x:kmx)
                END SELECT

                ! ... Top boundary conditions......................................
                ! a. Form wind stress term
                usurf = up(k1x,l)
                vsurf =(vp(k1v(i  ,j  ),        l  ) +                            &
                    vp(k1v(i  ,j-1),    lSC(l) ) +                            &
                    vp(k1v(i+1,j  ),    lEC(l) ) +                            &
                    vp(k1v(i+1,j-1),lSC(lEC(l))))/4.
                cwx = cdw(i,j)*rhoair*SQRT((vair(i,j)-vsurf)**2.+ &
                    (uair(i,j)-usurf)**2.)
                ! b. Modify [gg] matrix
                tausx     = cwx/rhopx(k1x)*uair(i,j)
                gg(  k1x) = gg(  k1x) + tausx*twodt1
                ! c. Modify [aa] matrix
                tausx     = cwx/rhopx(k1x)/hup(k1x,l)
                aa(2,k1x) = aa(2,k1x) + tausx*twodt1

                ! ... Bottom boundary conditions...................................
                ! a. Form bottom stress term
                taubx = cd*SQRT((uhpp(kmx,     l)*uhpp(kmx,        l)) +          &
                    &  ((0.25*(vhpp(kmx,lEC(l))+vhpp(kmx,        l)             &
                    &         +vhpp(kmx,lSC(l))+vhpp(kmx,lSC(lEC(l)))))**2))    &
                    &                /(hup(kmx,l)*hup(kmx,l))
                ! b. Modify [aa] matrix
                aa(2,kmx)   = aa(2,kmx) + taubx*twodt1

                ! .... Point sources and sinks ....................................
                IF ( iopss > 0 ) THEN

                    DO inn = 1, iopss
                        IF (i == ipss(inn)     .AND. &
                            j == jpss(inn)     ) THEN
                            DO k = k1x,kmx
                                IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                                ! ... Strength of Source - here it is assumed that
                                !     only half of the flow shows up in the control volume
                                !     used in the momentum equation -> factor 2 below
                                Usource = ABS(Qpss(k,inn))/(dx*dy*hup(k,l))*twodt1/2. ! not used in Jet model Chris same as follows
                                ! ... DEPENDING ON TYPE OF THE SOURCE/SINK
                                SELECT CASE (ptype(iodev(inn)))
                                    CASE (-2)
                                        Usource = 1.E2
                                        ! ... Velocity of the source in E direction (positive
                                        !     towards east if a source; negative or towards west if
                                        !     a sink) - idetr = 1 by default;
                                        Uhvalue = (Qpss(k,inn)*uEpss(iodev(inn))/dy)*idetr(iodev(inn)) !FJRPlume

                                    ! ... Chris jet model ...
                                    CASE ( 0)
                                        ! ... reconsider Usource to account for the total flux
                                        Usource = ABS(Qpss(k,inn)*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn))))/(dx*dy*hup(k,l))*twodt1/2.
                                        IF (uEpss(iodev(inn)) <=0) THEN         ! define artificial ambient entrainment (as a sink for the face)
                                            Uhvalue = (Qpss(k,inn)*uEpss(iodev(inn))/dy)
                                        ELSE
                                            Uhvalue = Qpss(k,inn)*(uEpss(iodev(inn))-sineagl(iodev(inn)))/dy*idetr(iodev(inn))*2/3   ! approximate um into a uniform u
                                        ENDIF
                                        ! ... overwrite when it is the pump
                                        IF (inn == iopss) THEN                                                   ! only account the pump in the direction normal to the jet discharge FCR
                                            IF (uEpss(iodev(inn)) <=0) THEN
                                                Uhvalue = Qpss(k,inn)                                                        &
                                                &*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn)))  &
                                                &/dy/2.                                                          ! pump as a sink here
                                            ELSE
                                                Uhvalue = 0;
                                            ENDIF
                                        ENDIF
                                    CASE DEFAULT
                                        Uhvalue = (Qpss(k,inn)*uEpss(iodev(inn))/dy)*idetr(iodev(inn))
                                END SELECT
                                aa(2,k) = aa(2,k) + Usource
                                gg(  k) = gg(  k) + Usource * Uhvalue
                            ENDDO
                        ENDIF
                        IF (i == ipss(inn)-1    .AND. &
                            j == jpss(inn)    ) THEN
                            DO k = k1x,kmx
                                IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                                ! ... Strength of Source - here it is assumed that
                                !     only half of the flow shows up in the control volume
                                !     used in the momentum equation -> factor 2 below
                                Usource = ABS(Qpss(k,inn))/(dx*dy*hup(k,l))*twodt1/2.
                                ! ... DEPENDING ON TYPE OF THE SOURCE/SINK
                                SELECT CASE (ptype(iodev(inn)))
                                    CASE (-2)
                                        Usource = 1.E2
                                        ! ... Velocity of the source in W direction (positive
                                        !     towards north if a source; negative or towards south if
                                        !     a sink) - idetr = 1 by default;
                                        Uhvalue = -(Qpss(k,inn)*uWpss(iodev(inn))/dy)*idetr(iodev(inn)) !FJRPlume

                                    ! ... Chris jet model ...
                                    CASE ( 0)
                                        ! ... reconsider Usource to account for the total flux
                                        Usource = ABS(Qpss(k,inn)*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn))))/(dx*dy*hup(k,l))*twodt1/2.
                                        IF (uWpss(iodev(inn)) <=0) THEN ! ... define artificial ambient entrainment (as a sink for the face)
                                            Uhvalue = -(Qpss(k,inn)*uWpss(iodev(inn))/dy)
                                        ELSE
                                            Uhvalue = -Qpss(k,inn)*(uWpss(iodev(inn))-sineagl(iodev(inn)))/dy*idetr(iodev(inn))*2/3  ! approximate um into a uniform u
                                        ENDIF
                                        ! ... overwrite when it is the pump
                                        IF (inn == iopss) THEN                                                   ! only account the pump in the direction normal to the jet discharge FCR
                                            IF (uWpss(iodev(inn)) <=0) THEN
                                                Uhvalue = - Qpss(k,inn)                                                        &
                                                &*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn)))    &
                                                &/dy/2.                                                          ! pump as a sink here
                                            ELSE
                                                Uhvalue = 0;
                                            ENDIF
                                        ENDIF
                                    CASE DEFAULT
                                        Uhvalue = -(Qpss(k,inn)*uWpss(iodev(inn))/dy)*idetr(iodev(inn))
                                END SELECT
                                aa(2,k) =   aa(2,k) + Usource
                                gg(  k) =   gg(  k) + Usource * Uhvalue
                            ENDDO
                        ENDIF
                    ENDDO
                ENDIF

                !.....Solve tridiagonal system for  [ag]  and  [ar]  arrays........
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        ag(k1x) = gg(k1x)/aa(2,k1x)
                        ar(k1x) = hh(k1x)/aa(2,k1x)
                    CASE (2:)   ! Two or more wet layers
                        CALL trid ( aa, gg, hh, ag, ar, k1x, kmx, kmx+1, nwlayers )
                END SELECT

                !.....Save  [ag]  and  [ar]  arrays and sum them over
                !     depth for use in solution of continuity equation.............
                agx(k1x:kmx,l) = ag(k1x:kmx)
                arx(k1x:kmx,l) = ar(k1x:kmx)
                eagx(i,j) = SUM(ag(k1x:kmx))
                earx(i,j) = SUM(ar(k1x:kmx))

            !.....End loop over u-pts.....
            END DO

        !                -----Y-momentum equation-----
        CASE (2)

            !.....Loop over interior v-pts
            !eagy = 0.0; eary = 0.0;
            DO l = 1, lm

                ! ... Map 2D-l into 3D-(i,j) indexes
                i = l2i(l); j = l2j(l);

                ! ... Skip if North Column is dry
                IF (.NOT.mask2d(i,j+1)) CYCLE

                ! ... Compute layer numbers of wet v-pt
                kmy = MIN(kmz(i,j+1), kmz(i,j))
                k1y =                 k1v(i,j)
                nwlayers = (kmy-k1y) + 1
                IF(nwlayers < 1) CYCLE

                ! .... Compute eddy viscosity at interfaces between v-pts)
                Avy = 0.0
                DO k = k1y, kmy
                    Avy(k) = 0.5*(Av(k,lNC(l))+Av(k,l))
                ENDDO

                ! .... Define average layer density at v-pts (in kg/m**3) .........
                rhopy(k1y:kmy) = 1000. ! Neglect vertical density variations
 
                !.....Compute explicit part of water surface slope term ...........
                wsy0 = rhopy(k1y) *  gdtdy  *(spp(i,j+1) - spp(i,j))

                !.... Compute baroclinic term .....................................
                SELECT CASE (ibc)
                    CASE (0)    ! No baroclinic term
                        bclncy(k1:kmy) = 0.0
                    CASE (1:)
                        DO k = k1y, kmy
                            hvpdrho(k) = gdtdy*hvp(k,l)*(rhop(k,lNC(l))-rhop(k,l))
                          ! IF(hp(k,l)<=ZERO .OR. hp(k,lNC(l))<=ZERO) hvpdrho(k) = 0.0
                        ENDDO
                        bclncy(k1y) = hvpdrho(k1y)
                        IF (kmy > k1y) THEN   ! Two or more wet layers
                            ! Recompute bottom layer baroclinic term along horizontal plane
                            CALL bclnc_km (l, kmy, 2, hvpdrho(kmy) )
                            DO k = k1y+1, kmy
                                bclncy(k) = bclncy(k-1) + hvpdrho(k-1) + hvpdrho(k)
                            END DO
                        END IF
                END SELECT

                ! ... Compute explicit portion of vertical diffusion term ........
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        vdiffy(k1y) = 0.0
                    CASE (2:)   ! Two or more wet layers (hvpp->hvp)
                        DO k = k1y+1 , kmy
                            deltaz(k) = hvp(k-1,l) + hvp(k,l)
                            Avydvdz(k)= Avy(k)*(vpp(k-1,l)-vpp(k,l))/deltaz(k)
                        ENDDO
                        Avydvdz(k1y)      = 0.0  ! Set value at free surface to zero
                        Avydvdz(kmy+1)    = 0.0  ! Set value at bottom boundary to zero
                        vdiffy(k1y:kmy) = twodt*(Avydvdz(k1y:kmy) - Avydvdz(k1y+1:kmy+1)) &
                            *2.*(1.-theta) ! This factor accounts for semi-implicitness
                END SELECT


                !.....Form  [hh]  matrix...........................................
                hh   (k1y:kmy) = hvp(k1y:kmy,l)/rhopy(k1y:kmy)

                !.....Compute  [gg]  matrix .......................................
                gg(k1y:kmy) = ex(k1y:kmy,l)                                       &
                    - hh(k1y:kmy  ) * (bclncy(k1y:kmy)+wsy0) * tz         &
                    + vdiffy(k1y:kmy)                        * tz

                !.....Form  [aa]  matrix...........................................
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        aa(2,k1y) = 1.0
                    CASE (2:)   ! Two or more wet layers (hv->hvp)
                        ! Define upper diagonal terms
                        aa(3,k1y:kmy-1)=-twodt1*Avy(k1y+1:kmy)/(hvp(k1y+1:kmy,l)*       &
                            & (hvp(k1y:kmy-1,l)+hvp(k1y+1:kmy,l)))*2.*theta
                        aa(3,kmy)      = 0.0
                        ! Define lower diagonal terms
                        aa(1,k1y+1:kmy)=-twodt1*Avy(k1y+1:kmy)/(hvp(k1y:kmy-1,l)*       &
                            & (hvp(k1y:kmy-1,l)+hvp(k1y+1:kmy,l)))*2.*theta
                        aa(1,k1y)      = 0.0
                        ! Define center diagonal terms
                        aa(2,k1y:kmy)  = 1.0                                            &
                            -(hvp(k1y-1:kmy-1,l)/hvp(k1y:kmy,l))*aa(1,k1y:kmy)  &
                            -(hvp(k1y+1:kmy+1,l)/hvp(k1y:kmy,l))*aa(3,k1y:kmy)
                END SELECT

                ! ... Top boundary conditions .....................................
                ! a. Form wind stress term
                vsurf = vp(k1y,l)
                usurf =(up(k1u(i  ,j  ),        l  ) +                            &
                    up(k1u(i-1,j  ),    lWC(l) ) +                            &
                    up(k1u(i  ,j+1),    lNC(l) ) +                            &
                    up(k1u(i-1,j+1),lWC(lNC(l))))/4. 
                cwy = cdw(i,j)*rhoair*SQRT((uair(i,j)-usurf)**2.+ &
                    (vair(i,j)-vsurf)**2.)
                ! b. Modify [gg] matrix
                tausy = cwy/rhopy(k1y)*vair(i,j)
                gg(k1y)   = gg(  k1y) + tausy*twodt1
                ! c. Modify [aa] matrix
                tausy = cwy/rhopy(k1y)/hvp(k1y,l)
                aa(2,k1y) = aa(2,k1y) + tausy*twodt1
             
                ! ... Bottom boundary conditions ..................................
                ! a. Form bottom stress term
                tauby = cd*SQRT((vhpp(kmy,l)*vhpp(kmy,l)) +                       &
                    &  ((0.25*(uhpp(kmy,lNC(l))+uhpp(kmy,lWC(lNC(l)))           &
                    &         +uhpp(kmy,lWC(l))+uhpp(kmy,l)))**2))              &
                    &         /(hvp(kmy,l)*hvp(kmy,l))
                ! b. Modify [aa] matrix
                aa(2,kmy) = aa(2,kmy) + tauby*twodt1

                ! .... Point sources and sinks ....................................
                IF( iopss > 0 ) THEN

                    DO inn = 1, iopss

                        IF (i == ipss(inn)     .AND. &
                            j == jpss(inn)     ) THEN
                            DO k = k1y,kmy
                                IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                                ! ... Strength of Source - here it is assumed that
                                !     only half of the flow shows up in the control volume
                                !     used in the momentum equation -> factor 2 below
                                Vsource = ABS(Qpss(k,inn))/(dx*dy*hvp(k,l))*twodt1/2.
                                ! ... DEPENDING ON TYPE OF THE SOURCE/SINK
                                SELECT CASE (ptype(iodev(inn)))
                                    CASE (-2)
                                        Vsource = 1.E2
                                        ! ... Velocity of the source in N direction (positive
                                        !     towards north if a source; negative or towards south if
                                        !     a sink) - idetr = 1 by default;
                                        Vhvalue = (Qpss(k,inn)*vNpss(iodev(inn))/dx)*idetr(iodev(inn)) !FJRPlume

                                    ! ... Chris jet model ...
                                    CASE ( 0)
                                        ! ... reconsider Usource to account for the total flux
                                        Vsource = ABS(Qpss(k,inn)*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn))))/(dx*dy*hvp(k,l))*twodt1/2.
                                        IF (vNpss(iodev(inn)) <=0) THEN         ! ... define artificial ambient entrainment (as a sink for the face)
                                            Vhvalue = (Qpss(k,inn)*vNpss(iodev(inn))/dx)
                                        ELSE
                                            Vhvalue = Qpss(k,inn)*(vNpss(iodev(inn))-sineagl(iodev(inn)))/dx*idetr(iodev(inn))*2/3   ! approximate um into a uniform u
                                        ENDIF
                                        ! ... overwrite when it is the pump
                                        IF (inn == iopss) THEN                                                   ! only account the pump in the direction normal to the jet discharge FCR
                                            IF (vNpss(iodev(inn)) <=0) THEN
                                                Vhvalue = Qpss(k,inn)                                                        &
                                                &*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn)))  &
                                                &/dx/2.                                                          ! pump as a sink here
                                            ELSE
                                                Vhvalue = 0;
                                            ENDIF
                                        ENDIF
                                    CASE DEFAULT
                                        Vhvalue = (Qpss(k,inn)*vNpss(iodev(inn))/dx)*idetr(iodev(inn))
                                END SELECT
                                aa(2,k) = aa(2,k) + Vsource
                                gg(  k) = gg(  k) + Vsource * Vhvalue
                            ENDDO
                        ENDIF
                        IF (i == ipss(inn)      .AND. &
                            j == jpss(inn)-1    ) THEN
                            DO k = k1y,kmy
                                IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                                ! ... Strength of Source - here it is assumed that
                                !     only half of the flow shows up in the control volume
                                !     used in the momentum equation -> factor 2 below
                                Vsource = ABS(Qpss(k,inn))/(dx*dy*hvp(k,l))*twodt1/2.
                                ! ... DEPENDING ON TYPE OF THE SOURCE/SINK
                                SELECT CASE (ptype(iodev(inn)))
                                    CASE (-2)
                                        Vsource = 1.E2
                                        ! ... Velocity of the source in S direction (negative
                                        !     towards south if a source; positive or towards north if
                                        !     a sink) - idetr = 1 by default;
                                        Vhvalue = -(Qpss(k,inn)*vSpss(iodev(inn))/dx)*idetr(iodev(inn)) !FJRPlume

                                    ! ... Chris jet model ...
                                    CASE ( 0)
                                        ! ... reconsider Usource to account for the total flux
                                        Vsource = ABS(Qpss(k,inn)*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn))))/(dx*dy*hvp(k,l))*twodt1/2.
                                        IF (vSpss(iodev(inn)) <=0) THEN ! ... define artificial ambient entrainment (as a sink for the face)
                                            Vhvalue = -(Qpss(k,inn)*vSpss(iodev(inn))/dx)
                                        ELSE
                                            Vhvalue = -Qpss(k,inn)*(vSpss(iodev(inn))-sineagl(iodev(inn)))/dx*idetr(iodev(inn))*2/3  ! approximate um into a uniform u
                                        ENDIF
                                        ! ... overwrite when it is the pump
                                        IF (inn == iopss) THEN                                                   ! only account the pump in the direction normal to the jet discharge FCR
                                            IF (vSpss(iodev(inn)) <=0) THEN
                                                Vhvalue = - Qpss(k,inn)                                                      &
                                                &*(uWpss(iodev(inn))+uEpss(iodev(inn))+vNpss(iodev(inn))+vSpss(iodev(inn)))  &
                                                &/dx/2.                                                          ! pump as a sink here
                                            ELSE
                                                Vhvalue = 0;
                                            ENDIF
                                        ENDIF
                                    CASE DEFAULT
                                        Vhvalue = -(Qpss(k,inn)*vSpss(iodev(inn))/dx)*idetr(iodev(inn))
                                END SELECT
                                aa(2,k) = aa(2,k) + Vsource
                                gg(  k) = gg(  k) + Vsource * Vhvalue
                            ENDDO
                        ENDIF
                    ENDDO

                ENDIF

                !.....Solve tridiagonal system for  [ag]  and  [ar]  arrays........
                SELECT CASE (nwlayers)
                    CASE (1)    ! Single wet layer
                        ag(k1y) = gg(k1y)/aa(2,k1y)
                        ar(k1y) = hh(k1y)/aa(2,k1y)
                    CASE (2:)   ! Two or more wet layers
                        CALL trid ( aa, gg, hh, ag, ar, k1y, kmy, km1, nwlayers )
                END SELECT

                !.....Save  [ag]  and  [ar]  arrays and sum them over
                !     depth for use in solution of continuity equation..............
                agy(k1y:kmy,l) = ag(k1y:kmy)
                ary(k1y:kmy,l) = ar(k1y:kmy)
                eagy(i,j) = SUM(ag(k1y:kmy))
                eary(i,j) = SUM(ar(k1y:kmy))
 
            !.....End loop over v-pts.....
            END DO
    END SELECT

    !.....Compute CPU time spent in subroutine.....
    etime = TIMER(0.0)
    t_matmom = t_matmom + (etime - btime)

END SUBROUTINE matmom

!***********************************************************************
SUBROUTINE bclnc_km (l, kb, ieq, huvpdrho )
!***********************************************************************
!
!  Purpose: To recompute the baroclinic term for a bottom layer of
!           variable depth using densities interpolated (or
!           extrapolated) to locations on a horizontal surface.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  18/09/00          P.E. Smith        Original f90 code
!
!-----------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN   ) :: l, kb, ieq
   REAL   , INTENT(INOUT) :: huvpdrho

   !.....Local variables.....
   REAL, PARAMETER :: eps = EPSILON(0.0)
   REAL :: dz1, dz2, rhop_e, rhop_w, rhop_n, rhop_s
   REAL :: hup_plus, hup_minus, hvp_plus, hvp_minus
   LOGICAL :: condition_e, condition_w, condition_n, condition_s

   !.....Choose x- or y-momentum equation...............................
   SELECT CASE ( ieq )


   !            -----Recompute the x-momentum term-----


   CASE (1)
      
      !.....Check if the depth at the press-pt on the east side of the 
      !     control volume is not equal to the depth at the u-pt.....
      hup_plus = hup(kb,l)+eps;  hup_minus = hup(kb,l)-eps
      IF((hp(kb,lEC(l)) > hup_plus) .OR. (hp(kb,lEC(l)) < hup_minus)) THEN
         condition_e = .TRUE.
         dz1 = 0.5*(hp(kb-1,lEC(l)) + hp (kb,lEC(l)))
         dz2 = 0.5*(hp(kb-1,lEC(l)) + hup(kb,    l ))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the u-pt
         rhop_e=rhop(kb-1,lEC(l))+(dz2/dz1)*(rhop(kb,lEC(l))-rhop(kb-1,lEC(l)))
      ELSE
         condition_e = .FALSE.
         rhop_e = rhop(kb,lEC(l))
      END IF

      !.....Check if the depth at the press-pt on the west side of the 
      !     control volume is not equal to the depth at the u-pt.....
      IF((hp(kb,l) > hup_plus) .OR. (hp(kb,l) < hup_minus)) THEN
         condition_w = .TRUE.
         dz1 = 0.5*(hp(kb-1,l) + hp (kb,l))
         dz2 = 0.5*(hp(kb-1,l) + hup(kb,l))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the u-pt
         rhop_w=rhop(kb-1,l)+(dz2/dz1)*(rhop(kb,l)-rhop(kb-1,l))
      ELSE
         condition_w = .FALSE.
         rhop_w = rhop(kb,l)
      END IF

      !.....If necessary, recompute the x-direction
      !     baroclinic term on a horizontal plane.....
      IF ( condition_e .OR. condition_w )  THEN
         huvpdrho = gdtdx*hup(kb,l)*(rhop_e-rhop_w)
      END IF


   !            -----Recompute the y-momentum term-----


   CASE (2)
      
      !.....Check if the depth at the press-pt on the north side of 
      !     the control volume is not equal to the depth at the v-pt.....
      hvp_plus = hvp(kb,l)+eps;  hvp_minus = hvp(kb,l)-eps
      IF((hp(kb,lNC(l)) > hvp_plus) .OR. (hp(kb,lNC(l)) < hvp_minus)) THEN
         condition_n = .TRUE.
         dz1 = 0.5*(hp(kb-1,lNC(l)) + hp (kb,lNC(l)))
         dz2 = 0.5*(hp(kb-1,lNC(l)) + hvp(kb,    l ))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the v-pt
         rhop_n=rhop(kb-1,lNC(l))+(dz2/dz1)*(rhop(kb,lNC(l))-rhop(kb-1,lNC(l)))
      ELSE
         condition_n = .FALSE.
         rhop_n = rhop(kb,lNC(l))
      END IF

      !.....Check if the depth at the press-pt on the south side of 
      !     the control volume is not equal to the depth at the v-pt.....
      IF((hp(kb,l) > hvp_plus) .OR. (hp(kb,l) < hvp_minus)) THEN
         condition_s = .TRUE.
         dz1 = 0.5*(hp(kb-1,l) + hp (kb,l))
         dz2 = 0.5*(hp(kb-1,l) + hvp(kb,l))
         ! Interpolate (or extrapolate) for the density at the
         ! horizontal level of the v-pt
         rhop_s=rhop(kb-1,l)+(dz2/dz1)*(rhop(kb,l)-rhop(kb-1,l))
      ELSE
         condition_s = .FALSE.
         rhop_s = rhop(kb,l)
      END IF

      !.....If necessary, recompute the y-direction
      !     baroclinic term on a horizontal plane.....
      IF ( condition_n .OR. condition_s )  THEN
         huvpdrho = gdtdy*hvp(kb,l)*(rhop_n-rhop_s)
      END IF

   END SELECT

END SUBROUTINE bclnc_km

!***********************************************************************
SUBROUTINE matcon
!***********************************************************************
!
!  Purpose: To calculate the matrix coefficients for solving the
!           continuity equation for zeta.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL :: cx, cy, dt1, dtdx1, dtdy1, rho4sx, rho4sy
   INTEGER :: i, j, k, l, k1s, kms, inn

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Constants.....
   cx     = gdt2dx2*tz*tz; 
   cy     = gdt2dy2*tz*tz;
   dtdx1  = dtdx*tz; 
   dtdy1  = dtdy*tz;
   dt1    = dt*tz; 
   rho4sx = 1000. ! Neglect density variations
   rho4sy = 1000. ! Neglect density variations

   !.....Calculate  [sx] matrix at u-pts & [sy] matrix at v-pts ....
   sx = 0.0; sy = 0.0;
   DO l = 1, lm

     ! ... Map 2D-l into 3D-(i,j) indexes 
     i = l2i(l); j = l2j(l);

     ! ... u-pts
     IF (mask2d(i+1,j)) THEN 
        sx(i,j)= cx * rho4sx * earx(i,j)
     ENDIF

     ! ... v-pts
     IF (mask2d(i,j+1)) THEN
        sy(i,j)= cy * rho4sy * eary(i,j)
     ENDIF

   ENDDO
	  
   !.....Calculate  [dd], [qq], and [rr]  matrices at zeta-pts.....
   dd = 0.0; qq = 0.0; rr = 1.0
   DO l = 1, lm   

     ! ... Map 2D-l into 3D-(i,j) indexes 
     i = l2i(l); j = l2j(l);

     ! ... Top & bottom cells
     k1s = k1
     kms = kmz(i,j)

     ! ... Form matrices
     dd(i,j) = dt1*SUM((uhpp(k1s:kms,l) -         &
             &          uhpp(k1s:kms,lWC(l)))/dx  &
             &        +(vhpp(k1s:kms,l) -         &
             &          vhpp(k1s:kms,lSC(l)))/dy)
     qq(i,j) = spp(i,j) - (dtdx1)*(eagx(i,j)-eagx(i-1,j))   &
             &          - (dtdy1)*(eagy(i,j)-eagy(i,j-1))   &
             &          - dd(i,j)
     rr(i,j) = 1 + sx(i,j) + sx(i-1,j) + sy(i,j) + sy(i,j-1)

   END DO;

   ! .... Modify qq matrices to incorporate sources/sinks
   IF ( iopss > 0 ) THEN
     DO inn = 1, iopss
       IF ( ptype(iodev(inn)) >= 0) CYCLE ! Commented by Chris  for surface elevation only. Bubble plume model will skip the step.
       i = ipss(inn); 
       j = jpss(inn); 
       qq(i,j) = qq(i,j) +  SUM(Qpss(:,inn))/(dx*dy)*twodt*tz
     ENDDO
   ENDIF

   !.....Adjust qq & sx/sy matrices for open boundary conditions.....
   IF (nopen > 0) CALL MODqqddrr4openBC

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_matcon = t_matcon + (etime - btime)

END SUBROUTINE matcon

!***********************************************************************
SUBROUTINE solverBLOCK
!***********************************************************************
!
!  Purpose: To solve the system matrix for zeta using the
!           preconditioned conjugate gradient method.
!
!  More details:
!           The unknowns are ordered from bottom to top of the
!           mesh moving from left to right. The system coefficient
!           matrix is stored in symmetric diagonal format in 'coef'.
!           Each column of 'coef' contains a diagonal of the system
!           matrix and 'jcoef' contains a non-negative integer for
!           each diagonal indicating its distance from the main
!           diagonal. Only the main diagonal, and non-zero diagonals
!           appearing in the upper triangular part of the system
!           matrix, are stored.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   EXTERNAL mic2, cg, si
   INTEGER :: nw1, inw1, maxnz1, i, j, m, ier, nr, nc, istat, nfirst=0

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Set matrix solution variables and workspace
   !     parameters on first entry into subroutine.....
   IF ( nfirst == 0 ) THEN
      nfirst = 1
      ! Compute the number of rows and columns of active cells in the grid
      nr = jlast-jfirst + 1;  nc = ilast-ifirst + 1
      ! Compute number of equations to be solved
      ndim = nr * nc
      ! Define the matrix bandwidth
      ibdwd = nr
      ! Define column width of 'coef' matrix
      maxnz = mdim
      ! Identify which diagonals are stored in the 'coef' matrix
      jcoef(1) = 0;  jcoef(2) = 1;  jcoef(3) = ibdwd
      ! Liberally estimate workspace requirements
      nw = 20*ndim;  inw = 5*ndim
   END IF

   !.....Reset workspace parameters.....
   nw1 = nw+(2*ndim); inw1 = inw; maxnz1 = maxnz

   !.....Allocate arrays.....
   ALLOCATE (coef(ndim,mdim), rhs(ndim), zeta(ndim), wksp(nw1),    &
          &  iwksp(inw1), STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 23 )

   !.....Define  coef  and  rhs  arrays in
   !     preparation for calling .....
   DO i = ifirst, ilast
      DO j = jfirst, jlast
         m = (i-ifirst)*ibdwd + (j-jfirst)+1
         coef(m,1) = rr(i,j)
         coef(m,2) = -sy(i,j)
         coef(m,3) = -sx(i,j)
         rhs(m)    = qq(i,j)
         ! initial guess at zeta
         zeta(m)   = sp(i,j)
      END DO
   END DO

   !.....Adjust coef & rhs for open boundary conditions.....
   IF (nopen > 0) CALL MODcoef4openBC

   !.....Set parameter defaults.....
   CALL dfault ( iparm, rparm )

   !.....Reset some default parameter values.....
   iparm(1) = 2       ! Use the default for DCC runs
   iparm(2) = 500     ! Limit maximum number of iterations to 500 ed. by Chris (200)
   iparm(3) = 1       ! Warning messages and minimum output
   iparm(4) = i6      ! Define fortran unit number for output
   iparm(21)= 0       ! Use scalar algorithm for matrix factorization
   rparm(1) = 1.E-6   ! Try default stopping test value Chris convergence value?? (1.E-6)

   !.....Solve for zeta.....
   WRITE (UNIT=i6,FMT='("**Enter nspcg   n = ", I6)') n
   CALL nspcg (mic2,cg,ndim,mdim,ndim,maxnz1,coef,jcoef,jp,ip,zeta,   &
             & ubar,rhs,wksp,iwksp,nw1,inw1,iparm,rparm,ier)
   WRITE (UNIT=i6,FMT='("**Exit nspcg    n = ", I6)') n
   IF(ier /= 0) WRITE(UNIT=i6,FMT='("**ERROR", I5, " from nspcg")') ier

   !.....STOP program execution if a fatal error is encountered in nspcg.....
   IF (ier < 0 ) THEN
      PRINT *, " "
      PRINT '(" Fatal error in matrix solution on time step = ", I7)', n 
      PRINT '(" **ERROR", I5, " from nspcg")', ier
      PRINT '(" Time = ", F10.4, " hours")', thrs
      PRINT *, " "
      PRINT *, " "
      PRINT *, " ****STOPPING si3d due to fatal error in matrix solution"
      WRITE (UNIT=i6,FMT='(" ****STOPPING si3d due to fatal matrix error")' )
      WRITE (UNIT=i6,FMT='(" Time = ", F10.4, " hours")') thrs 
      STOP
   END IF 

   !.....Load matrix solution into  zeta  array.....
   DO m = 1, ndim
      i = (m-1)/ibdwd + ifirst
      j = m + (jfirst-1) - (i-ifirst)*ibdwd
      s(i,j) = zeta(m)
   END DO

   ! ... Calculate surface layer thickness at time n+1
   CALL slayer_h

   !.....Save the workspace parameters, slightly over-dimensioning it
   !     Otherwise one gets some errors. 
   nw = FLOOR(nw1*1.5); inw = FLOOR(inw1*1.5); maxnz = maxnz1

   !.....Deallocate arrays.....
   DEALLOCATE ( coef, rhs, zeta, wksp, iwksp )

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_solver = t_solver + (etime - btime)

END SUBROUTINE solverBLOCK

!***********************************************************************
SUBROUTINE solverSPARSE
!***********************************************************************
!
!  Purpose: To solve the system matrix for zeta using the
!           preconditioned conjugate gradient method. It uses 
!           Storage format 1 (i.e. ELLPACK or iparm(12) = 1; 
!           in this manner a considerable amount of time is saved
!           as we do not have to store an imxjm by imxjm matrix, 
!           pentadiagonal but very large; instead we only store 
!           lmxlm matrix; this is extremely useful in sparse 
!           bathymetries such as in rivers. Each row in the matrix 
!           has at most 5 non-zero elements, which are stored in coeffA; 
!           the location of the coefficients in the matrix is stored 
!           in jcoefA - see instructions.  
!           
!-----------------------------------------------------------------------

   !.....Local variables.....
   EXTERNAL mic1, jac1, cg, si
   INTEGER :: nw1a, inw1a, maxnz1a, i, j, m, & 
              ier, nrA, ncA, istat, nfirstA=0, ios
   INTEGER, SAVE:: i895 = 895
   CHARACTER(LEN=25):: solverfile="SolverAMODE.txt"

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Set matrix solution variables and workspace
   !     parameters on first entry into subroutine.....
   IF ( nfirstA == 0 ) THEN
      nfirstA = 1
      ! Compute number of equations to be solved
      ndimA = lm
      ! Define column width of 'coef' matrix
      maxnzA = mdimA
      ! Liberally estimate workspace requirements
      nwA = 20*ndimA;  inwA = 5*ndimA
      ! Open file 
      !OPEN (UNIT=i895, FILE=solverfile, IOSTAT=ios)
   END IF

   !.....Reset workspace parameters.....
   nw1a = nwA+(2*ndimA); inw1a = inwA; maxnz1a = maxnzA

   !.....Allocate arrays.....
   ALLOCATE (coeffA(ndimA,mdimA), &
          &  jcoefA(ndimA,mdimA), & 
          &  rhs   (ndimA)      , &
          &  zeta  (ndimA)      , &
          &  wksp  (nw1A )      , &
          &  iwksp (inw1A)      , STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 23 )

   !.....Define  coef, jcoef1  and  rhs  arrays in
   !     preparation for calling .....
   DO m = 1, lm
     i = l2i(m); 
     j = l2j(m);
     coeffA (m,1) =  rr(i,j)
     coeffA (m,2) = -sx(i-1,j  )
     coeffA (m,3) = -sy(i  ,j-1)
     coeffA (m,4) = -sy(i  ,j  )
     coeffA (m,5) = -sx(i  ,j  )
     jcoefA (m,1) =  m
     jcoefA (m,2) =  lWC(m) 
     jcoefA (m,3) =  lSC(m)   
     jcoefA (m,4) =  lNC(m)  
     jcoefA (m,5) =  lEC(m) 
     IF(jcoefA(m,2)==lm1) jcoefA(m,2)=0;
     IF(jcoefA(m,3)==lm1) jcoefA(m,3)=0;
     IF(jcoefA(m,4)==lm1) jcoefA(m,4)=0;
     IF(jcoefA(m,5)==lm1) jcoefA(m,5)=0;
     rhs(m)       =  qq(i,j)
     ! initial guess at zeta
     zeta(m)      = sp(i,j)
   END DO

   !.....Adjust coef & rhs for open boundary conditions.....
   IF (nopen > 0) CALL MODcoefA4openBC

   !.....Set parameter defaults.....
   CALL dfault ( iparm, rparm )

   !.....Reset some default parameter values.....
   iparm(1) = 2       ! Use the default for DCC runs
   iparm(2) = 300     ! Limit maximum number of iterations to 500
   iparm(3) = 1       ! Warning messages and minimum output
   iparm(4) = i6      ! Define fortran unit number for output
   !iparm(21)= 0       ! Use scalar algorithm for matrix factorization
   rparm(1) = 1.E-9   ! Try default stopping test value
   iparm(12) = 1;     ! Storage mode use (1 = Primary format)

   !.....Solve for zeta.....
   WRITE (UNIT=i6,FMT='("**Enter nspcg   n = ", I6)') n
   CALL nspcg (mic1,cg,ndimA,mdimA,ndimA,maxnz1A,coeffA,jcoefA,jp,ip,zeta,   &
             & ubar,rhs,wksp,iwksp,nw1A,inw1A,iparm,rparm,ier)
   WRITE (UNIT=i6,FMT='("**Exit nspcg    n = ", I6)') n
   IF(ier /= 0) WRITE(UNIT=i6,FMT='("**ERROR", I5, " from nspcg")') ier

   !.....STOP program execution if a fatal error is encountered in nspcg.....
   IF (ier < 0 ) THEN
      PRINT *, " "
      PRINT '(" Fatal error in matrix solution on time step = ", I7)', n 
      PRINT '(" **ERROR", I5, " from nspcg")', ier
      PRINT '(" Time = ", F10.4, " hours")', thrs
      PRINT *, " "
      PRINT *, " "
      PRINT *, " ****STOPPING si3d due to fatal error in matrix solution"
      WRITE (UNIT=i6,FMT='(" ****STOPPING si3d due to fatal matrix error")' )
      WRITE (UNIT=i6,FMT='(" Time = ", F10.4, " hours")') thrs 
      STOP
   END IF 

   !.....Load matrix solution into  zeta  array.....
   DO m = 1, lm
     i = l2i(m); 
     j = l2j(m);
     s(i,j) = zeta(m)
   END DO

   ! ... Calculate surface layer thickness at time n+1
   CALL slayer_h

   !.....Save the workspace parameters, slightly over-dimensioning it
   !     Otherwise one gets some errors. - FJR: I really do not get this
   nwA = FLOOR(nw1A*1.5); inwA = FLOOR(inw1A*1.5); maxnzA = maxnz1A

   !.....Deallocate arrays.....
   DEALLOCATE ( coeffA,jcoefA, rhs, zeta, wksp, iwksp )

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_solver = t_solver + (etime - btime)

END SUBROUTINE solverSPARSE

!***********************************************************************
SUBROUTINE slayer_h
!***********************************************************************
!
!  Purpose: To recompute the new values for the surface layer thicknesses
!           (h, hu, hv) after the zeta array is redefined. Note that
!           cells may become dry at time n+1, but no new cells will appear
!           at this time (wetting). The indexes for the surface layers
!           are not updated at this time, since they are used in subr. vel
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l


   ! ... Initialize layer thickness for time n+1 
   h  = hp ; 
   hu = hup; 
   hv = hvp;
   
   ! ... Redo calculations for surface cells
   !     1.- If drying occurs redo calcs. at cell k1s+1 
   !     2.- If cell k1s becomes thicker than its nominal size 
   !         just ignore. Wetting is not done at this time 
   !     Arrays storing surface cells are not modified at this
   !     time since they are used in subr. vel

   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);
      
      ! ... At zeta-points
      k = k1z(i,j)
      h(k,l)=AMIN1(zlevel(k+1),hhs(i,j)) + s(i,j)        
      IF(h (k,l) <= HMIN) THEN
         h  (k,l) = ZERO; 
         k = k + 1; 
         h (k,l) = AMIN1(zlevel(k+1),hhs(i,j)) + s(i,j)
         IF (h  (k,l) <= HMIN) h  (k,l)  = ZERO
      ENDIF

      ! ... At u-points    
      IF (mask2d(i+1,j)) THEN
        k = k1u(i,j)
        hu(k,l)=AMIN1(zlevel(k+1),hhu(i,j)) +             &
        &        MAX(s(i,j),s(i+1,j))
        IF(hu(k,l) <= HMIN) THEN
           hu (k,l) = ZERO; 
           k = k + 1; k = k1u(i,j)
           hu (k,l)=AMIN1(zlevel(k+1),hhu(i,j)) +         &
           &        MAX(s(i,j),s(i+1,j))
           IF(hu (k,l) <= HMIN) hu (k,l)  = ZERO;
        ENDIF
      ENDIF
      
      ! ... At v-points    
      IF (mask2d(i,j+1)) THEN
        k = k1v(i,j)
        hv(k,l)=AMIN1(zlevel(k+1),hhv(i,j)) +             &
        &        MAX(s(i,j),s(i,j+1))
        IF(hv (k,l) <= HMIN) THEN
           hv (k,l) = ZERO; 
           k = k + 1; 
           hv (k,l)=AMIN1(zlevel(k+1),hhv(i,j)) +         &
           &        MAX(s(i,j),s(i,j+1))
           IF(hv (k,l) <= HMIN) hv (k,l)  = ZERO;
        ENDIF
      ENDIF

   ENDDO  

END SUBROUTINE slayer_h

!***********************************************************************
SUBROUTINE layer_h
!***********************************************************************
!
!  Purpose: To recompute the new values for the surface layer thicknesses
!           (h, hu, hv) after the zeta array is redefined. For a linear
!           problem (ilin=1), the surface layer thicknesses are left as
!           constants. (Note: hu and hv on closed boundaries should be
!           considered undefined.)
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kms, kmx, kmy

   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);
      
      ! ... At zeta-points
      kms = kmz(i,j)
      DO k = k1, kms
        h (k,l)=AMIN1(zlevel(k+1),hhs(i,j)) -            &
        &       AMAX1(zlevel(  k),-s (i,j))
        IF(h (k,l) <= HMIN) h (k,l) = ZERO;
      ENDDO
  
      ! ... At u-points    
      IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(i,j),kmz(i+1,j))
        DO k = k1, kmx
          hu (k,l)=AMIN1(zlevel(k+1),hhu(i,j)) -            &
          &        AMAX1(zlevel(  k),-MAX(s(i,j),s(i+1,j)))
          IF(hu(k,l) <= HMIN) hu(k,l) = ZERO;
        ENDDO
      ENDIF
      
      ! ... At v-points    
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(i,j),kmz(i,j+1))
        DO k = k1, kmy
          hv (k,l)=AMIN1(zlevel(k+1),hhv(i,j)) -            &
          &        AMAX1(zlevel(  k),-MAX(s(i,j),s(i,j+1)))
          IF(hv (k,l) <= HMIN) hv (k,l) = ZERO;
        ENDDO
      ENDIF

   ENDDO  

END SUBROUTINE layer_h

!***********************************************************************
SUBROUTINE layer_hp
!***********************************************************************
!
!  Purpose: To recompute the old values for the surface layer thicknesses
!           (hp, hup, hvp) after the zeta array is smoothed. For a linear
!           problem (ilin=1), the surface layer thicknesses are left as
!           constants. (Note: hup and hvp on closed boundaries should be
!           considered undefined.)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, k1s, kms, nwlup, nwlvp, nwlsp

   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);
      
      ! ... At zeta-points
      kms = kmz(i,j)
      k1s = k1z(i,j)
      DO k = k1, kms
          hp (k,l)=AMIN1(zlevel(k+1),hhs(i,j)) -            &
          &        AMAX1(zlevel(  k),-sp(i,j))
          IF(hp(k,l) <= HMIN) hp(k,l) = ZERO;
      ENDDO

      ! ... Compute the array of vertical distances from the
      !     free surface to the top of each layer - Chris
      zfromt(k1s,l) = 0.0
      DO k = k1s+1, kms+1
          zfromt(k,l)    = zfromt(k-1,l) + hp(k-1,l);
      ENDDO
  
      ! ... At u-points    
      IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(i,j),kmz(i+1,j))
        DO k = k1, kmx
          hup(k,l)=AMIN1(zlevel(k+1),hhu(i,j)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(i,j),sp(i+1,j)))
          IF(hup(k,l) <= HMIN) hup(k,l) = ZERO;
        ENDDO
      ENDIF
      
      ! ... At v-points    
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(i,j),kmz(i,j+1))
        DO k = k1, kmy
          hvp(k,l)=AMIN1(zlevel(k+1),hhv(i,j)) -            &
          &        AMAX1(zlevel(  k),-MAX(sp(i,j),sp(i,j+1)))
          IF(hvp(k,l) <= HMIN) hvp(k,l) = ZERO;
        ENDDO
      ENDIF

   ENDDO  

END SUBROUTINE layer_hp

!***********************************************************************
SUBROUTINE TopLayerIndexp
!***********************************************************************
!
!  Purpose: To determine the top layer index given that values of hp, 
!           hup and hvp have been calculated previously
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, nwlup, nwlvp, nwlsp

   k1z = km1; 
   k1u = km1; 
   k1v = km1;
   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);
      
      ! ... At zeta-points
      kms = kmz(i,j)
      DO k = k1, kms
        IF(hp (k,l) > ZERO) THEN
          k1z(i,j) = k
          EXIT
        ENDIF
      ENDDO
  
      ! ... At u-points    
      IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(i,j),kmz(i+1,j))
        DO k = k1, kmx
          IF(hup(k,l) > ZERO) THEN
            k1u(i,j) = k
            EXIT
          ENDIF
        ENDDO
      ENDIF
      
      ! ... At v-points    
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(i,j),kmz(i,j+1))
        DO k = k1, kmy
          IF(hvp(k,l) > ZERO) THEN
            k1v(i,j) = k
            EXIT
          ENDIF
        ENDDO
      ENDIF

   ENDDO  

END SUBROUTINE TopLayerIndexp

!***********************************************************************
SUBROUTINE vel
!***********************************************************************
!
!  Purpose: To solve the momentum equations explicitly for velocity.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL, DIMENSION(:,:), POINTER :: cxx, cyy
   REAL :: gthx1, gthy1, rho4cxx, rho4cyy, vvtemp, uutemp
   INTEGER :: i, j, k, l, istat, kmx, kmy, kms, k0x, k0y, k0s, k1s, k1ss

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Constants.....
   gthx1 = gdtdx*tz; gthy1 = gdtdy*tz

   !                -----X-momentum equation-----

   !.....Allocate cxx pointer array for temporary storage.....
   ALLOCATE (cxx(im1,jm1), STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 24 )

   ! ... Define constants: Ignore density variations in z 
   rho4cxx = 1000. 
   rho4cyy = 1000. 

   ! ... Loop over cells
   DO l = 1, lm
     
      ! .... Map l into 2D-xy space
      i = l2i(l); j = l2j(l);

      ! ... Skip if East Column is dry
      IF(.NOT.mask2d(i+1,j)) CYCLE
                 
      !.....Top & Bottom wett u-points
      kmx = MIN(kmz(i+1,j),kmz(i,j))
      k0x = k1u(i,j)

      !.....Solve for the water surface slope portion of the
      !     x-mom eq and save the result in the cxx array.....
      cxx(i,j) = gthx1 * rho4cxx * (s(i+1,j)-s(i,j))

      !.....Solve the x-momentum equation for uh.....
      DO k = k0x, kmx
        uh(k,l) =  agx(k,l)    - cxx(i,j)*arx(k,l)
      ENDDO

      ! ... Redo near surface flux calcs. if Drying occurs
      !     hu is calculated after the solution of zeta
      !     in subr. solver. The top most layer during
      !     time n (n+1/2) remains the same through the 
      !     calculations to predict n+1 from n or n+1/2
      k = k0x;
      IF (hu(k  ,l) <= ZERO) THEN
        ! ... Update fluxes
        uh(k+1,l) = uh(k,l)+uh(k+1,l)
        uh(k  ,l) = 0.0
        ! ... Update surface array
        k1u(i,j) = k0x+1
      ENDIF
  	      
   END DO

   !.....Deallocate pointer cxx
   DEALLOCATE ( cxx )

   
   !                -----Y-momentum equation-----
   
   
   !.....Allocate cyy pointer array for temporary storage.....
   ALLOCATE (cyy(im1,jm1), STAT=istat)
   IF (istat /= 0) CALL allocate_error ( istat, 25 )

   ! ... Loop over cells
   DO l = 1, lm
     
      ! .... Map l into 2D-xy space
      i = l2i(l); j = l2j(l);

      ! ... Skip if North Column is dry
      IF(.NOT.mask2d(i,j+1)) CYCLE

      !.....Top & Bottom wett v-points .....
      kmy = MIN(kmz(i,j+1),kmz(i,j))
      k0y = k1v(i,j)
                        
      !.....Solve for the water surface slope portion of the
      !     y-mom eq and save the result in the cyy array.....
      cyy(i,j) = gthy1 * rho4cyy * (s(i,j+1)-s(i,j))
         
      !.....Solve the y-momentum equation for vh.....
      DO k = k0y, kmy
        vh(k,l) =  agy(k,l)    - cyy(i,j)*ary(k,l)
      ENDDO 

      ! ... Redo near surface flux calcs. if Drying occurs
      !     hu is calculated after the solution of zeta
      !     in subr. solver. The top most layer during
      !     time n (n+1/2) remains the same through the 
      !     calculations to predict n+1 from n or n+1/2
      k = k0y;
      IF (hv(k  ,l) <= ZERO) THEN
        ! ... Update fluxes
        vh(k+1,l) = vh(k,l)+vh(k+1,l)
        vh(k  ,l)   = 0.0
        ! ... Update surface array
        k1v(i,j) = k0y+1
      ENDIF
       
   END DO

   !.....Deallocate pointer cyy
   DEALLOCATE ( cyy )

   !.... Calculate vertical velocities from uh & uhpp values
   !     to be used in the solution of scalar transport eq.
   CALL continuity(1)

   ! ... Update surface array at s-points if drying occurs
   DO l = 1, lm
      ! .... Map l into 2D-xy space
      i = l2i(l); j = l2j(l);
      ! ... Top wett s-points
      k0s = k1z(i,j);
      ! ... Recalculate surface cell if water column
      !     is wett. 
      IF (k0s<km1 .AND. h(k0s,l)<=ZERO) k1z(i,j)=k0s+1
   ENDDO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_vel = t_vel + (etime - btime)

END SUBROUTINE vel

!***********************************************************************
SUBROUTINE vel2
!***********************************************************************
!
!  Purpose: To recompute velocities at new time layer if wetting occurs
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL :: vvtemp, uutemp
   INTEGER :: i, j, k, l

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Recalculate layer thicknesses including wetting & drying
   CALL layer_h

   ! ... Loop over cells
   DO l = 1, lm
     
     ! .... Map l into 2D-xy space
     i = l2i(l); j = l2j(l);

     ! ... Skip if East Column is dry
     IF(mask2d(i+1,j)) THEN       
       k = k1u(i,j);
       IF (hu(k-1,l) >  ZERO) THEN
         uutemp    = uh(k,l)/(hu(k,l)+hu(k-1,l))
         uh(k  ,l) = uutemp * hu(k  ,l);
         uh(k-1,l) = uutemp * hu(k-1,l);
       ENDIF
     ENDIF
     
     ! ... Skip if North Column is dry
     IF(mask2d(i,j+1)) THEN
       k = k1v(i,j);
       IF (hv(k-1,l) >  ZERO) THEN
         vvtemp    = vh(k,l)/(hv(k,l)+ hv(k-1,l))
         vh(k  ,l) = vvtemp * hv(k  ,l);
         vh(k-1,l) = vvtemp * hv(k-1,l);
       ENDIF	               
     ENDIF
	       
   END DO

   !.....No need to recalculate vertical velocities since
   !     they are calculated in either save or setrap routines for
   !     next iteration or next time setp

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_vel = t_vel + (etime - btime)

END SUBROUTINE vel2


!***********************************************************************
SUBROUTINE continuity (ist)
!***********************************************************************
!
!  Purpose: To compute wp from horizontal components of the velocity 
!           field
!
!-----------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT(IN):: ist
  
   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, inn, nn

   SELECT CASE (ist)

   CASE (1) ! Compute wp from uh, uhpp, vh and vhpp

     DO l = 1, lm

       ! .... Map l into 2D-xy space
       i = l2i(l); j = l2j(l);

       ! ... Bottom wett s-points
       kms = kmz(i,j);
       k1s = k1z(i,j); 

       ! ... Cycle if water column is dry
       IF( k1s > kms ) CYCLE

       ! .... Loop over cells in water colum to estimate vertical 
       !      velocities which are consistent with the formulation of 
       !      continuity (mass conservation). These velocities are then
       !      used for scalar transport calculations
       wp(:,l) = 0.0
       DO k = kms,k1s,-1
         wp(k,l) = wp  (k+1,l)                         &
             &   -(uh  (k  ,l)-uh  (k  ,lWC(l))+       &
                   uhpp(k  ,l)-uhpp(k  ,lWC(l)))/twodx &
             &   -(vh  (k  ,l)-vh  (k  ,lSC(l))+       &
                   vhpp(k  ,l)-vhpp(k  ,lSC(l)))/twody
       ENDDO

       ! ... Correct wp estimates for surface cell, due to
       !     advective flux from neighbouring cells above it
       DO k = k1,k1s-1
         wp(k1s,l) = wp  (k1s,l)                         &
               &   -(uh  (k  ,l)-uh  (k  ,lWC(l))+       &
                     uhpp(k  ,l)-uhpp(k  ,lWC(l)))/twodx &
               &   -(vh  (k  ,l)-vh  (k  ,lSC(l))+       &
                     vhpp(k  ,l)-vhpp(k  ,lSC(l)))/twody
       ENDDO

     ENDDO

     ! .... Modify w estimates to incorporate sources/sinks (PSS)
     ! .... Modify wp to account for water jet model (Chris)
     IF ( iopss > 0 ) THEN
       DO inn = 1, iopss
         i = ipss(inn) 
         j = jpss(inn)
         l = ij2l(i,j); 
         kms = kmz(i,j); 
         k1s = k1z(i,j);
         nn  = iodev(inn);
         DO k = kms,k1s,-1
             wp(k,l) = wp(k,l) + SUM(Qpss(k:kms,inn))/(dx*dy)         ! Chris
         ENDDO
         IF (ptype(nn) == 0 .AND. inn < iopss ) THEN                  ! determine if it is the jet nozzle source ksrc=kms-jetcell
                                                                      ! exclude the pump location added on 22/01/15 Chris
             ! ... east exit of the jet cell
             IF (uEpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lEC(l)) =                                   &
                         &          wp(k,lEC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... west exit of the jet cell
             IF (uWpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lWC(l)) =                                   &
                         &          wp(k,lWC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... north exit of the jet cell
             IF (vNpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lNC(l)) =                                   &
                         &          wp(k,lNC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... south exit of the jet cell
             IF (vSpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lSC(l)) =                                   &
                         &          wp(k,lSC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... vertical entrainment
             ! Lower boundary
             wp(kms-jetcell-jettopc+1:kms-jetcell,l) =                 &
                 &        wp(kms-jetcell-jettopc+1:kms-jetcell,l)      & ! flow from the cell below (if wp is positive, flow is going upward)
                 &        + wLmtn(nn)                                  & ! coefficient of flow momentum into jet sources (normally negative value)
                 &        *Qpss(kms-jetcell,inn)/(dx*dy)                 ! one source two jets, use 2 here
             ! Upper boundary
             wp(kms-jetcell-jettopc+1:kms-jetcell,l) =                 &
                 &        wp(kms-jetcell-jettopc+1:kms-jetcell,l)      & ! flow from the cell below
                 &        + wUmtn(nn)                                  & ! coefficient of flow into jet sources (one source two jets, use 2 here)
                 &        *Qpss(kms-jetcell,inn)/(dx*dy)


         ENDIF   
       ENDDO
     ENDIF

   ! ... Calculate wp in terms of uhp and vhp values - for computations
   !     of velocities at next time step
   CASE (2)

     DO l = 1, lm 

       ! .... Map l into 2D-xy space
       i = l2i(l); j = l2j(l);

       ! ... Bottom wett s-points
       kms = kmz(i,j);
       k1s = k1z(i,j); 
       
       ! ... Cycle if water column is dry
       IF( k1s > kms ) CYCLE

       ! .... Loop over cells in water colum to estimate vertical 
       !      velocities 
       wp(:,l) = 0.0
       DO k = kms,k1s,-1
         wp(k,l) = wp(k+1,l)-(uhp(k,l)-uhp(k,lWC(l)))/dx    &
                  &         -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

       ! ... Correct wp estimates for surface cell, due to
       !     advective flux from neighbouring cells above it
       DO k = k1,k1s-1
         wp(k1s,l) = wp(k1s,l)-(uhp(k,l)-uhp(k,lWC(l)))/dx    &
                  &           -(vhp(k,l)-vhp(k,lSC(l)))/dy
       END DO

     END DO

     ! .... Modify w estimates to incorporate sources/sinks (PSS)
     ! .... Modify wp to account for water jet model (Chris)
     IF ( iopss > 0 ) THEN
       DO inn = 1, iopss
         i = ipss(inn) 
         j = jpss(inn)
         l = ij2l(i,j); 
         kms = kmz(i,j); 
         k1s = k1z(i,j);
         nn  = iodev(inn);
         DO k = kms,k1s,-1
             wp(k,l) = wp(k,l) + SUM(Qpss(k:kms,inn))/(dx*dy)         ! Chris
         ENDDO
         IF (ptype(nn) == 0 .AND. inn < iopss ) THEN                  ! determine if it is the jet nozzle source ksrc=kms-jetcell
                                                                      ! exclude the pump location added on 22/01/15 Chris
             ! ... east exit of the jet cell
             IF (uEpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lEC(l)) =                                   &
                         &          wp(k,lEC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... west exit of the jet cell
             IF (uWpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lWC(l)) =                                   &
                         &          wp(k,lWC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... north exit of the jet cell
             IF (vNpss(nn) > 0) THEN
                 DO k = kms-jetcell, kms-jetcell-jettopc+1, -1
                     wp(k,lNC(l)) =                                   &
                         &          wp(k,lNC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... south exit of the jet cell
             IF (vSpss(nn) > 0) THEN
                 DO k = kms-jetcell, k1s, -1
                     wp(k,lSC(l)) =                                   &
                         &          wp(k,lSC(l))                      &
                         &        + sineagl(nn)                       & ! flow angle
                         &        *Qpss(kms-jetcell,inn)/(dx*dy)
                 ENDDO
             ENDIF

             ! ... vertical entrainment
             ! Lower boundary
             wp(kms-jetcell-jettopc+1:kms-jetcell,l) = &
                 &        wp(kms-jetcell-jettopc+1:kms-jetcell,l)      & ! flow from the cell below (if wp is positive, flow is going upward)
                 &        + wLmtn(nn)                                  & ! coefficient of flow momentum into jet sources (normally negative value)
                 &        *Qpss(kms-jetcell,inn)/(dx*dy)                 ! one source two jets, use 2 here
             ! Upper boundary
             wp(kms-jetcell-jettopc+1:kms-jetcell,l) =                 &
                 &        wp(kms-jetcell-jettopc+1:kms-jetcell,l)      & ! flow from the cell below
                 &        + wUmtn(nn)                                  & ! coefficient of flow into jet sources (one source two jets, use 2 here)
                 &        *Qpss(kms-jetcell,inn)/(dx*dy)


         ENDIF     
       ENDDO
     ENDIF

   END SELECT

   ! ... Modify velocity estimates near the boundaries to 
   !     account for open boundaries.
   IF (nopen > 0) CALL MODvel4openBC

 END SUBROUTINE continuity

!***********************************************************************
SUBROUTINE exsal
!***********************************************************************
!
!  Purpose: To evaluate the explicit terms (advection and horizontal
!           diffusion) in the salinity equation. The sum of these
!           terms are saved in the 3-d array  ex(k,l)  which is the
!           primary output from this subroutine.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL    :: adv, hd, hti, twodt1, advk1s, scC,   & 
              scE, scW, scN, scS, scU, scD,        &
              uE , uW , vN , vS , wU , wD  
   INTEGER :: i, j, k, l, istat, kms, k1s, nwlayers

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Solve using Flux-Limiters if requested. In that case just
   !     return from this subroutine.
   IF (itrsca == 4) THEN 
     CALL exsalFL; 
     RETURN
   ENDIF

   ! ... Constants used in solution
   twodt1 = twodt*tz

   ! ... Calculate hdxpp & hdypp arrays for diffusion terms & 
   !               th    & th1   arrays for advection terms
   hdxpp = 0.0; hdypp = 0.0; th = 0.0; th1 = 0.0;
   DO l = 1, lm

      ! ... 3D-(i,j) indexes for l
      i = l2i(l); j = l2j(l);

      ! ... Retrieve top & bottom wet sal-pts .................
      kms = kmz(i,j)
      k1s = k1z(i,j)

      ! ... Calculate hdxpp & hdypp array at u-&v- pts ........
      !     Note that interfaces connecting wett & dry cells  
      !     DO NOT have diffussive transport  
      DO k = k1, kms
        hdxpp(k,l) = Ax0*hupp(k,l)
        hdypp(k,l) = Ay0*hvpp(k,l)
      ENDDO
                      
      !.....Calculate weighting arrays for vertical advection 
      DO k = k1s+1, kms
         th (k,l) = hp(k-1,l)/(hp(k-1,l)+hp(k,l))
         th1(k,l) = 1.-th(k,l)
      ENDDO

   END DO

   !.....Calculate the explicit part of the scalar transport equation.....
   ex = 0.0;
   DO l = 1, lm
            
      ! ... Map l- into 3D-(i,j) indexes 
      i = l2i(l); j = l2j(l);
            
      !..... Top & bottom s-pts.....
      kms = kmz(i,j); 
      k1s = k1z(i,j);

      ! ... Cycle DRY cells
      IF ( k1s > kms ) CYCLE
     
      ! ... Loop over cells in water column
      DO k = k1s, kms
  
        ! ... Define velocities (2 x vel) at EWNS faces 
        uE = (uh(k,    l )+uhpp(k,    l ))
        uW = (uh(k,lWC(l))+uhpp(k,lWC(l)))
        vN = (vh(k,    l )+vhpp(k,    l ))
        vS = (vh(k,lSC(l))+vhpp(k,lSC(l)))

        ! ... Define velocities at UD faces
        wU = wp (k  ,l); IF(k == k1s) wU = 0.0;
        wD = wp (k+1,l); IF(k == kms) wD = 0.0;
        
        ! ... Calculate advective term
        SELECT CASE (itrsca)

        CASE (1) ! Centered scheme 

          ! ... Scalars at EWNS faces
          scE = (salp(k,lEC(l)) + salp(k,l)) 
          scW = (salp(k,lWC(l)) + salp(k,l)) 
          scN = (salp(k,lNC(l)) + salp(k,l)) 
          scS = (salp(k,lSC(l)) + salp(k,l)) 

          ! ... Scalars at UD faces
          scU = th (k  ,l)*salp(k  ,l)+th1(k  ,l)*salp(k-1,l) 
          scD = th (k+1,l)*salp(k+1,l)+th1(k+1,l)*salp(k  ,l)  
            
          ! ... Advective term 
          adv = ( uE * scE - uW * scW ) / fourdx +           &
                ( vN * scN - vS * scS ) / fourdy +           &
                ( wU * scU - wD * scD )      

        CASE (2) ! Upwind scheme 

          adv  = ( (uE+ABS(uE)) * salpp(k,    l ) +          &
                   (uE-ABS(uE)) * salpp(k,lEC(l)) -          &
                   (uW+ABS(uW)) * salpp(k,lWC(l)) -          &
                   (uW-ABS(uW)) * salpp(k,    l ) ) / fourdx &
                +( (vN+ABS(vN)) * salpp(k,    l ) +          &
                   (vN-ABS(vN)) * salpp(k,lNC(l)) -          &
                   (vS+ABS(vS)) * salpp(k,lSC(l)) -          &
                   (vS-ABS(vS)) * salpp(k,    l ) ) / fourdy & 
                 -((wD+ABS(wD)) * salpp(k+1,l) +             &
                   (wD-ABS(wD)) * salpp(k  ,l)) / 2.         &
                 +((wU+ABS(wU)) * salpp(k  ,l) +             &
                   (wU-ABS(wU)) * salpp(k-1,l)) / 2.

        CASE (3) ! Centered in all faces but in vertical flux
                 ! faces located next below the free surface 

          ! ... Scalars at EWNS faces
          scE = (salp(k,lEC(l)) + salp(k,l)) 
          scW = (salp(k,lWC(l)) + salp(k,l)) 
          scN = (salp(k,lNC(l)) + salp(k,l)) 
          scS = (salp(k,lSC(l)) + salp(k,l)) 
  
          ! ... Scalars at UD faces
          scU = th(k  ,l)*salp(k  ,l)+th1(k  ,l)*salp(k-1,l) 
          scD = th(k+1,l)*salp(k+1,l)+th1(k+1,l)*salp(k  ,l)  
              
          ! ... Horizontal advection 
          adv= ( uE * scE - uW * scW ) / fourdx +            &
               ( vN * scN - vS * scS ) / fourdy            
  
          ! ... Vertical advection - Note that 1.E-1 is an
          !     an 'arbitrary' small value. One can choose 
          !     a priori any other small positive value
          IF      (k == k1s   .AND. h(k1s,l) < SMALL) THEN
            adv=adv-((wD+ABS(wD)) * salpp(k+1,l) +           &
                     (wD-ABS(wD)) * salpp(k  ,l)) / 2.
          ELSE IF (k == k1s+1 .AND. h(k1s,l) < 1.E-1) THEN 
            adv=adv+((wU+ABS(wU)) * salpp(k  ,l) +           &
                     (wU-ABS(wU)) * salpp(k-1,l)) / 2.       &
                   -  wD*scD
          ELSE
            adv=adv+wU*scU-wD*scD 
          ENDIF
                                   
        END SELECT

        !.....Horizontal diffusion.....
        hd= (hdxpp(k,    l )*(salpp(k,lEC(l)) - salpp(k,    l ))       &
          & -hdxpp(k,lWC(l))*(salpp(k,    l ) - salpp(k,lWC(l))))/dxdx &
          &+(hdypp(k,    l )*(salpp(k,lNC(l)) - salpp(k,    l ))       &
          & -hdypp(k,lSC(l))*(salpp(k,    l ) - salpp(k,lSC(l))))/dydy

        !.....Sum all terms (Limiting the value of advection)
        ex(k,l) = hpp(k,l)*salpp(k,l)/twodt1-adv+hd*ihd

      END DO

      ! .... Modify ex(k1s,l) term for horizontal fluxes from cells above k1s
      advk1s = 0.0E0
      DO k = k1, k1s-1

        !.....Horizontal diffusion.....
        hd= (hdxpp(k,    l )*(salpp(k,lEC(l)) - salpp(k,    l ))       &
          & -hdxpp(k,lWC(l))*(salpp(k,    l ) - salpp(k,lWC(l))))/dxdx &
          &+(hdypp(k,    l )*(salpp(k,lNC(l)) - salpp(k,    l ))       &
          & -hdypp(k,lSC(l))*(salpp(k,    l ) - salpp(k,lSC(l))))/dydy

        ! ... Velocity at faces
        uE = uh(k,    l )+uhpp(k,    l )
        uW = uh(k,lWC(l))+uhpp(k,lWC(l))
        vN = vh(k,    l )+vhpp(k,    l )
        vS = vh(k,lSC(l))+vhpp(k,lSC(l))
        
        ! ... Calculate advective term
        SELECT CASE (itrsca)

        CASE (1) ! Centered scheme 

          ! ... Scalars at EWNS faces
          scE = (salp(k,lEC(l)) + salp(k,l)) 
          scW = (salp(k,lWC(l)) + salp(k,l)) 
          scN = (salp(k,lNC(l)) + salp(k,l)) 
          scS = (salp(k,lSC(l)) + salp(k,l)) 
          
          ! ... Advective term 
          advk1s = advk1s +                                  &
                   ( uE * scE - uW * scW ) / fourdx +        &
                   ( vN * scN - vS * scS ) / fourdy - hd *ihd 

        CASE (2) ! Upwind scheme 

          advk1s = advk1s +                                  &
                 ( (uE+ABS(uE)) * salpp(k,    l ) +          &
                   (uE-ABS(uE)) * salpp(k,lEC(l)) -          &
                   (uW+ABS(uW)) * salpp(k,lWC(l)) -          &
                   (uW-ABS(uW)) * salpp(k,    l ) ) / fourdx &
                +( (vN+ABS(vN)) * salpp(k,    l ) +          &
                   (vN-ABS(vN)) * salpp(k,lNC(l)) -          &
                   (vS+ABS(vS)) * salpp(k,lSC(l)) -          &
                   (vS-ABS(vS)) * salpp(k,    l ) ) / fourdy - hd *ihd

        CASE (3) ! Centered in all faces but in vertical flux
                 ! faces located next below the free surface 

          ! ... Scalars at EWNS faces
          scE = (salp(k,lEC(l)) + salp(k,l)) 
          scW = (salp(k,lWC(l)) + salp(k,l)) 
          scN = (salp(k,lNC(l)) + salp(k,l)) 
          scS = (salp(k,lSC(l)) + salp(k,l)) 
               
          ! ... Horizontal advection 
          advk1s = advk1s +                                   &
               ( uE * scE - uW * scW ) / fourdx +             &
               ( vN * scN - vS * scS ) / fourdy - hd *ihd           
                                  
        END SELECT

      ENDDO

      ! ... Modify ex- values for surface cells
      ex(k1s,l) = ex(k1s,l) - advk1s

   ENDDO

   ! ... Modify explicit term to account for flow boundary conditions
   CALL MODexsal4openbc

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exsal = t_exsal + (etime - btime)

END SUBROUTINE exsal

!***********************************************************************
SUBROUTINE exsalFL
!***********************************************************************
!
!  Purpose: To evaluate the explicit terms (advection) in the salinity 
!           equation using flux limiter methods. The sum of these
!           terms are saved in the array  ex(k,l)  which is the
!           primary output from this subroutine.
!
!-----------------------------------------------------------------------

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, gamma1, istat
   REAL    :: vel, ratio, C_f, delz, twodt1, hd
   REAL, DIMENSION (4          ) :: ss

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !!.....Allocate hdxpp, hdypp, th, and th1 arrays.....
   !ALLOCATE (hdxpp(km1,lm1), hdypp(km1,lm1), STAT=istat)
   !IF ( istat /= 0) CALL allocate_error ( istat, 26 )

   ! ... Constants used in solution
   twodt1 = twodt*tz

   ! ... Calculate hdxpp & hdypp arrays for diffusion terms & 
   hdxpp = 0.0; hdypp = 0.0;
   DO l = 1, lm

      ! ... 3D-(i,j) indexes for l
      i = l2i(l); j = l2j(l);

      ! ... Retrieve top & bottom wet sal-pts .................
      kms = kmz(i,j)
      k1s = k1z(i,j)

      ! ... Calculate hdxpp & hdypp array at u-&v- pts ........
      !     Interfaces connecting wett & dry cells will not 
      !     have diffussive transport in present formulation 
      DO k = k1, kms
        hdxpp(k,l) = Ax0*hupp(k,l)
        hdypp(k,l) = Ay0*hvpp(k,l)
      ENDDO
                        
   END DO

   ! ... Initialize ex & flux arrays to zeros
   ex = 0.0; fluxX= 0.0; fluxY = 0.0; fluxZ = 0.0;

   DO l = 1, lm; 

     ! ... Map l- into (i,j)-indexes .........................
     i = l2i(l); j = l2j(l);
	 
     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(i,j)
     k1s = k1z(i,j)
 
     DO k = k1s, kms;
	    
       ! ... EW fluxes .......................................
       IF (hup(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2  
         vel  = (uhpp(k,l) + uh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = salpp(k,    l ); 
         ss(3)  = salpp(k,lEC(l)); 
         IF (hpp(k,    lWC(l) )<=ZERO) THEN; ss(1) = ss(2); 
          ELSE; ss(1)=salpp(k,    lWC(l) ); ENDIF;
         IF (hpp(k,lEC(lEC(l)))<=ZERO) THEN; ss(4) = ss(3); 
          ELSE; ss(4)=salpp(k,lEC(lEC(l))); ENDIF;

         ! ... Calculate Cf for flux computation
         C_f    = 0.0; 
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF
 
         ! ... Calculate fluxes
         fluxX(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.

       ENDIF

       ! ... NS fluxes .......................................
       IF (hvp(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2  
         vel  = (vhpp(k,l) + vh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = salpp(k,        l  ); 
         ss(3)  = salpp(k,    lNC(l) ); 
         IF (hpp(k,    lSC(l) )<= ZERO) THEN; ss(1) = ss(2);  
           ELSE; ss(1)  = salpp(k,    lSC(l) ); ENDIF;
         IF (hpp(k,lNC(lNC(l)))<= ZERO) THEN; ss(4) = ss(3);  
           ELSE; ss(4)  = salpp(k,lNC(lNC(l))); ENDIF;

         ! ... Calculate Cf for flux computation
         C_f    = 0.0
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF
 
         ! ... Calculate fluxes
         fluxY(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.

       ENDIF

       ! ... UD fluxes .......................................
       IF (hp(k-1,l) > ZERO) THEN

         ! ... Velocities at time n + 1
         vel  = wp(k,l); IF (k == k1s) vel = 0.0;

         ! ... Define stencil for scalar transport
         ss(2) = salpp(k  ,l); 
         ss(3) = salpp(k-1,l); 
         IF (hpp(k-2,l)<=ZERO) THEN ; ss(4)=ss(3);
            ELSE; ss(4)=salpp(k-2,l); ENDIF;
         IF (hpp(k+1,l)<=ZERO) THEN ; ss(1)=ss(2);
            ELSE; ss(1)=salpp(k+1,l); ENDIF;  

         ! ... Define C_f for flux computations
         C_f   = 1.0 ! Default method is Lax-Wendroff
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))  
           ! MC flux limiter (VanLeer, 1977)   
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes
         delz = (hp(k,l) + hp(k-1,l))/2.
         fluxZ(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/delz*C_f)*(ss(3)-ss(2))/2.

       ENDIF

     ENDDO; 
   ENDDO; 

   ! ... Update ex array with x-flux divergence
   DO l = 1, lm; 

     ! ... Map l- into (i,j)-indexes .........................
     i = l2i(l); j = l2j(l);
	 
     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(i,j)
     k1s = k1z(i,j)
 
     DO k = k1s, kms;	   

        !.....Horizontal diffusion.....
        hd= (hdxpp(k,    l )*(salpp(k,lEC(l)) - salpp(k,    l ))       &
          & -hdxpp(k,lWC(l))*(salpp(k,    l ) - salpp(k,lWC(l))))/dxdx &
          &+(hdypp(k,    l )*(salpp(k,lNC(l)) - salpp(k,    l ))       &
          & -hdypp(k,lSC(l))*(salpp(k,    l ) - salpp(k,lSC(l))))/dydy

        !.....Sum all terms 
        ex(k,l) =   hpp(k,l)*salpp(k,l)/twodt1        &  
                - (fluxX(k,l) - fluxX(k,lWC(l))) / dx &
                - (fluxY(k,l) - fluxY(k,lSC(l))) / dy &
                - (fluxZ(k,l) - fluxZ(k+1,l   )) + hd * ihd

     ENDDO; 

   ENDDO

   ! ... Modify explicit term to account for flow boundary conditions
   CALL MODexsal4openbc

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exsal = t_exsal + (etime - btime)

END SUBROUTINE exsalFL

!***********************************************************************
SUBROUTINE imsal
!***********************************************************************
!
!  Purpose: To solve for active scalar concentration.
!
!  29-may-2009	(F.J.Rueda)	 Include tpload (temp. load associated to rwps)
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   REAL :: twodt1, Tsource, Qsource
   INTEGER :: i, j, k, l, k1s, kms, kt, nwlayers, inn, nn, kk, noc
   REAL, DIMENSION (1:km1) :: hn
   REAL, DIMENSION (3,1:km1) :: aa

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over interior sal-pts to solve for
   !     matrix from the active scalar equation.....
   DO l = 1, lm
            
      ! ... 3D-(i,j) indexes for l
      i = l2i(l); j = l2j(l);
            
      !.....Compute top & bottom layer numbers & No. of layers ....
      kms = kmz(i,j);
      k1s = k1z(i,j);
      nwlayers = kms-k1s+1

      ! ... Define layer thikness at time n - The corrections for
      !     surface and recently submerged cells are needed to 
      !     keep mass conservation -  The test used to check
      !     mass conservation is that of a surface seiche with 
      !     an equilibrium water surface level at the level of 
      !     where grids cells change from level k to k+1 
      hn(k1s+1:kms) = h(k1s+1:kms,l)
      hn(k1s      ) = twodt1*wp(k1s,l)+hpp(k1s,l)
      IF (hpp(k1s,l)<= ZERO) THEN
        hn(k1s+1) = hpp(k1s+1,l)
      ENDIF

      SELECT CASE (nwlayers)
      !.....Calculate active scalar for case of a single layer.....
      CASE (1)   

        ! ... Use h(k1s,l) instead of hn(k1s) - 
        aa( 2,k1s) = hn(k1s)/twodt1
        ds(   k1s) = ex(k1s,l) + HeatSource(k1s,l) 
        sal(k1s,l) = ds(k1s  )/aa(2,k1s)

        ! ... For one-layer columns that become dry - The
        !     value of the threshold 1.E-5 is completely 
        !     arbitrary but small - it is used to avoid the
        !     occurrence of errors in scalar conc. etimates
        !     arising from errors in dividing ds by aa - 
        !     Note that the error allowed in estimating zeta
        !     is of O(10-6) - see SUB. SOLVER
        IF (h(k1s,l) < 1.E-2) sal(k1s,l) = salpp(k1s,l) 

      !.....Calculate active scalar for case of two or more layers.....
      CASE (2:)

         !.....Form coefficient matrix [aa]
         ! Define upper diagonal terms
         aa(3,k1s:kms-1) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(3,kms)       =  0.0
         ! Define lower diagonal terms
         aa(1,k1s+1:kms) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(1,k1s)       =  0.0
         ! Define center diagonal terms
         aa(2,k1s:kms)   =  hn(k1s:kms)/twodt1-aa(1,k1s:kms)-aa(3,k1s:kms)
     
         !.....form r.h.s. matrix [ds]..... 
         DO k = k1s, kms
            ds(k) = ex(k,l) + HeatSource(k,l) ! ex -> temperature term
         ENDDO

         ! ... Modify matrices of temperature term to take into accout mixing action of sinks/sources
         IF ( iopss > 0 ) THEN

             DO inn = 1, iopss

                 ! ... Temperature modif for adjacent cell to the jet cell - Chris
                 IF (ptype(iodev(inn)) == 0 .AND. inn /= iopss) THEN                                 !  add heat source to jet only when it's jet model
                     nn  = iodev(inn);
                     IF  ( (l == lEC(ij2l(ipss(inn),jpss(inn))) .AND. uEpss(nn) > 0) .OR.  &         !  east exit of the jet cell
                         & (l == lWC(ij2l(ipss(inn),jpss(inn))) .AND. uWpss(nn) > 0) .OR.  &
                         & (l == lNC(ij2l(ipss(inn),jpss(inn))) .AND. vNpss(nn) > 0) .OR.  &
                         & (l == lSC(ij2l(ipss(inn),jpss(inn))) .AND. vSpss(nn) > 0)       ) THEN

                         DO k = k1s, kms
                             ds(k)=ds(k) + Qpss(k,inn)*Tpss(k,inn)/(dx*dy)*(1-1/(1+sineagl(nn)))     !  diffusion angle correction same coeff applied on wp
                         ENDDO
                     ENDIF
                 ENDIF


                 IF ( j /= jpss(inn) .OR. i /=ipss(inn) ) CYCLE
                 SELECT CASE (ptype(iodev(inn)))
                     ! ... account the mixing of jet model
                     CASE (0 )
                         DO k = k1s, kms
                             ! control the overall temperature diffusion (temperature source term) Chris
                             Qsource  = Qpss(k,inn)/(dx*dy)                                                       &
                                 & *(uWpss(nn)+uEpss(nn)+vNpss(nn)+vSpss(nn))
                             ! Qsource  = 0.0;
!                             IF  ((uEpss(nn) > 0) .OR. (uWpss(nn) > 0)) THEN
!                                 Tsource  = (salpp(k+1,l)+salpp(k-1,l)+salpp(k,lNC(l))+salpp(k,lSC(l)))/4! Tpss(k,inn)
!                             ELSE
!                                 Tsource  = (salpp(k+1,l)+salpp(k-1,l)+salpp(k,lEC(l))+salpp(k,lWC(l)))/4
!                             ENDIF
                             Tsource = Tpss(k,inn);
                             ds(k)    = ds(k)+Qsource*Tsource*SQRT(1-(2*sineagl(nn))**2)             ! offset the temperature source
                         ENDDO

                     ! ... account the mixing of plume
                     CASE (1:)
                         DO k = k1s, kms
                             IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                             Qsource  = Qpss(k,inn)/(dx*dy)
                             Tsource  = Tpss(k,inn)
                             ds(k)    = ds(k)+Qsource*Tsource
                         ENDDO
                 END SELECT
             ENDDO

         ENDIF

         !.....Solve tridiagonal system for the
         !     vertical distribution of active scalar.....
         CALL trid1 (aa, ds, sal1, k1s, kms, km1, nwlayers)
 
         !.....Define scalars at new time step....
         sal(k1s:kms  ,l) = sal1(1:nwlayers)
         sal(k1 :k1s-1,l) = sal1(1         ) 

         ! if it is the source and it is the water jet model then correct the source temperature - Chris 07/MAY/2015
         DO inn = 1, iopss
             IF ( j == jpss(inn) .AND. i ==ipss(inn) .AND. ptype(iodev(inn)) == 99 .AND. SUM(Qpss(:,inn))>0 ) THEN ! DISABLED
                 IF  ((uEpss(nn) > 0) .OR. (uWpss(nn) > 0)) THEN
                     sal(kms-jetcell-jettopc+1:kms-jetcell,l)= (salpp(kms-jetcell-jettopc+1:kms-jetcell,lEC(l)) &
                         & +salpp(kms-jetcell-jettopc+1:kms-jetcell,lWC(l)))/2;
                 ELSE
                     sal(kms-jetcell-jettopc+1:kms-jetcell,l)= (salpp(kms-jetcell-jettopc+1:kms-jetcell,lNC(l)) &
                         & +salpp(kms-jetcell-jettopc+1:kms-jetcell,lSC(l)))/2;
                 ENDIF
             ENDIF
         ENDDO

        
      END SELECT

   !.....End loop over scalar-pts.....
   END DO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_salin = t_salin + (etime - btime)

END SUBROUTINE imsal

!***********************************************************************
PURE FUNCTION densty_s ( temperature, salinity )
!***********************************************************************
!
!  Purpose: To compute density (in kg/m**3) from active scalars
!           It uses UNESCO Eq.of state for density of freshwater 
!           taken from Gill(1982) - Atmosphere-Ocean Dynamics, Appendix 3
!           However, at this point pressure (depth) effects are not
!           included in the calculation of water density.
!           This function is based on the original function written by 
!           P.E. Smith in which the first arg. was salinity and the 2nd temp.
!           Here, we use temp. as first argument, as it is the first arg. whose
!           transport equation is solved in the code. These changes
!           were made as temp. is the main active scalar in Stockton Channel.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

    ! ... Io variables
	REAL, INTENT(IN) :: temperature, salinity
	REAL             :: densty_s
	
    densty_s =999.842594                &
      +6.793952e-2*temperature          &
      -9.095290e-3*temperature**2.      & 
      +1.001685e-4*temperature**3.      &
      -1.120083e-6*temperature**4.      &
      +6.536332e-9*temperature**5.
  
END FUNCTION densty_s

!***********************************************************************
SUBROUTINE trid ( acoef, g, r, ag, ar, k1, km, km1, n )
!***********************************************************************
!            
!  Purpose: Tridiagonal matrix solver for the momentum equation using
!           the double-sweep method
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   7/1/98           P.E. Smith        Original f90 code
!  2/15/99           P.E. Smith        Reset kmax from 20 to 60 layers
!  6/15/99           P.E. Smith        Reset kmax to 200 layers
!
!-----------------------------------------------------------------------

   !.....Dimensioning parameter.....
   INTEGER, PARAMETER :: kmax = 500

   !.....Arguments.....
   INTEGER, INTENT(IN) :: k1, km, km1, n
   REAL, DIMENSION(km1), INTENT(IN)    :: g, r
   REAL, DIMENSION(km1), INTENT(INOUT) :: ag, ar
   REAL, DIMENSION(3,km1), INTENT(IN)  :: acoef

   !.....Local variables.....
   INTEGER, AUTOMATIC :: k, kk
   REAL, AUTOMATIC, DIMENSION(kmax) :: a, b, c, d, e, c1, d1, e1, d2, e2

   !.....Timing.....
!   REAL, EXTERNAL :: TIMER
!   REAL :: btime, etime
!   btime = TIMER(0.0)
    n_trid = n_trid + 1

   !.....Load diagonals of coefficient matrix into
   !     1-d arrays and define r.h.s. vectors.....
   k = 0
   DO kk = k1, km
      k = k + 1
      a(k) = acoef(1,kk)
      b(k) = acoef(2,kk)
      c(k) = acoef(3,kk)
      d(k) = g(kk)
      e(k) = r(kk)
   END DO

   !.....Forward sweep--transform coefficient
   !     matrix into upper bidiagonal form.....
   c1(1) = c(1)/b(1)
   d1(1) = d(1)/b(1)
   e1(1) = e(1)/b(1)
   DO k = 2, n
      c1(k) = c(k)/(b(k) - a(k)*c1(k-1))
      d1(k) = (d(k) - a(k)*d1(k-1))/(b(k) - a(k)*c1(k-1))
      e1(k) = (e(k) - a(k)*e1(k-1))/(b(k) - a(k)*c1(k-1))
   END DO

   !.....Backward sweep--transform coefficient
   !     matrix into diagonal form.....
   d2(n) = d1(n)
   e2(n) = e1(n)
   DO k = n-1, 1, -1
      d2(k) = d1(k) - c1(k)*d2(k+1)
      e2(k) = e1(k) - c1(k)*e2(k+1)
   END DO

   !.....Load r.h.s. solution vectors into  [ag]  and  [ar]
   !     arrays for passing back to the calling program.....
   k = 0
   DO kk = k1, km
      k = k + 1
      ag(kk) = d2(k)
      ar(kk) = e2(k)
   END DO

END SUBROUTINE trid

!***********************************************************************
SUBROUTINE trid1 ( acoef, ds, sal, k1, km, km1, n )
!***********************************************************************
!
!  Purpose: Tridiagonal matrix solver for the salinity equation using
!           the double-sweep method
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   7/1/98           P.E. Smith        Original f90 code
!   5/1/99           P.E. Smith        Reset kmax from 20 to 60 layers
!  6/15/99           P.E. Smith        Reset kmax to 200 layers
!
!-----------------------------------------------------------------------

   !.....Dimensioning parameter.....
   INTEGER, PARAMETER :: kmax = 500

   !.....Arguments.....
   INTEGER, INTENT(IN) :: k1, km, km1, n
   REAL, DIMENSION(km), INTENT(INOUT)  :: ds
   REAL, DIMENSION(n),  INTENT(INOUT)  :: sal
   REAL, DIMENSION(3,km1), INTENT(IN)  :: acoef

   !.....Local variables.....
   INTEGER :: k, kk
   REAL, DIMENSION(kmax) :: a, b, c, d, e, f

   ! Note: n = number of 3-d layers (also number of unknowns)

   !.....Load diagonals of coefficient matrix into
   !     1-d arrays and rename r.h.s. vector.....
   k = 0
   DO kk = k1, km
      k = k + 1
      a(k) = acoef(1,kk)
      b(k) = acoef(2,kk)
      c(k) = acoef(3,kk)
      d(k) = ds(kk)
   END DO

   !.....Initialize for forward sweep.....
      e(1) = -c(1)/b(1)
      f(1) =  d(1)/b(1)

   !.....Forward sweep (solve for  e  and  f  vectors).....
   IF(n == 2) GO TO 1
   DO k = 2, n-1
      e(k) = -c(k)/(b(k)+a(k)*e(k-1))
      f(k) = (d(k)-a(k)*f(k-1))/(b(k)+a(k)*e(k-1))
   END DO

   !.....Compute salinity in bottom layer.....
 1 sal(n) = (d(n)-a(n)*f(n-1))/(b(n)+a(n)*e(n-1))

   !.....Backward sweep (solve for salinity vector).....
   DO k = n-1, 1, -1
      sal(k) = e(k)*sal(k+1) + f(k)
   END DO

END SUBROUTINE trid1

!***********************************************************************
SUBROUTINE init
!***********************************************************************
!
!  Purpose: To define initial conditions for the simulation.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios,     &
              kmx  , kmy  , kms,   &
              nwlsp, nwlup, nwlvp, js
   REAL :: x, rhoz, salz, zeta0, amp, deltZ
   INTEGER, PARAMETER :: InitProc = 0

   ! ... Initialize time step counter, time in seconds and hours from &
   !     start of simulations (these are global variables defined in 
   !     si3d_types.
   n = 0; its = 0; thrs = 0.0E0

   !.....Define initial water surface elevations.....
   SELECT CASE (InitProc)
   CASE (0)
     s    = zetainit; 
     sp   = zetainit; 
     spp  = zetainit;
     !s  (1:46,:) = -15.;
     !sp (1:46,:) = -15.;
     !spp(1:46,:) = -15.;
 
   CASE(1)

     amp = 0.05;
     DO i = 1, im1
       ! x = REAL(i*idx) - 1.5*dx ! idt real
       x = i*idx - 1.5*dx ! idt real
       zeta0 = amp*COS(pi*x/xl)
       deltZ = zlevel(k1+1); ! Move initial free surface to first interface
       DO j = 1, jm1
         sp  (i,j) = zeta0  - deltZ; 
         spp (i,j) = zeta0  - deltZ;
         s   (i,j) = zeta0  - deltZ;
       END DO
     END DO
   END SELECT

   ! ... Define thickness of cells a s-,u- & v-points
   hup = ZERO; 
   hvp = ZERO; 
   hp  = ZERO; 
   k1z = km1; 
   k1u = km1; 
   k1v = km1;
   DO l = 1, lm

      ! ... Map 3D-(i,j) from 2D-l indexes
      i = l2i(l); j = l2j(l);
      
      ! ... At zeta-points
      kms = kmz(i,j)
      nwlsp = 0 
      DO k = k1, kms
          hp (k,l)=AMIN1(zlevel(k+1),hhs(i,j)) -            &
          &        AMAX1(zlevel(  k),-sp(i,j))
          IF(hp(k,l) > HMIN) THEN
            nwlsp = nwlsp + 1;
            IF(nwlsp==1) k1z(i,j) = k
          ELSE
            hp(k,l)=ZERO;
          ENDIF
      ENDDO

      ! Set zeta = hhs(i,j) for columns with mask2d = TRUE (i.e.
      ! potentially wett) but intitially dry (k1z = km1).
      IF (k1z(i,j) == km1) THEN
         s  (i,j) = -hhs(i,j)+HMIN;
         sp (i,j) = -hhs(i,j)+HMIN;
         spp(i,j) = -hhs(i,j)+HMIN;
      ENDIF
  
      ! ... At u-points    
      IF (mask2d(i+1,j)) THEN
        kmx = MIN(kmz(i,j),kmz(i+1,j))
        nwlup = 0
        DO k = k1, kmx
          hup(k,l)=AMIN1(zlevel(k+1),hhu(i,j)) -            &
          &        AMAX1(zlevel(  k),-(sp(i,j)+sp(i+1,j))/2.)
          IF(hup(k,l) > HMIN) THEN
            nwlup = nwlup + 1;
            IF(nwlup==1) k1u(i,j) = k
          ELSE
            hup(k,l)=ZERO;
          ENDIF
        ENDDO
      ENDIF
      
      ! ... At v-points    
      IF (mask2d(i,j+1)) THEN
        kmy = MIN(kmz(i,j),kmz(i,j+1))
        nwlvp = 0
        DO k = k1, kmy
          hvp(k,l)=AMIN1(zlevel(k+1),hhv(i,j)) -            &
          &        AMAX1(zlevel(  k),-(sp(i,j)+sp(i,j+1))/2.)
          IF(hvp(k,l) > HMIN) THEN
            nwlvp = nwlvp + 1;
            IF(nwlvp==1) k1v(i,j) = k
          ELSE
            hvp(k,l)=ZERO;
          ENDIF
        ENDDO
      ENDIF
 
   ENDDO  

   hupp = hup; hu = hup; 
   hvpp = hvp; hv = hvp; 
   hpp  = hp ; h  = hp

   !.....Check if the cells hosting diffusers have other cells on both sides
   !.....Check if the cells hosting diffusers have other cells on both sides
   !IF (iopss > 0) THEN
   !  DO js = 1, iopss
   !     i = ipss(js)
   !     j = jpss(js)
   !     l = ij2l(i,j)
   !     k = kmz(i,j)-1 
   !     IF(ptype(js)==1) THEN
   !       IF(hup(k,l)<=ZERO .OR. hup(k,lWC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !         PRINT *, kmz(i,j), kmz(i+1,j), kmz(i-1,j)
   !         STOP
   !       ENDIF
   !     ELSEIF(ptype(js)==2) THEN
   !       IF(hvp(k,l)<=ZERO .OR. hvp(k,lSC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !         PRINT *, kmz(i,j), kmz(i,j+1), kmz(i,j-1)
   !         STOP
   !       ENDIF
   !     ELSEIF(ptype(js)==3) THEN
   !       IF(hup(k,l)<=ZERO .OR. hup(k,lWC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !         PRINT *, kmz(i,j), kmz(i+1,j), kmz(i-1,j)
   !         STOP
   !       ENDIF
   !       IF(hvp(k,l)<=ZERO .OR. hvp(k,lSC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !         PRINT *, kmz(i,j), kmz(i,j+1), kmz(i,j-1)
   !         STOP
   !       ENDIF
   !     ENDIF
   !  ENDDO
   !ENDIF

   !IF (iopss > 0) THEN
   !  DO js = 1, iopss
   !     i = ipss(js)
   !     j = jpss(js)
   !     l = ij2l(i,j)
   !     k = kmz(i,j)-3 
   !     IF(ptype(js)==1) THEN
   !       IF(hup(k,l)<=ZERO .OR. hup(k,lWC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !         PRINT *, kmz(i,j), kmz(i+1,j), kmz(i-1,j)
   !       ENDIF
   !     ENDIF
   !     IF(ptype(js)==2) THEN
   !       IF(hvp(k,l)<=ZERO .OR. hvp(k,lSC(l))<=ZERO) THEN
   !         PRINT *, 'Cell (i,j)=(', i,',', j,') has one of the sides closed'
   !       ENDIF
   !     ENDIF
   !  ENDDO
   !ENDIF
 
   !.....Initialize 1-d arrays.....
   uout  = 0.0; 
   vout  = 0.0; 
   wout  = 0.0; 
   uhout = 0.0
   Avout = 0.0; 
   Dvout = 0.0; 
   sal1 = 0.0; 
   ds = 0.0;  
 
   !.....Initialize solution arrays.....
   eagx = 0.0; 
   eagy = 0.0; 
   earx = 0.0; 
   eary = 0.0
   sx = 0.0; 
   sy = 0.0; 
   dd = 0.0; 
   qq = 0.0; 
   rr = 1.0; 

   !.....Initialize velocity arrays.....
   u=u0; up=u0; upp=u0; uh=u0h0; uhp=u0h0; uhpp=u0h0
   v=v0; vp=v0; vpp=v0; vh=v0h0; vhp=v0h0; vhpp=v0h0 
   wp=w0

   ! ... Initialize eddy coefficient arrays ..... 
   Av=0.0; 
   Dv=0.0; 
   
   !.....Initialize arrays used in the soln of the matrix mom eq.....
   ex=0.0; agx = 0.0; arx = 0.0; agy = 0.0; ary = 0.0 
   
   !.....Initialize 3D active scalar and density arrays.....
   CALL InitializeScalarFields

   ! ... Inialize turbulence quanties for Turbulence Model
   CALL InitializeTurbulenceModel

   ! ... Initialize variables used in Point Source-Sint models
   IF (iopss > 0) THEN
       Qpss = 0.0E0
       Tpss = 0.0E0
       Rpss = 0.0E0
   ENDIF

   ! .... Initialize other models
   SELECT CASE (ecomod)
       CASE (-1) ! Tracer Cloud Releases
           CALL InitTracerCloud
       CASE (1) ! Water quality routines
           CALL WQinit
       CASE (4) ! FABM
           PRINT *, " FABM initialization..."
           CALL init_fabm
   END SELECT

   ! ... Deallocate vector holding scalar concs.
   DEALLOCATE(ScalarProfile)


END SUBROUTINE init

!************************************************************************
SUBROUTINE InitializeScalarFields
!************************************************************************
!
!  Purpose: To read the salinity initial condition from a file. (Note:
!           The salinity initial condition file will normally have been
!           prepared by a pre-processing program.)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios, imm1, jmm1, kmm1, ncols, ncols1, nc, &
              nsets, ia, ib, nn, ntr1
   INTEGER :: nprof, npf, ils, ile, jls, jle, nci ! Two-basin initalization
   REAL    :: Vamp, rhoamp, Ts, Tb,   &  ! Used to initialize IW-problem
              NBV, meandepth, length, &
              rhohere, x, z, rhos, rhob
   CHARACTER(LEN=18)  :: initfmt
   INTEGER, PARAMETER :: InitProc =  0

   SELECT CASE (InitProc)
  
   ! ... OPTION -1 -  Use two profiles to initialize  ---------------------
   CASE (-1)

     !.....Open salinity initial condition file.....
     sal_ic_file = 'si3d_init.txt'
     OPEN (UNIT=i4, FILE='si3d_init.txt', STATUS="OLD", FORM="FORMATTED", IOSTAT=ios)
     IF(ios /= 0) CALL open_error ( "Error opening "//sal_ic_file, ios )
 
     ! .... Skip over first five header records in init file  
     READ (UNIT=i4, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 13 )

     ! .... Read no. of profiles used to intialize the model
     READ (UNIT = i4, FMT='(14X,I20)', IOSTAT=ios) nprof
     IF (ios /=0) CALL input_error ( ios, 13 )
     nprof = MAX(nprof,1)

     ! .... Allocate space for working variables
     ALLOCATE ( ScalarProfile (km1, rdtr+1), STAT = ios )
     IF (ios /= 0) THEN; PRINT *, 'Error alloc. init. arrays'; STOP; ENDIF
	    
     ! .....Write the format of the data records into an internal file
     WRITE (UNIT=initfmt, FMT='("(10X,",I3,"G11.2)")') rdtr+1

     ! ... Loop over profiles for different regions in initial field
     nci  = 0; 
     salp = 0.0E0; IF (rdtr>1) tracer = 0.0E0;
     DO npf = 1, nprof

       ! ... Read indexes defining limits of sub-domain
       READ (UNIT=i4,FMT='(/(14X,I20))',IOSTAT=ios) ils
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) ile
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) jls
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) jle
       READ (UNIT=i4,FMT='(/        ) ',IOSTAT=ios) 
 
       ! ... Read data array 
       DO k = 1, km1
         READ (UNIT=i4, FMT=initfmt, IOSTAT=ios) &
              (ScalarProfile(k,nn), nn = 1, rdtr+1)
         IF (ios /= 0) CALL input_error ( ios, 14 )
       END DO

       ! ... Loop over columns within the sub-domain
       DO i = i1, im; DO j = j1, jm;
         IF (.NOT.mask2d(i,j)) CYCLE
         IF ( (i .GE. ils) .AND. &
         &    (i .LT. ile) .AND. &
         &    (j .GE. jls) .AND. &
         &    (j .LT. jle) ) THEN 		   

           ! .... map (i,j) into l-index & count
           l = ij2l(i,j); nci = nci + 1;

           ! .... Initialize active scalar fields 
           DO k = 1, km1; 
             salp(k,l) = ScalarProfile(k,1)   
           END DO;

           ! ... Initialize Non-active scalar fields 
           IF (rdtr > 0) THEN
             DO  nn = 1, rdtr
               DO k = 1, km1
                 tracer(k,l,nn) = ScalarProfile(k,nn+1)  
               ENDDO
             ENDDO ! ... End loop over tracers	  
           ENDIF
         ENDIF
       ENDDO; ENDDO; 

     ENDDO
     ! ... Check if all the wett domain has been initialized
     IF (nci .NE. lm) THEN
	   PRINT *, 'ERROR: Some columns have not been initialized'
	   STOP
     ENDIF
     ! ... Initialize variables at other time steps
     sal = salp; 
     salpp = salp;
     tracerpp = tracer;

     PRINT *, '..... Done initializing scalar fields from ASCII file'

     ! ... Close io file
     CLOSE (i4)
 

   ! ... OPTION 0 -  Initialize from file ---------------------
   CASE (0)

     !.....Open salinity initial condition file.....
     sal_ic_file = 'si3d_init.txt'
     OPEN (UNIT=i4, FILE='si3d_init.txt', STATUS="OLD", FORM="FORMATTED", IOSTAT=ios)
     IF(ios /= 0) CALL open_error ( "Error opening "//sal_ic_file, ios )
 
     ALLOCATE ( ScalarProfile (km1, rdtr+1), STAT = ios )
     IF (ios /= 0) THEN; PRINT *, 'Error alloc. init. arrays'; STOP; ENDIF

     ! Skip over first five header records in open boundary condition file  
     READ (UNIT=i4, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 13 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=initfmt, FMT='("(10X,",I3,"G11.2)")') rdtr+1

     ! Read data array and store it in memory
     ntr1 = rdtr; IF (ecomod < 0) rdtr = 0;
     DO k = 1, km1
       READ (UNIT=i4, FMT=initfmt, IOSTAT=ios) &
            (ScalarProfile(k,nn), nn = 1, ntr1+1)
       IF (ios /= 0) CALL input_error ( ios, 14 )
     END DO

     ! ... Initialize the active scalar field (allways) 
     salp = 0.0
     DO k = 1, km1; 
        salp(k,:) = ScalarProfile(k,1)
     END DO;
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (rdtr > 0) THEN
       tracer = 0.0;
       IF (ecomod < 0 ) THEN
         CALL InitTracerCloud
       ELSE
         DO  nn = 1, rdtr
           DO k = 1, km1
             tracer(k,:,nn) = ScalarProfile(k,nn+1)  
           ENDDO
         END DO ! ... End loop over tracers	  
       ENDIF
       tracerpp = tracer;
     ENDIF

     ! ... Close io file
     CLOSE (i4)
 
   ! ... OPTION 1 - Use analytical solution to excite IW ----------
   CASE (1) 

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Parameters used to define the solution
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom 
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! FLOAT(km-k1+1)*ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (rdtr > 0) THEN
       DO nn = 1, rdtr;
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   ! ... OPTION 2 - Use analytical solution in half a closed basin to test 
   !                the nesting algorithms nesting. All variables defining the basin
   !                & the IW need to be the same in the fine & coarse grid - 
   !                In the fine grid we only modify the length and x. 
   CASE (2) 

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Make sure these constants are as in CASE (1)
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx; length = length * 2.; 
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom 
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx; x = x + length/2.; 
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! z = FLOAT(k) * ddz - 1.5 * ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (rdtr > 0) THEN
       DO nn = 1, rdtr;
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   ! ... OPTION 2 - Use analytical solution in half a closed basin to test 
   !                the nesting algorithms nesting. All variables defining the basin
   !                & the IW need to be the same in the fine & coarse grid - 
   !                In the fine grid we only modify the length and x. - EAST sims
   CASE (3) 

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Make sure these constants are as in CASE (1)
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx; length = length * 2.; 
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom 
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx;    
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! z = FLOAT(k) * ddz - 1.5 * ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (rdtr > 0) THEN
       DO nn = 1, rdtr
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF


   END SELECT

    ! tracer injection at the beginning Chris now disabled
!   IF (iinj > 0) THEN
!       CALL do_tracer_injection ()
!   ENDIF

   ! ... Initialize density field at time n-1 & n
   DO l = 1, lm1; DO k = k1, km1; 
      rhop(k,l) = densty_s ( salp(k,l), t0 ) - 1000.
   END DO; END DO

END SUBROUTINE InitializeScalarFields

!************************************************************************
SUBROUTINE InitializeScalarFieldsTWO
!************************************************************************
!
!  Purpose: To read the salinity initial condition from a file. (Note:
!           The salinity initial condition file will normally have been
!           prepared by a pre-processing program.)
!
!-------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios, imm1, jmm1, kmm1, ncols, ncols1, nc, &
              nsets, ia, ib, nn
   INTEGER :: nprof, npf, ils, ile, jls, jle, nci
   REAL    :: Vamp, rhoamp, Ts, Tb,   &  ! Used to initialize IW-problem
              NBV, meandepth, length, &
              rhohere, x, z, rhos, rhob
   CHARACTER(LEN=18)  :: initfmt
   INTEGER, PARAMETER :: InitProc = 0  

   SELECT CASE (InitProc)
  
   ! ... OPTION 0 -  Initialize from file ---------------------
   CASE (0)

     !.....Open salinity initial condition file.....
     sal_ic_file = 'si3d_init.txt'
     OPEN (UNIT=i4, FILE='si3d_init.txt', STATUS="OLD", FORM="FORMATTED", IOSTAT=ios)
     IF(ios /= 0) CALL open_error ( "Error opening "//sal_ic_file, ios )
 
     ! .... Skip over first five header records in init file  
     READ (UNIT=i4, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 13 )

     ! .... Read no. of profiles used to intialize the model
     READ (UNIT = i4, FMT='(14X,I20)', IOSTAT=ios) nprof
     IF (ios /=0) CALL input_error ( ios, 13 )
     nprof = MAX(nprof,1)

     ! .... Allocate space for working variables
     ALLOCATE ( ScalarProfile (km1, ntr+1), STAT = ios )
     IF (ios /= 0) THEN; PRINT *, 'Error alloc. init. arrays'; STOP; ENDIF
	    
     ! .....Write the format of the data records into an internal file
     WRITE (UNIT=initfmt, FMT='("(10X,",I3,"G11.2)")') ntr+1

     ! ... Loop over profiles for different regions in initial field
     nci  = 0; 
     salp = 0.0E0; IF (ntr>1) tracer = 0.0E0;
     DO npf = 1, nprof

       ! ... Read indexes defining limits of sub-domain
       READ (UNIT=i4,FMT='(/(14X,I20))',IOSTAT=ios) ils
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) ile
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) jls
       READ (UNIT=i4,FMT='(  14X,I20) ',IOSTAT=ios) jle
       READ (UNIT=i4,FMT='(/        ) ',IOSTAT=ios) 
 
       ! ... Read data array 
       DO k = 1, km1
         READ (UNIT=i4, FMT=initfmt, IOSTAT=ios) &
              (ScalarProfile(k,nn), nn = 1, ntr+1)
         IF (ios /= 0) CALL input_error ( ios, 14 )
       END DO

       ! ... Loop over columns within the sub-domain
       DO i = i1, im; DO j = j1, jm;
         IF (.NOT.mask2d(i,j)) CYCLE
         IF ( (i .GE. ils) .AND. &
         &    (i .LT. ile) .AND. &
         &    (j .GE. jls) .AND. &
         &    (j .LT. jle) ) THEN 		   

           ! .... map (i,j) into l-index & count
           l = ij2l(i,j); nci = nci + 1;

           ! .... Initialize active scalar fields 
           DO k = 1, km1; 
             salp(k,l) = ScalarProfile(k,1)   
           END DO;

           ! ... Initialize Non-active scalar fields 
           IF (ntr > 0) THEN
             DO  nn = 1, ntr
               DO k = 1, km1
                 tracer(k,l,nn) = ScalarProfile(k,nn+1)  
               ENDDO
             ENDDO ! ... End loop over tracers	  
           ENDIF
         ENDIF
       ENDDO; ENDDO; 

     ENDDO
     ! ... Check if all the wett domain has been initialized
     IF (nci .NE. lm) THEN
	   PRINT *, 'ERROR: Some columns have not been initialized'
	   STOP
     ENDIF
     ! ... Initialize variables at other time steps
     sal = salp; 
     salpp = salp;
     tracerpp = tracer;

     PRINT *, '..... Done initializing scalar fields from ASCII file'

     ! ... Close io file
     CLOSE (i4)
 
   ! ... OPTION 1 - Use analytical solution to excite IW ----------
   CASE (1) 

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Parameters used to define the solution
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom 
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! FLOAT(km-k1+1)*ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (ntr > 0) THEN
       DO nn = 1, ntr; 
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   ! ... OPTION 2 - Use analytical solution in half a closed basin to test 
   !                the nesting algorithms nesting. All variables defining the basin
   !                & the IW need to be the same in the fine & coarse grid - 
   !                In the fine grid we only modify the length and x. 
   CASE (2) 

     Vamp = 0.10; Ts = 25.00; Tb = 15.0; ! Make sure these constants are as in CASE (1)
     meandepth = zl; ! FLOAT(km-k1+1)*ddz
     length    = FLOAT(im-i1+1)*dx; length = length * 2.; 
     rhos = 1028.*(1.-1.7E-4*(Ts-10.));		! Surface density
     rhob = 1028.*(1.-1.7E-4*(Tb-10.));		! Bottom  density
     drho = rhos - rhob;		! Change in density from top to bottom 
     NBV=SQRT(-g/rhos*drho/meandepth);
     rhoamp=rhos*Vamp*NBV/g;
     DO l = 1, lm
       i = l2i(l); j = l2j(l);
       x = FLOAT(i) * dx - 1.5 * dx; x = x + length/2.; 
       DO k = k1, km
          z = zlevel(k+1) - 0.5 * hp(k,l); ! z = FLOAT(k) * ddz - 1.5 * ddz
          rhohere = rhos -z*drho/meandepth+rhoamp*COS(pi*x/length)*SIN(pi*z/meandepth);
          salp(k,l)= 10.-((rhohere-1028.)/1028.)/1.7E-4;
       ENDDO
     END DO
     sal = salp; 
     salpp = salp;

     ! ... Initialize Non-active scalar fields 
     IF (ntr > 0) THEN
       DO nn = 1, ntr; 
         tracer(:,:,nn) = sal;
       ENDDO
       tracerpp = tracer;
     ENDIF

   END SELECT

   ! ... Initialize density field at time n-1 & n
   DO l = 1, lm1; DO k = k1, km1; 
      rhop(k,l) = densty_s ( salp(k,l), t0 ) - 1000.
   END DO; END DO

END SUBROUTINE InitializeScalarFieldsTWO

!***********************************************************************
SUBROUTINE smooth
!***********************************************************************
!
!  Purpose: To smooth the solution from the leapfrog step with the
!           Asselin time filter (Mon. Weather Rev., v. 100, 1972,
!           p. 487-490). Smoothing is only performed if ismooth>=1.
!           The degree of smoothing is determined from the parameter
!           beta. Beta=0.05 is recommended. Values as high as 1.0
!           can be used. The choices for ismooth are:
!                If ismooth = 0 --> no smoothing
!                If ismooth = 1 --> smooth zeta and velocity
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kmx, kmy, kms, k1s, k1x, k1y
   REAL    :: wght, wghtpp, scC, scCpp

   !.....Smooth zeta (and then recalculate hp, hup, and hvp).....
   DO l = 1, lm; 
     i = l2i(l); j = l2j(l);
     sp(i,j) = sp(i,j) + (beta2)*(s(i,j)-2.*sp(i,j)+spp(i,j))
   ENDDO
   CALL layer_hp 

   !.....Smooth horizontal velocity components.....
   DO l = 1, lm; 
     ! ... Map l- into (i,j)-indexes
     i = l2i(l); j = l2j(l);
     ! ... At u-points     
     IF(mask2d(i+1,j)) THEN
       kmx = MIN(kmz(i,j),kmz(i+1,j))
       DO k = k1, kmx
         IF(hup(k,l)<=0.0) CYCLE
         uhp(k,l)= uhp(k,l)+beta2*(uh(k,l)-2.*uhp(k,l)+uhpp(k,l))
         up (k,l)= uhp(k,l)/hup(k,l)
       ENDDO
     ENDIF
     ! ... At v-points
     IF(mask2d(i,j+1)) THEN
       kmy = MIN(kmz(i,j),kmz(i,j+1))
       DO k = k1, kmy
         IF(hvp(k,l)<=0.0) CYCLE
         vhp(k,l)= vhp(k,l)+beta2*(vh(k,l)-2.*vhp(k,l)+vhpp(k,l))
         vp (k,l)= vhp(k,l)/hvp(k,l)
       ENDDO
     ENDIF
   ENDDO

   !.....No need to recalculate vertical velocity components
   !     since these values should be stored in wpp either in save or 
   !     in settrap, which are not used in any computations - 

END SUBROUTINE smooth

!***********************************************************************
SUBROUTINE settrap
!***********************************************************************
!
!  Purpose: To setup the arrays at the  n  and  n+1/2  time levels
!           for use in the first iteration of the trapezoidal step.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables..... 
   INTEGER :: i, j, k, l, kmx, kmy, kms, k1x, k1y, k1s
   REAL    :: uutemp, vvtemp, wght, wghtpp, scC, scCpp
   
   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)
   

   !....Zeta array.....
   spp  = sp
   sp   = 0.5*(s + spp)

   ! ... Save layer thickness at time n 
   hpp  = hp; 
   hupp = hup; 
   hvpp = hvp; 

   ! ... Define layer thickness at time n+1/2 & 
   !     recompute top layer index
   CALL layer_hp
   CALL TopLayerIndexp

   ! ... Define variable values at time n+1/2. 
   !     Only define values at cells that at n+1/2 
   !     are wett. 
   DO l = 1, lm

     ! ... Map 3D-(i,j) from 2D-l indexes
     i = l2i(l); j = l2j(l);
      
     ! ... At s-points
     kms = kmz(i,j) 
     DO k = k1, kms
       salpp(k,l) = salp(k,l);     
       salp (k,l)=(sal(k,l)+salpp(k,l))/2.
       rhop (k,l)=densty_s(salp(k,l),t0)-1000.
     ENDDO
   
     ! ... At u-points    
     IF (mask2d(i+1,j)) THEN
       kmx = MIN(kmz(i,j),kmz(i+1,j))
       DO k = k1, kmx
         uhpp(k,l) = uhp(k,l)
         upp (k,l) = up(k,l)
         IF (hup(k,l)>ZERO) THEN
           uhp (k,l) = 0.5*(uh(k,l) + uhpp(k,l))
           up  (k,l) = uhp(k,l)/hup(k,l) 
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs. at n+1/2
       k = k1u(i,j);
	  
       ! a. Wetting occurs from n+1/2 to n+1
       IF (hu  (k-1,l) > ZERO) THEN
         uhp(k,l) = uhp(k,l)+uh(k-1,l)/2.
         up (k,l) = uhp(k,l) / hup(k,l)
       ENDIF
	 
       ! b. Drying occurs from n to n+1/2
       IF (hupp(k-1,l) > ZERO) THEN
         uhp (k  ,l) = uhp (k,l)+uhpp(k-1,l)/2.
         up  (k  ,l) = uhp(k,l) / hup(k,l)
       ENDIF	 	    	  

     ENDIF
      
     ! ... At v-points    
     IF (mask2d(i,j+1)) THEN
       kmy = MIN(kmz(i,j),kmz(i,j+1))
       DO k = k1, kmy
         vhpp(k,l) = vhp(k,l)
         vpp (k,l) = vp(k,l)
         IF (hvp(k,l)>ZERO) THEN
           vhp (k,l) = 0.5*(vh(k,l) + vhpp(k,l))
           vp  (k,l) = vhp(k,l)/hvp(k,l)   
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0            
         ENDIF	          
       ENDDO

       ! ... Redo near surface flux calcs. at n+1/2 
       k = k1v(i,j);
	  
       ! a. Wetting occurs from n+1/2 to n+1
       IF (hv  (k-1,l) > ZERO) THEN
         vhp(k,l) = vhp(k,l)+vh(k-1,l)/2.
         vp (k,l) = vhp(k,l) / hvp(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hvpp(k-1,l) > ZERO) THEN
         vhp (k,l) = vhp (k,l)+vhpp(k-1,l)/2.
         vp  (k,l) = vhp (k,l) /  hvp (k,l)
       ENDIF
	
     ENDIF

   ENDDO  
          
   !.....Recalculate vertical velocity at n+1/2 -  used in
   !     computing horizontal fluxes at n+1 
   CALL continuity(2)

   ! ... Save bndry. variables from n-1 into n
   IF (nopen > 0) THEN
    uhEBpp = uhEBp; huEBpp = huEBp;
    uhWBpp = uhWBp; huWBpp = huWBp;
    vhNBpp = vhNBp; hvNBpp = hvNBp;
    vhSBpp = vhSBp; hvSBpp = hvSBp;
    uhEBp  = (uhEBpp + uhEB)/2.; huEBp = (huEBp+huEB)/2.;
    uhWBp  = (uhWBpp + uhWB)/2.; huWBp = (huWBp+huWB)/2.;
    vhNBp  = (vhNBpp + vhNB)/2.; hvNBp = (hvNBp+hvNB)/2.;
    vhSBp  = (vhSBpp + vhSB)/2.; hvSBp = (hvSBp+hvSB)/2.;
   ENDIF

   ! ... Work with turbulence quantities (TurbModel)
   IF (iturb>0) CALL settrap_2EqTVars

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_settrap = t_settrap + (etime - btime)

END SUBROUTINE settrap

!***********************************************************************
SUBROUTINE save
!***********************************************************************
!
!  Purpose: Save solution for next time step
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, kms, k1s
   REAL    :: uutemp, vvtemp

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   SELECT CASE (istep)

   !                   -----After a trapezoidal step (istep=2)-----
   CASE (2)

     !.....Save zeta.....
     sp = s

     ! ... Save layer thicknesses (h is calculated from s in solver)
     hp  = h ; 
     hup = hu;
     hvp = hv; 

     ! ... Retrieve index for surface layer for next step
     CALL TopLayerIndexp

     ! ... Save state variables     
     DO l = 1, lm;   

       ! ... Map 2D-l into 3D-(i,j) indexes
       i = l2i(l); j = l2j(l)

       DO k = k1, km;
         ! ... At s-points  
         salp (k,l) = sal(k,l)
         rhop (k,l) = densty_s ( salp(k,l), t0 ) - 1000.
       ENDDO

       DO k = k1, km;
         ! ... At u-points
         IF (hup(k,l)>ZERO)THEN
           uhp (k,l) = uh (k,l) ! For horiz. scalar & momentum advection 
           up  (k,l) = uhp(k,l)/hup(k,l) ! For horiz. momentum advection
           u   (k,l) = up (k,l) ! For output purposes
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
           u   (k,l) = 0.0
         ENDIF
       ENDDO

       DO k = k1, km;
         ! ... At v-points
         IF (hvp(k,l)>ZERO)THEN
           vhp (k,l) = vh (k,l) ! For horiz. scalar & momentum advection
           vp  (k,l) = vhp(k,l)/hvp(k,l) ! For horiz. momentum advection
           v   (k,l) = vp (k,l) ! For output purposes
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
           v   (k,l) = 0.0
         ENDIF
       ENDDO

     END DO;

     !.....Recalculate vertical velocity wp to be used 
     !     in calcuation of velocities at next time step
     CALL continuity(2)

     ! ... Save bndry. variables from n-1 into n 
     IF (nopen > 0) THEN
       uhEBp  = uhEB ; huEBp = huEB ;
       uhWBp  = uhWB ; huWBp = huWB ;
       vhNBp  = vhNB ; hvNBp = hvNB ;
       vhSBp  = vhSB ; hvSBp = hvSB ;
     ENDIF

     ! ... Save Turbulence variables
     IF (iturb>0) CALL save_2EqTVars

     ! ... Save tracers
     IF (ntr > 0) THEN
      tracerpp = tracer; 
     ENDIF


   !                   -----After a leapfrog step (istep=1)-----
   CASE (1)       

     !.....Save zeta.....
     spp   = sp
     sp    = s

     !.....Save layer thickness at time n
     hpp  = hp
     hupp = hup    
     hvpp = hvp    

     ! ... Save layer thickness at time n+1
     hp  = h
     hup = hu
     hvp = hv

     ! ... Retrieve index for surface layer for next step
     CALL TopLayerIndexp

     ! .... Save other variables
     DO l = 1, lm; 
 
       ! ... Maks l- into (i,j)-indexes
       i = l2i(l); j = l2j(l);

       DO k = k1, km;
         ! ... At s-points  
         salpp(k,l) = salp(k,l)
         salp (k,l) = sal (k,l); 
         rhop (k,l) = densty_s ( salp(k,l), t0 ) - 1000.
       ENDDO

       DO k = k1, km;
         ! ... At u- and v- points
         uhpp(k,l) = uhp(k,l)
         vhpp(k,l) = vhp(k,l)
         upp (k,l) = up (k,l) 
         vpp (k,l) = vp (k,l)
       ENDDO

       DO k = k1, km;
         ! ... At u-points
         IF (hup(k,l)>ZERO)THEN
           uhp (k,l) = uh (k,l)
           up  (k,l) = uhp(k,l)/hup(k,l)
           u   (k,l) = up (k,l) ! For output purposes  
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0
           u   (k,l) = 0.0
         ENDIF
       ENDDO   

       DO k = k1, km;
         ! ... At v-points
         IF (hvp(k,l)>ZERO)THEN
           vhp (k,l) = vh (k,l)
           vp  (k,l) = vhp(k,l)/hvp(k,l)
           v   (k,l) = vp (k,l) ! For output purposes 
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0
           v   (k,l) = 0.0
         ENDIF
       END DO

     END DO

     !.....Recalculate vertical velocity wp to be used 
     !     in calcuation of velocities at next time step
     CALL continuity(2)

     ! ... Save bndry. variables  
     IF (nopen > 0) THEN
       uhEBpp = uhEBp; huEBpp = huEBp;
       uhWBpp = uhWBp; huWBpp = huWBp;
       vhNBpp = vhNBp; hvNBpp = hvNBp;
       vhSBpp = vhSBp; hvSBpp = hvSBp;
       uhEBp  = uhEB ; huEBp  = huEB ;
       uhWBp  = uhWB ; huWBp  = huWB ;
       vhNBp  = vhNB ; hvNBp  = hvNB ;
       vhSBp  = vhSB ; hvSBp  = hvSB ;
     ENDIF 

     ! ... Save Turbulence variables
     IF (iturb>0) CALL save_2EqTVars

     ! ... Save tracers
     IF (ntr > 0) THEN
      tracerpp = tracer; 
     ENDIF

   CASE DEFAULT

      PRINT *, "Invalid value of ISTEP in subroutine save"

   END SELECT

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_save = t_save + (etime - btime)

END SUBROUTINE save

!***********************************************************************
SUBROUTINE settrap2
!***********************************************************************
!
!  Purpose: To setup the arrays at the n+1/2 time level for use in
!           the second and subsequent iterations of the trapezoidal
!           step. Do not use smoothing as in the version of the code
!           by P.E. Smith.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i,j,k,l,kms,kmy,kmx
   REAL    :: wght, wghtpp

   !....Zeta array.....
   sp = 0.5*(s + spp)

   ! ...Define layers thickness at time n+1/2
   CALL layer_hp
   CALL TopLayerIndexp

   ! ... 3d arrays
   DO l = 1, lm

     ! ... Map 3D-(i,j) from 2D-l indexes
     i = l2i(l); j = l2j(l);
      
     ! ... At s-points
     kms = kmz(i,j)
     DO k = k1, kms
       salp (k,l)= (sal(k,l)+salpp(k,l))/2.
       rhop (k,l)=densty_s(salp(k,l),t0)-1000.
     ENDDO

     ! ... At u-points    
     IF (mask2d(i+1,j)) THEN
       kmx = MIN(kmz(i,j),kmz(i+1,j))
       DO k = k1, kmx
         IF (hup(k,l)>ZERO) THEN
           uhp (k,l) = 0.5*(uh(k,l) + uhpp(k,l))
           up  (k,l) = uhp(k,l)/hup(k,l) 
         ELSE
           uhp (k,l) = 0.0
           up  (k,l) = 0.0 
         ENDIF
       ENDDO

       ! ... Redo near surface flux calcs. 
       k = k1u(i,j);
	  
       ! a. Wetting occurs from n+1/2 to n+1
       IF (hu  (k-1,l) > ZERO) THEN
         uhp(k,l) = uhp(k,l)+uh(k-1,l)/2.
         up (k,l) = uhp(k,l) / hup(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hupp(k-1,l) > ZERO) THEN
         uhp (k,l) = uhp (k,l)+uhpp(k-1,l)/2.
         up  (k,l) = uhp(k,l) / hup(k,l)
       ENDIF	    	  

     ENDIF
      
     ! ... At v-points    
     IF (mask2d(i,j+1)) THEN
       kmy = MIN(kmz(i,j),kmz(i,j+1))
       DO k = k1, kmy
         IF (hvp(k,l)>ZERO) THEN
           vhp (k,l) = 0.5*(vh(k,l) + vhpp(k,l))
           vp  (k,l) = vhp(k,l)/hvp(k,l)   
         ELSE
           vhp (k,l) = 0.0
           vp  (k,l) = 0.0            
         ENDIF	          
       ENDDO

       ! ... Redo near surface flux calcs. 
       k = k1v(i,j);
	  
       ! a. Wetting occurs from n+1/2 to n+1
       IF (hv  (k-1,l) > ZERO) THEN
         vhp(k,l) = vhp(k,l)+ vh(k-1,l)/2.
         vp (k,l) = vhp(k,l) / hvp(k,l)
       ENDIF

       ! b. Drying occurs from n to n+1/2
       IF (hvpp(k-1,l) > ZERO) THEN
         vhp (k,l) = vhp (k,l)+vhpp(k-1,l)/2.
         vp  (k,l) = vhp (k,l) /  hvp (k,l)
       ENDIF
	  
     ENDIF

   ENDDO  

   !.....Recalculate vertical velocity components to be used in 
   !     calculating horizontal velocity at next iteration
   CALL continuity(2)
  
   ! ... Work with turbulence quantities (TurbModel)
   IF (iturb>0) CALL settrap2_2EqTVars

END SUBROUTINE settrap2

!***********************************************************************
SUBROUTINE outr
!***********************************************************************
!
!  Purpose: To output model run parameters & performance measures to file.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER(LEN=12) :: output_file="si3d_out.txt"
   CHARACTER :: date*8, time*10, zone*5
   INTEGER, DIMENSION(8) :: values
   INTEGER :: ios

   !.....Open output file.....
   OPEN (UNIT=i6, FILE=output_file, IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//output_file, ios )

   !.....Get date and time of run.....
   CALL date_and_time ( date, time, zone, values )

   !.....Output run parameters.....
   WRITE (UNIT=i6, FMT='(A)') title
   WRITE (UNIT=i6, FMT='("Run number = ",A8,A4,",  Start date of run:  ",I2, &
                       & "/",I2,"/",I4," at ",I4.4," hours")')date,time(1:4),&
                       & imon,iday,iyr,ihr

   !.....Output space-time domains, cell size & time step .................
   WRITE (UNIT=i6,FMT=1) xl,yl,zl,idx,idy,idz,tl,idt, dzmin,zetainit 

   ! ... Read parameters controlling solution algorithm .................
   WRITE (UNIT=i6,FMT=2) itrap,niter,ismooth,                            &
       & beta, iturb, Av0, Dv0, iadv, itrmom, ihd, Ax0, Ay0, f, theta,   &
       & ibc,isal, itrsca, nopen, cd, ifsurfbc, dtsurfbc, cw, wa, phi, ntr, id_tracer_out(1),id_tracer_out(2),id_tracer_out(3)

 1 FORMAT (/              &
     &' xl=   '  , F9.1,/ &
     &' yl=   '  , F9.1,/ &
     &' zl=   '  , F9.3,/ & ! idt real
     &' idx=  '  , F9.3,/ & ! idt real
     &' idy=  '  , F9.3,/ & ! idt real
     &' idz=  '  , F9.3,/ & ! idt real
     &' tl=   '  , G9.2,/ &
     &' idt=  '  , G9.2,/ & ! idt real
     &' dzmin='  , F9.2,/ &
     &' zeta0='  , F9.2 / )

 2 FORMAT(/              &
     &' itrap= ' ,  I9,/ &
     &' niter= ' ,  I9,/ &
     &' smooth=' ,  I9,/ &
     &' beta=  ' ,F9.4,/ &
     &' iturb= ' ,  I9,/ &
     &' Av0=   ' ,F9.4,/ &
     &' Dv0=   ' ,F9.4,/ &
     &' iadv=  ' ,  I9,/ &
     &' trmom= ' ,  I9,/ &
     &' ihd=   ' ,  I9,/ &
     &' Ax0=   ' ,F9.4,/ &
     &' Ay0=   ' ,F9.4,/ &
     &' f=     ' ,F9.4,/ &
     &' theta= ' ,F9.4,/ &
     &' ibc=   ' ,  I9,/ &
     &' isal=  ' ,  I9,/ &
     &' trsal= ' ,  I9,/ &
     &' nopen= ' ,  I9,/ &
     &' cd=    ' ,F9.4,/ &
     &' isbc=  ' ,  I9,/ &
     &' tsbc=  ' ,F9.2,/ &
     &' cw=    ' ,F9.2,/ &
     &' wa=    ' ,F9.2,/ &
     &' phi=   ' ,F9.2,/ &
     &' ntr=   ' ,  I9  ,/ &
     &' idtr(1)= ' ,  I9,/ &
     &' idtr(2)= ' ,  I9,/ &
     &' idtr(3)= ' ,  I9 )

END SUBROUTINE outr

!***********************************************************************
SUBROUTINE outt
!***********************************************************************
!
!  Purpose: To write output to timefile(s). A separate timefile is
!           opened for each node where output is requested.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER :: date*8, time*10, zone*5
   CHARACTER(LEN=9)  :: nodeno    ="         "
   CHARACTER(LEN=15) :: filenm    ="               "
   REAL :: qu, stidal, tdays
   INTEGER, DIMENSION(8)     :: values
   INTEGER :: nn, i, j, k, l, kkk, itdays, ios, nchar, it
   INTEGER, SAVE :: i10, i30, i60
   LOGICAL, SAVE :: first_entry = .TRUE.

   INTEGER :: kmzij, nk
   REAL    :: hhsij

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Open timefiles on first entry into the subroutine.....
   IF( first_entry ) THEN
      first_entry = .FALSE.
      DO nn = 1, nnodes

         i = inode(nn)
         j = jnode(nn)

         ! Convert node numbers to a character variable
         CALL nodech ( i, j, nodeno, nchar )

         ! Name the standard si3d timefile
         filenm = 'tf'//nodeno(1:nchar)//'.txt'

         ! Open the timefile
         i60 = i6 + nn    ! Use file numbers 61-80
         OPEN ( UNIT=i60, FILE=filenm, IOSTAT=ios )
         IF(ios /= 0) CALL open_error ( "Error opening "//filenm, ios )

         !.....Get date and time of run.....
         CALL date_and_time ( date, time, zone, values )

         !.....Output run title and column headings for standard si3d format.....
         WRITE (UNIT=i60, FMT='(A)') title
         WRITE (UNIT=i60, FMT='("Run number = ", A8, A4,                      &
                 & ",  Start date of run:  ",I2,                              &
                 & "/",I2,"/",I4," at ",I4.4," hours")') date,time(1:4),      &
                 & imon,iday,iyr,ihr
         IF ( i > 0 .AND. j > 0) THEN
           kmzij = (kmz(i,j)-k1+1); 
           hhsij =  hhs(i,j)      ; 
         ELSE
           kmzij = km-k1+1; 
           hhsij = (km-1) * ddz; 
         ENDIF
         IF (idt .GE. 0.01 .AND. ddz .GE. 0.01) THEN ! idt real
           WRITE (UNIT=i60, FMT=1) i,j,kmzij,idt,hhsij,ddz,         &
                                    & iexplt,itrap,cd,ismooth,beta,niter,iextrp, &
                                    & f,tramp,iupwind
         ELSE ! idt real
           WRITE (UNIT=i60, FMT=8) i,j,kmzij,idt,hhsij,ddz,         &
                                    & iexplt,itrap,cd,ismooth,beta,niter,iextrp, &
                                    & f,tramp,iupwind
         ENDIF ! idt real
       1 FORMAT( "i =",I4,"  j =",I4, "  km =", I4,"   dt =",F5.2," sec",     & ! idt real
                  & "  hhs =", F6.3," m","   dz =", F5.2," m"/                & ! idt real
                  & "iexplt =",I2, "  itrap =", I2, "  cd = ", F7.4, 2X,      &
                  & "ismooth =", I2, "  beta =", F6.3," niter =", I2/         &
                  & "iextrp =",  I2, "  f =", F7.4, "  tramp=", F9.1, 2X,     &
                  & "iupwind =", I2 )
       8 FORMAT( "i =",I4,"  j =",I4, "  km =", I4,"   dt =",F7.4," sec",     & ! idt real
                  & "  hhs =", F6.3," m","   dz =", F7.4," m"/                & ! idt real
                  & "iexplt =",I2, "  itrap =", I2, "  cd = ", F7.4, 2X,      &
                  & "ismooth =", I2, "  beta =", F6.3," niter =", I2/         &
                  & "iextrp =",  I2, "  f =", F7.4, "  tramp=", F9.1, 2X,     &
                  & "iupwind =", I2 )
         WRITE (UNIT=i60, FMT=2)     
       2 FORMAT( 1X,"   time     ","  step     ","  zeta ","       SOD      ","       kt       ", &
                    "   depth   ","    u       ","  v      "," w       ",  &
                    "   Av       ","        Dv      ","  scalar    "," Tracers-> " )
         WRITE (UNIT=i60, FMT=3)
       3 FORMAT( 1X,"    hrs     ","   no      ","   cm       ","  mg/m2/s      ","      m/s      ", &
                    " m    ","   cm/s   "  ,"   cm/s   " ,"  cm/s      " , &
                    " cm2/s      ","     cm2/s     ","  oC        ","  g/l   -> " )
         
      END DO
   END IF
 
   !.....Output values at time step  n  to timefile(s).....
   DO nn = 1, nnodes 

      i = inode(nn); 
      j = jnode(nn); 
      i60 = i6 + nn; 

      IF ( i > 0 .AND. j > 0) THEN
        kmzij = kmz(i,j); 
      ELSE
        kmzij = km; 
      ENDIF

      ! ... Map (i,j)- into l-index
      l = ij2l(i,j)

      IF (i > 0 .AND. j > 0) THEN
        ! ... Initialize output variables to -99.  
        uout  = -99.0E-2 ! 10-2 since the output is cm /s
        vout  = -99.0E-2
        wout  = -99.0E-2
        Avout = -99.0E-4 ! 10-4 since the output is cm2/s
        Dvout = -99.0E-4
        uhout = -99.0
        scout = -99.0
        trout = -99.0   
        DO k  = k1, kmzij
          IF (h(k,l)<=ZERO) CYCLE  
          !uout(k)  = 0.5 * (u  (k,l) + u  (k,lWC(l))) 
          !vout(k)  = 0.5 * (v  (k,l) + v  (k,lSC(l)))
          !wout(k)  = 0.5 * (wp (k,l) + wp (k+1,l   ))
          uout(k)  = u (k, l)
          vout(k)  = v (k, l)
          wout(k)  = wp(k, l)
          Avout(k) = 0.5 * (Av (k,l) + Av (k+1,l   ))
          Dvout(k) = 0.5 * (Dv (k,l) + Dv (k+1,l   ))
          uhout(k) = 0.5 * (uh (k,l) + uh (k,lWC(l)))
          scout(k) = sal(k,l)
          IF (ntr>0) THEN
            DO it = 1, nFABMtr_out
              trout(k,it) = tracer(k,l,id_tracer_out(it))
            ENDDO
          ENDIF
       END DO
     ELSE
        ! ... Initialize output variables to -99.
        uout  = 0.E0 
        vout  = 0.E0
        wout  = 0.E0        
        Dvout = 0.E0
        uhout = 0.E0
        scout = 0.E0
        trout = 0.E0   
        DO k  = k1, kmzij
          nk = 0; 
          DO l = 1,lm
            IF (h(k,l)<=ZERO) CYCLE 
            nk = nk + 1 
            uout(k)  = uout(k)+u (k, l)
            vout(k)  = vout(k)+v (k, l)
            wout(k)  = wout(k)+wp(k, l)
            Avout(k) = Avout(k)+0.5 * (Av (k,l) + Av (k+1,l   ))
            Dvout(k) = Dvout(k)+0.5 * (Dv (k,l) + Dv (k+1,l   ))
            uhout(k) = uhout(k)+0.5 * (uh (k,l) + uh (k,lWC(l)))
            scout(k) = scout(k)+sal(k,l)
            IF (ntr>0) THEN
              DO it = 1, nFABMtr_out
                trout(k,it) = trout(k,it)+tracer(k,l,id_tracer_out(it))
              ENDDO
            ENDIF
          ENDDO
          IF (nk > 0) THEN
            uout(k)  = uout(k)/nk
            vout(k)  = vout(k)/nk
            wout(k)  = wout(k)/nk
            Avout(k) = Avout(k)/nk
            Dvout(k) = Dvout(k)/nk
            uhout(k) = uhout(k)/nk
            scout(k) = scout(k)/nk
            IF (ntr>0) THEN
              DO it = 1, nFABMtr_out
                trout(k,it) = trout(k,it)/nk
              ENDDO
            ENDIF
          ELSE
            uout(k)  = -99.E-2
            vout(k)  = -99.E-2
            wout(k)  = -99.E-2
            Avout(k) = -99.E-4
            Dvout(k) = -99.E-4
            uhout(k) = -99.
            scout(k) = -99.
            IF (ntr>0) THEN
              DO it = 1, nFABMtr_out
                trout(k,it) = -99.
              ENDDO
            ENDIF
          ENDIF
        ENDDO
     ENDIF

     IF (thrs == 0) THEN
       k4sodout(i,j)=0.0
     ENDIF

     ! ... Write variables to output file
     IF (ntr <= 0) THEN
       WRITE (UNIT=i60, FMT=4) thrs, n, s(i,j), k4sodout(i,j),     &
            & kt(i,j), (zlevel(k+1), uout (k), vout(k),            &
            & wout(k), Avout(k), Dvout(k), scout(k),               &
            & k =k1,kmzij)
     ELSE
       WRITE (UNIT=i60, FMT=5) thrs, n, s(i,j), k4sodout(i,j),     &
            & kt(i,j), (zlevel(k+1), uout (k), vout(k), wout(k),   &
            & Avout(k), Dvout(k),                                  & !
            & scout(k  ),                                          & 
            & trout(k,1),                                          &
            & trout(k,2),                                          &
            & trout(k,3),                                          &
            & trout(k,4),                                          &
            & trout(k,5),                                          &
            & trout(k,6),                                          &
            & trout(k,7),                                          &
            & trout(k,8),                                          &
            & trout(k,9),                                          &
            & trout(k,10),                                         &
            & trout(k,11),                                         &
            & trout(k,12),                                         &
            & trout(k,13),                                         &
            & trout(k,14),                                         &
            & trout(k,15),                                         &
            & k = k1,kmzij)
     ENDIF
   4 FORMAT(1X,F10.3,I10,2PF9.2,1PE15.5,1PE15.5,0PF9.2,2(2PF10.2),2PF9.4,        &
                        &  2(4PF15.7),  0PE15.7 / ( 60X,0PF9.2,2(2PF10.2),2PF9.4,&
                           2(4PF15.7),  0PE15.7 ))
   5 FORMAT(1X,F10.3,I10,2PF9.2,1PE15.5,1PE15.5,0PF9.2,2(2PF10.2),2PF9.4,        &
                        & 2(4PF15.7),16(0PF15.7)/ ( 60X,0PF9.2,2(2PF10.2),2PF9.4,&
                        & 2(4PF15.7),16(0PF15.7)))
   END DO
    
   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_outt = t_outt + (etime - btime)

END SUBROUTINE outt


!***********************************************************************
SUBROUTINE outw
!***********************************************************************
!
!  Purpose: To write the wind field provided to the model as boundary 
!           condition 
!
!-----------------------------------------------------------------------

  !.....Local variables.....
  CHARACTER (LEN = 10) :: wind_file 
  INTEGER, PARAMETER   :: wind_id0 = 420 
  INTEGER :: wind_id
  INTEGER :: i, j, k, ios, istat, n_frames, k_out, m1,m2, kp, kb
  INTEGER, SAVE :: ipoints
  REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array 
  INTEGER :: year_out, day_out, mon_out 
  REAL    :: hour_out

  ! ... Return if meteorological field is uniform - 
  IF ( ifsurfbc < 10 ) RETURN

  !..... Do only on first entry....
  IF( n == 0 ) THEN

    ! Open spacefile on first entry and print specifications for reading the file 
    n_frames = nts/MAX(iop,1)
    wind_id = wind_id0 
    wind_file = 'wfield.bnr'
    OPEN(unit=wind_id,file=wind_file,FORM='UNFORMATTED',IOSTAT=ios)
    IF(ios /= 0) THEN; PRINT *, "Error opening wind file ", ios; STOP; ENDIF
    !... Write number of time slices to output file
    WRITE(wind_id) n_frames
    ! ... Find out the number of surface points & store it
    ipoints = 0
    DO j = j1, jm; DO i = i1, im
      ! Ignore dry cells
      IF (.NOT. mask2d(i,j)  ) CYCLE
      ipoints = ipoints + 1
    END DO; END DO
    WRITE(wind_id) ipoints  

    ! ... Time stamp
    year_out = iyr
    mon_out  = imon
    day_out  = iday
    hour_out = ihr
   
    ! ... Output sheets
    ALLOCATE( out_array ( ipoints, 4 ), STAT=istat )
    IF (istat /= 0) THEN; PRINT *, 'ERROR allocating space for out_array'; STOP; ENDIF
    k_out = 0;
    DO j = j1, jm; DO i = i1, im
      ! Ignore dry cells
      IF (.NOT. mask2d(i,j)) CYCLE
      ! Update counter
      k_out = k_out + 1;        
      ! Save wind stress components into output variable 
      out_array(k_out,1) = FLOAT(i)
      out_array(k_out,2) = FLOAT(j) 
      out_array(k_out,3) = 0.0 ! Wind stress component in the EW direction
      out_array(k_out,4) = 0.0 ! Wind stress component in the NS direction
    END DO; END DO

    ! ... Id # for plane file
    wind_id = wind_id0 
    ! ... Print time stamp followed by the records
    WRITE(wind_id) n,year_out,mon_out,day_out,hour_out,          &
    &            ((out_array(m1,m2),m2=1,4),m1=1,k_out)
    DEALLOCATE (out_array)
  
  ! ... On successive time steps  
  ELSE

    ! ... Time stamp
    year_out = iyr
    mon_out  = imon
    day_out  = iday
    hour_out = ihr 

    ! ... Output sheets
    ALLOCATE( out_array ( ipoints, 2 ), STAT=istat )
    IF (istat /= 0) THEN; PRINT *, 'ERROR allocating space for out_array'; STOP; ENDIF
   
    k_out = 0; 
    DO j = j1, jm; DO i = i1, im
      ! Ignore dry cells
      IF (.NOT. mask2d(i,j)) CYCLE
      ! Update counter
      k_out = k_out + 1;        
      ! Save wind stress components into output variable
      out_array(k_out,1) = uair(i,j) ! wind component in the EW direction
      out_array(k_out,2) = vair(i,j) ! wind component in the NS direction
    END DO; END DO

    ! ... Id # for plane file
    wind_id = wind_id0 
    ! ... Print time stamp followed by the records
    WRITE(wind_id) n,year_out,mon_out,day_out,hour_out,          &
    &            ((out_array(m1,m2),m2=1,2),m1=1,k_out)
    DEALLOCATE (out_array)
   
END IF

END SUBROUTINE outw

!***********************************************************************
SUBROUTINE outv
!***********************************************************************
!
!  Purpose: To write output @ a cross section to a binary file
!
!-----------------------------------------------------------------------


   !.....Local variables.....
   CHARACTER (LEN = 12) :: section_file 
   INTEGER, PARAMETER   :: section_id0 = 900 
   INTEGER :: section_id, year_out, mon_out, day_out, i, j, k, l
   INTEGER :: m1, m2, ios, istat, n_frames, ipoints, k_out
   REAL    :: hour_out
   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array 

   !.....Open on first entry and print initial conditions ....
   IF( n == 0 ) THEN
   n_frames = nts/MAX(iox,1)
   DO j = 1, n_sections
      section_id = section_id0 + j
      section_file = "section_    .sc"
      IF ( j < 10 ) WRITE ( section_file(9:12), FMT='(I1,"   ")' ) j
      IF(( j >= 10 ) .AND. ( j < 100)) &
      &             WRITE ( section_file(9:12), FMT='(I2,"  " )' ) j
      IF ( j > 100) WRITE ( section_file(9:12), FMT='(I3," "  )' ) j
      OPEN(unit=section_id,file=section_file,FORM='UNFORMATTED',IOSTAT=ios)
      IF(ios /= 0) THEN
        PRINT *, "Error opening section file = ", j, ios; STOP
      ENDIF 
      !... Write number of time slices to output file
      WRITE(section_id) n_frames
      ipoints = 0
      CALL CountSectionCells ( ipoints, j )
      interior_section_points (j) = ipoints
      WRITE(section_id) ipoints, km1
    END DO

   ! ... Time stamp
   year_out = iyr
   mon_out  = imon
   day_out  = iday
   hour_out = ihr

   ! ... Output sections
   DO k = 1, n_sections
     ipoints = interior_section_points (k)
     ALLOCATE( out_array ( ipoints, 11 ), STAT=istat )
     IF (istat /= 0) THEN
       PRINT *, 'ERROR in out_V_plane allocating space for out_array'
     STOP
     ENDIF
     k_out = 0
     DO m1 = 1, n_section_cells (k)
       i = xinode (k,m1)
       j = xjnode (k,m1)
       l = ij2l(i,j)
       IF ( .NOT. mask2D(i,j)) CYCLE 
       DO m2 = k1, kmz(i,j)
         k_out = k_out + 1
         out_array(k_out,1) = FLOAT(i )
         out_array(k_out,2) = FLOAT(j )
         out_array(k_out,3) = FLOAT(m1)
         out_array(k_out,4) = zlevel(m2+1)-hp(m2,l) ! FLOAT(m1)
         out_array(k_out,5) = up  (m2,l) 
         out_array(k_out,6) = vp  (m2,l) 
         out_array(k_out,7) = wp  (m2,l) 
         out_array(k_out,8) = salp(m2,l); 
           IF (hp(m2,l)<= 0.) & 
           out_array(k_out,8) = -99.;
         IF (ntr > 0) THEN
           out_array(k_out,9) = tracerpp(m2,l,LDO);
             IF (hp(m2,l)<= 0.) & 
             out_array(k_out,9) = -99.;
         ELSE            
           out_array(k_out,9) = 0.5 * (Av (m2,l) + Av (m2+1,l))
         ENDIF         
         out_array(k_out,10) = 0.5 * (Dv (m2,l) + Dv (m2+1,l))

         ! ... output Turbulent Kinetic Energy / Fluxes Chris 10/Mar/2017
         out_array(k_out,11) = fluxY (m2,l)
       END DO
     END DO
     ! ... Id # for plane file
     section_id = section_id0 + k
     ! ... Print time stamp followed by the records - mod Chris
     WRITE(section_id) n,year_out,mon_out,day_out,hour_out,        &
    &            ((out_array(m1,m2),m2=1,11),m1=1,ipoints)
     DEALLOCATE (out_array)
   END DO

   ELSE

   ! ... Time stamp
   year_out = iyr
   mon_out  = imon
   day_out  = iday
   hour_out = ihr

   ! ... Output sections
   DO k = 1, n_sections
     ipoints = interior_section_points (k)
     ALLOCATE( out_array ( ipoints, 7 ), STAT=istat )
     IF (istat /= 0) THEN
       PRINT *, 'ERROR in out_V_plane allocating space for out_array'
       STOP
     ENDIF
     k_out = 0
     DO m1 = 1, n_section_cells(k)
       i = xinode (k,m1)
       j = xjnode (k,m1)
       l = ij2l(i,j)
       IF ( .NOT. mask2D(i,j)) CYCLE
       DO m2 = k1, kmz(i,j)
         k_out = k_out + 1
         out_array(k_out,1) = up  (m2,l) 
         out_array(k_out,2) = vp  (m2,l) 
         out_array(k_out,3) = wp  (m2,l) 
         out_array(k_out,4) = salp(m2,l); 
         IF (hp(m2,l)<= 0.) &
           out_array(k_out,4) = -99.;
         IF (ntr > 0) THEN
           out_array(k_out,5) = tracer(m2,l,LDO); ! output DO only
           IF (hp(m2,l)<= 0.) &
             out_array(k_out,5) = -99.;
         ELSE            
           out_array(k_out,5) = 0.5 * (Av (m2,l) + Av (m2+1,l))
         ENDIF
         out_array(k_out,6) = 0.5 * (Dv (m2,l) + Dv (m2+1,l))

         ! ... output Turbulent Kinetic Energy / Fluxes Chris 10/Mar/2017
         out_array(k_out,7) = fluxY (m2,l)
       END DO
     END DO
     ! ... Id # for plane file
     section_id = section_id0 + k
     ! ... Print time stamp followed by the records - mod Chris
     WRITE(section_id) n,year_out,mon_out, day_out,hour_out,          &
     &            ((out_array(m1,m2),m2=1,7),m1=1,ipoints)
     DEALLOCATE (out_array)
   END DO
	
 END IF

END SUBROUTINE outv

!***********************************************************************
SUBROUTINE CountSectionCells ( i1, ksection )
!***********************************************************************
!
!  Purpose: To count wet points in a given X-section (ksection) 
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  02/22/00          F.J. Rueda        Original f90 code    
!  03/27/00          F.J. Rueda        Check for dry cells (not counted)  
!
!-----------------------------------------------------------------------
!

 INTEGER, INTENT (inout) :: i1
 INTEGER, INTENT (in)    :: ksection 

 ! ... Local variables
 INTEGER :: count, i, j, k, m

 count = 0
 DO  m = 1, n_section_cells (ksection)
    i = xinode(ksection,m)
    j = xjnode(ksection,m)
    IF( .NOT. mask2D(i,j) ) THEN
      PRINT *, '****************** WARNING **********************'
      PRINT *, 'Section # ', ksection, ': DRY CELL (',i,j,')'
      PRINT *, '*************************************************'
      CYCLE
    END IF
    count = count + kmz(i,j) - 1
 END DO
 i1 = count
 
END SUBROUTINE CountSectionCells

!***********************************************************************
SUBROUTINE outNB
!***********************************************************************
!
!  Purpose: To write transport and scalar values at a cross section to a 
!           binary file to be used as boundary conditions in nesting 
!           procedure.
!  FJR - 
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER (LEN = 12) :: xfile, Ifile
   INTEGER :: nboid, m1, m2, ios, istat, ipts, ko, iboid
   INTEGER :: i , j , k , l, kmx, kmy, js, je, is, ie, icl 
   INTEGER :: ii, jj, kk, iis, iie, jjs, jje, kks, kke, nn
   REAL    :: dzi, uflow, uflowi, vflow, vflowi, dflow
   REAL, ALLOCATABLE, DIMENSION(:,:) :: outvar 
   INTEGER, PARAMETER:: iboid0 = 20

   ! ... Return if interior boundaries are not requested
   IF (nxNBO <= 0 .OR. ioNBO<= 0 .OR. ioNBTOGGLE <= 0) RETURN

   !.....Open files on first entry & write dims. of the problem
   IF( n == 0 ) THEN
     nfrNBO = nts/MAX(ioNBO,1)
     ntrNBO = ntr
     DO nn = 1, nxNBO
       ! ... Find No. cells in nested boundary and set value of iptNBO
       CALL FindCellsNBO(ipts, nn)
       iptNBO(nn) = ipts
       ! ... Open files
       nboid = nboid0 + nn
       iboid = iboid0 + nn
       xfile = "nbofilex0    "
       Ifile = "nbofilei0    "
       IF ( nn < 10 ) WRITE ( xfile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( xfile( 9:12), FMT='(I2,"  ")' ) nn
       IF ( nn < 10 ) WRITE ( Ifile(10:12), FMT='(I1,"  ")' ) nn
       IF ( nn>= 10 ) WRITE ( Ifile( 9:12), FMT='(I2,"  ")' ) nn
       OPEN(unit=nboid,file=xfile,FORM='UNFORMATTED',IOSTAT=ios)
       OPEN(unit=iboid,file=Ifile,FORM='FORMATTED',IOSTAT=ios)

       IF(ios /= 0) THEN
         PRINT *, "Error opening xfile = ", nn, ios; STOP
       ENDIF 
       !... Write side of nested boundary (i.e. N, S, E or W)
       WRITE(nboid) isdNBO(nn)
       !... Write time information (no. of frames and time between frames) 
       WRITE(nboid) nfrNBO(nn)
       !... Write spatial information (no. of cells) & no. of tracers
       WRITE(nboid) iptNBO(nn), ntrNBO(nn)
     END DO
   ENDIF

   ! ... Assing output variables and write to output files
   DO nn = 1, nxNBO

     ! ... Allocate space for output variables
     ipts = iptNBO(nn)
     ALLOCATE( outvar ( ipts, 5+ntr ), STAT=istat )
     IF (istat /= 0) CALL allocate_error (istat,30)
   
     SELECT CASE (isdNBO(nn))

       ! East or West boundary
       CASE (1,3)

         ! Get i-, j- indexes for bdry. point
         i  = isbcNBO(nn); 
         js = jsbcNBO(nn); 
         je = jebcNBO(nn)

         ! ... Assign variable values to cells within nested grid 
         uflow = 0.0E0
         icl = 0
         DO j = js,je
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells in coarse-grid grid cell
           ii  = ((i-i1)+1)*xxNBO+i1-1 
           jjs = ((j-j1)  )*xxNBO+j1
           jje = ((j-j1)+1)*xxNBO+j1-1
           kmx = kmz(i,j)
           DO jj = jjs, jje 
             DO k = k1, kmx
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF 
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj) 
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = uh (k,l ) 
               outvar(icl,5) = sal(k,l ) 
               IF (ntr > 0) THEN
                 outvar(icl,6:5+nFABMtr_out) = tracer(k,l,id_tracer_out(1:nFABMtr_out))
               ENDIF
               uflow = uflow + uh(k,l)
             ENDDO
           ENDDO         
         ENDDO
         dzi = 0.0E0
         DO i = i1,isbcNBO(nn);
           DO j = j1, jm
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(i,j)
           ENDDO
         ENDDO

       ! North or South boundaries
       CASE (2,4)
    
         ! Get i-, j- indexes for bdry. point
         j  = jsbcNBO(nn); 
         is = isbcNBO(nn); 
         ie = iebcNBO(nn)
    
         ! ... Assign variable values to cells within nested grid 
         icl = 0; uflow = 0.0E0
         DO i = is,ie
           ! ... Get l- from (i,j)
           l   = ij2l(i,j)
           ! ... Indexes for fine-grid cells within coarse-grid grid cell 
           jj  = ((j-j1)+1)*xxNBO+j1-1 
           iis = ((i-i1)  )*xxNBO+i1
           iie = ((i-i1)+1)*xxNBO+i1-1
           kmy = kmz(i,j)
           DO ii = iis, iie
             DO k = k1, kmy
               icl = icl + 1
               IF (icl > ipts) THEN
                 PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
                 STOP
               ENDIF 
               outvar(icl,1) = FLOAT(ii)
               outvar(icl,2) = FLOAT(jj) 
               outvar(icl,3) = FLOAT(k )
               outvar(icl,4) = vh (k,l ) 
               outvar(icl,5) = sal(k,l ) 
               IF (ntr > 0) THEN
                 outvar(icl,6:5+nFABMtr_out) = tracer(k,l,id_tracer_out(1:nFABMtr_out))
               ENDIF
               uflow = uflow + vh(k,l)
             ENDDO
           ENDDO         
         ENDDO

         ! ... Water Surface Elevation --> mass
         dzi = 0.0E0
         DO j = j1,jsbcNBO(nn);
           DO i = i1, im
             IF (.NOT. mask2d(i,j)) CYCLE
             ! ... Get l- from (i,j)
             l   = ij2l(i,j)
             ! ... Add water surface displacements
             dzi = dzi + s(i,j)
           ENDDO
         ENDDO

     END SELECT
 
     ! ... Double check whether icl = ipts
     IF ( icl .NE. ipts) THEN
       PRINT *, ' STOP - Counters in SUBROUTINE outNB do not agree'
       STOP
     ENDIF

     ! ... Id # for nesting boundary file
     nboid = nboid0 + nn
     iboid = iboid0 + nn

     ! ... Write variables to nesting boundary files
     IF((MOD(n,MAX(ioNBO,1)) == 0)) THEN 
       IF (n == 0 ) THEN     
         WRITE(nboid) thrs,((outvar(m1,m2),m2=1,5+ntr),m1=1,ipts)
       ELSE
         WRITE(nboid) thrs,((outvar(m1,m2),m2=4,5+ntr),m1=1,ipts)
       END IF
     ENDIF

     DEALLOCATE (outvar)

     ! ... Mass Balance Check
     WRITE (UNIT=iboid, FMT='(3E20.11)') thrs, uflow, dzi*dx*dy 

   ENDDO

END SUBROUTINE outNB

!***********************************************************************
SUBROUTINE FindCellsNBO ( ix, nn )
!***********************************************************************
!
!  Purpose: To count wett cells ix in a nested grid within a given 
!           X-section nn
!
!-----------------------------------------------------------------------

 ! ... Arguments
 INTEGER, INTENT (inout) :: ix
 INTEGER, INTENT (in)    :: nn 

 ! ... Local variables
 INTEGER :: count, i, j, k, l, kmx, kmy, js, je, is, ie 

 SELECT CASE (isdNBO(nn))

   ! East or West boundary
   CASE (1,3)

     ! Get i-, j- indexes for bdry. point
     i  = isbcNBO(nn); 
     js = jsbcNBO(nn); 
     je = jebcNBO(nn);

     ! Make sure i,j locations are wett cells
     DO j = js, je   
       IF ( .NOT. mask2d(i,j) ) THEN
         PRINT *, "  "
         PRINT *, " ****STOP - EW nesting bdry. DRY"
         STOP
       END IF
     ENDDO

     ! ... Count of No.cells at boundary of nested grid 
     count = 0;
     DO j = js,je
       kmx = kmz(i,j)
       count = count + (kmx-k1+1)
     ENDDO
     count = count * xxNBO

   ! North or South boundary
   CASE (2,4)

     ! ... Get i-, j- indexes for bdry. point
     j  = jsbcNBO(nn); 
     is = isbcNBO(nn); 
     ie = iebcNBO(nn);
     
     ! ... Make sure i,j locations are wett cells
     DO i = is, ie   
       IF ( .NOT. mask2d(i,j) ) THEN
         PRINT *, "  "
         PRINT *, " ****STOP - NS nesting bdry. DRY "
         STOP 
       END IF
     ENDDO

     ! ... Count of No.cells at boundary of nested grid 
     count = 0;
     DO i = is, ie
       kmy = kmz(i,j)
       count = count + (kmy-k1+1) 
     ENDDO
     count = count * xxNBO   

 END SELECT  
 ix = count

 
END SUBROUTINE FindCellsNBO

!***********************************************************************
SUBROUTINE outh
!***********************************************************************
!
!  Purpose: To write output at a specific layer in binary format
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER (LEN = 10) :: plane_file 
   INTEGER, PARAMETER   :: plane_id0 = 800 
   INTEGER :: plane_id
   INTEGER :: i, j, k, l, ios, istat, n_frames, ipoints , k_out, m1,m2, kp, kb
   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array 
   INTEGER :: year_out, day_out, mon_out 
   REAL    :: hour_out

   !.....Open spacefile on first entry and print initial conditions ....
   IF( n == 0 ) THEN

     ! ... Determine No. of time slices to output
     n_frames = nts/MAX(iop,1)

     DO j = 1, n_planes
       ! ... Check that plane no. is below 2
       IF ( p_out(j) < k1 ) THEN
          PRINT *, 'ERROR: Output plane is not in DOMAIN'
          STOP
       END IF
       ! ... Open file to output solution
       plane_id = plane_id0 + j
       plane_file = "plane_    " 
       IF ( p_out(j) < 10 ) WRITE ( plane_file(7:10), FMT='(I1,"   ")' ) p_out(j)
       IF(( p_out(j) >= 10 ) .AND. ( p_out(j) < 100)) &
       &                    WRITE ( plane_file(7:10), FMT='(I2,"  " )' ) p_out(j)
       IF ( p_out(j) > 100) WRITE ( plane_file(7:10), FMT='(I3," "  )' ) p_out(j)
       OPEN(unit=plane_id,file=plane_file,FORM='UNFORMATTED',IOSTAT=ios)
       IF(ios /= 0) PRINT *, "Error opening plane file = ", p_out(j), ios 
       !... Write number of time slices to output file
       WRITE(plane_id) n_frames
       !... Find & write number of wett cells in the plane
       ipoints = 0;
       CALL CountPlaneCells ( ipoints, p_out(j) )
       interior_plane_points (j) = ipoints
       WRITE(plane_id) ipoints
     END DO

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr
  
     ! ... Output planes
     DO k = 1, n_planes

       ipoints = interior_plane_points (k)
       ALLOCATE( out_array ( ipoints, 8 ), STAT=istat )
       IF (istat /= 0) THEN
         PRINT *, 'ERROR allocating space for out_array'
         STOP
       ENDIF
       kp = p_out(k)
       k_out = 0

       ! ... Bottom cell plane 
       IF ( kp > km ) THEN
         DO j = j1, jm; DO i = i1, im
           ! ... Ignore dry cells
           IF (.NOT. mask2d(i,j)) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kmz (i,j)
           l  = ij2l(i,j)
           ! ... Define output variables 
           k_out = k_out + 1
           out_array(k_out,1) = FLOAT(i)
           out_array(k_out,2) = FLOAT(j) 
           out_array(k_out,3) = up (kb,l) 
           out_array(k_out,4) = vp (kb,l)
           out_array(k_out,5) = wp (kb,l)
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           !out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,7) = s(i,j)
           !out_array(k_out,8) = hhs(i,j)
         END DO; END DO
       ! ... Interior plane 
       ELSE 
         DO j = j1, jm; DO i = i1, im
           ! Ignore dry cells
           IF (.NOT. mask2d(i,j)) CYCLE
           IF ( kmz(i,j) < kp   ) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kp
           l  = ij2l(i,j)
           ! ... Define output variables 
           k_out = k_out + 1
           out_array(k_out,1) = FLOAT(i)
           out_array(k_out,2) = FLOAT(j) 
           out_array(k_out,3) = up (kb,l) 
           out_array(k_out,4) = vp (kb,l)
           out_array(k_out,5) = wp (kb,l)
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           !out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,7) = s(i,j)
           !out_array(k_out,8) = hhs(i,j)
         END DO; END DO
       END IF

       ! ... Id # for plane file
       plane_id = plane_id0 + k
       ! ... Print time stamp followed by the records
       WRITE(plane_id) n,year_out,mon_out,day_out,hour_out,          &
       &            ((out_array(m1,m2),m2=1,8),m1=1,ipoints)
       DEALLOCATE (out_array)

     END DO

   ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ! ... Output planes
     DO k = 1, n_planes

       ipoints = interior_plane_points (k)
       ALLOCATE( out_array ( ipoints, 8 ), STAT=istat )
       IF (istat /= 0) THEN
         PRINT *, 'ERROR allocating space for out_array'
         STOP
       ENDIF
       kp = p_out(k)
       k_out = 0

       ! ... Bottom cell plane 
       IF ( kp > km ) THEN
         DO j = j1, jm; DO i = i1, im
           ! ... Ignore dry cells
           IF (.NOT. mask2d(i,j)) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kmz (i,j)
           l  = ij2l(i,j)
           ! ... Define output variables 
           k_out = k_out + 1
           out_array(k_out,3) = up (kb,l) !+ up(kb,lWC(l)))/2.
           out_array(k_out,4) = vp (kb,l) !+ vp(kb,lSC(l)))/2.
           out_array(k_out,5) = wp (kb,l) !+ wp(kb+1,  l ))/2.
           out_array(k_out,6) = salp(kb,l)
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           !out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,7) = s(i,j)
           !out_array(k_out,8) = hhs(i,j)        
         END DO; END DO
       ! ... Interior plane 
       ELSE 
         DO j = j1, jm; DO i = i1, im
           ! Ignore dry cells
           IF (.NOT. mask2d(i,j)) CYCLE
           IF ( kmz(i,j) < kp   ) CYCLE
           ! ... Determine k- & l- indexes for output cell
           kb = kp
           l  = ij2l(i,j)
           ! ... Define output variables 
           k_out = k_out + 1
           out_array(k_out,3) = up (kb,l) !+ up(kb,lWC(l)))/2.
           out_array(k_out,4) = vp (kb,l) !+ vp(kb,lSC(l)))/2.
           out_array(k_out,5) = wp (kb,l) !+ wp(kb+1,  l ))/2.
           out_array(k_out,6) = salp(kb,l); 
           IF(hp(kb,l)<=ZERO) out_array(k_out,6) = -99.
           !out_array(k_out,7) = 0.5*(Av(kb,l)+Av(kb+1,l))
           out_array(k_out,8) = 0.5*(Dv(kb,l)+Dv(kb+1,l))
           out_array(k_out,7) = s(i,j)
           !out_array(k_out,8) = hhs(i,j)
         END DO; END DO
       END IF

       ! ... Id # for plane file
       plane_id = plane_id0 + k
       ! ... Print time stamp followed by the records
       WRITE(plane_id) n,year_out,mon_out,day_out, hour_out,         &
       &            ((out_array(m1,m2),m2=3,8),m1=1,ipoints)
       DEALLOCATE (out_array)
     END DO
   END IF

END SUBROUTINE outh

!***********************************************************************
SUBROUTINE CountPlaneCells ( count, kplane )
!***********************************************************************
!
!  Purpose: To count wet points in a given kplane layer 
!
!-----------------------------------------------------------------------

 ! ... Arguments
 INTEGER, INTENT (inout) :: count
 INTEGER, INTENT (in)    :: kplane 

 ! ... Local variables
 INTEGER :: i, j

 ! ... Code
 count = 0
 IF ( km < kplane ) THEN
   DO j = j1, jm; DO i = i1, im
     ! Ignore dry cells
     IF (.NOT. mask2d(i,j)  ) CYCLE
       count = count + 1
   END DO; END DO
 ELSE
   DO j = j1, jm; DO i = i1, im
     ! Ignore dry cells
     IF (.NOT. mask2d(i,j)  ) CYCLE
     IF ( kmz(i,j) < kplane ) CYCLE; 
       count = count + 1
   END DO; END DO
 END IF
 
END SUBROUTINE CountPlaneCells

!***********************************************************************
SUBROUTINE outz
!***********************************************************************
!
!  Purpose: To write tracer concentrations in computational domain
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER (LEN = 11) :: tracer_file 
   CHARACTER (LEN = 13) :: tracerbc_file 
   INTEGER, PARAMETER   :: tracer_id0 = 1200 
   INTEGER, PARAMETER   :: tracerbc_id0 = 1400
   INTEGER :: tracer_id, tracerbc_id
   INTEGER :: i, j, k, l, k_t,ios, istat, n_frames,  k_out, m1,m2
   INTEGER, SAVE :: ipoints
   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array 
   INTEGER :: year_out, day_out, mon_out 
   REAL    :: hour_out
   INTEGER:: jj
   CHARACTER(LEN=24) :: fmt

   IF ( ntr <= 0 ) RETURN

   !.....Open spacefile on first entry and print initial conditions ....
   IF( n == 0 ) THEN 
    
     ipoints = 0
     DO l = 1, lm
        i = l2i(l); j = l2j(l)
        DO k = k1, kmz(i,j)
           ipoints = ipoints + 1; 
        ENDDO
     ENDDO
 
     DO j = 1, nFABMtr_out
       n_frames = nts/MAX(iotr,1)
       tracer_id = tracer_id0 + j
       tracer_file = "tracer_    "
       IF ( j < 10 ) WRITE ( tracer_file(8:11), FMT='(I1,"   ")' ) j
       IF ((j >=10 ) .AND. (j < 100)) WRITE ( tracer_file(8:11), FMT='(I2,"  ")' ) j
       IF ( j >100 ) WRITE ( tracer_file(8:11), FMT='(I3," ")' ) j
       OPEN(unit=tracer_id,file=tracer_file,FORM='UNFORMATTED',IOSTAT=ios)
       IF(ios /= 0) THEN; PRINT *, "Error opening tracer file = ", j, ios;STOP;ENDIF 
       !... Write number of time slices to output file
       WRITE(tracer_id) n_frames
       WRITE(tracer_id) ipoints
     END DO

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ALLOCATE( out_array ( ipoints, 4 ), STAT=istat )
     IF (istat /= 0) THEN;
       PRINT *, 'ERROR allocating space in output_tracer'
       STOP
     ENDIF
  
     ! ... Output tracer concentrations
     DO k_t = 1, nFABMtr_out
        k_out = 0
        ! ... Assign values to the output array 
        DO j = j1, jm; DO i = i1, im; 
             IF (.NOT. mask2d(i,j)) CYCLE
             l = ij2l(i,j)
             DO k = k1, kmz(i,j)
             k_out = k_out + 1
             out_array(k_out,1) = FLOAT(i)
             out_array(k_out,2) = FLOAT(j) 
             out_array(k_out,3) = FLOAT(k)
             out_array(k_out,4) = tracer(k,l,id_tracer_out(k_t)) !* h(k,l)
             END DO; 
        END DO; END DO
        ! ... Id # for plane file
        tracer_id = tracer_id0 + k_t
        ! ... Print time stamp followed by the records
        WRITE(tracer_id) n,year_out,mon_out, day_out,hour_out,  &
       &            ((out_array(m1,m2),m2=1,4),m1=1,ipoints)
     END DO
     DEALLOCATE (out_array)

   ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ! ... Allocate space
     ALLOCATE( out_array ( ipoints, 1 ), STAT=istat )
     IF (istat /= 0) THEN;
       PRINT *, 'ERROR allocating space in output_tracer'
       STOP
     ENDIF
  
     ! ... Output tracer concentrations
     DO k_t = 1, nFABMtr_out
        k_out = 0
        ! ... Assign values to the output array 
        DO j = j1, jm; DO i = i1, im; 
             IF (.NOT. mask2d(i,j)) CYCLE
             l = ij2l(i,j)
             DO k = k1, kmz(i,j)
             k_out = k_out + 1
             out_array(k_out,1) = tracer(k,l,id_tracer_out(k_t)) !* h(k,l) Chris
             END DO; 
        END DO; END DO
        ! ... Id # for plane file
        tracer_id = tracer_id0 + k_t
        ! ... Print time stamp followed by the records
        WRITE(tracer_id) n,year_out,mon_out, day_out,hour_out,   &
        &            ((out_array(m1,1)),m1=1,ipoints)
     END DO
     DEALLOCATE (out_array)

   END IF

   10 FORMAT( 1X,2I7,F10.3,F40.10)

END SUBROUTINE outz

!***********************************************************************
SUBROUTINE outp
!***********************************************************************
!
!  Purpose: To write the complete solution in the computational domain
!           The resulting file is used to drive particle tracking
!           simulations with PTRACK-TOOL. 
!-----------------------------------------------------------------------

   !.....Local variables.....
   CHARACTER (LEN = 16) :: ptrack_file 
   INTEGER, PARAMETER   :: ptrack_id = 1002 
   INTEGER              :: i, j, k, l, ios, istat, Noframes, m1, m2,  &
                           kout, year_out, day_out, mon_out
   REAL                 :: hour_out
   INTEGER, SAVE        :: ipoints
   REAL, ALLOCATABLE, DIMENSION(:,:) :: out_array 

   IF( n == 0 ) THEN

     ! ... Determine No. of frames to output & No. of interior points
     Noframes = nts/MAX(apxml,1) 
     ipoints = 0
     DO l = 1, lm
        i = l2i(l); j = l2j(l)
        DO k = k1, kmz(i,j)
           ipoints = ipoints + 1; 
        ENDDO
     ENDDO

     ! ... Open output file & print data & initial conditions
     ptrack_file = "ptrack_hydro.bnr"
     OPEN(unit=ptrack_id,file=ptrack_file,FORM='UNFORMATTED',IOSTAT=ios)
     IF(ios /= 0) THEN
       PRINT *, "Error opening ptrack hydro file = ", ios
       STOP
     ENDIF
     !... Write number of time slices to output file
     WRITE(ptrack_id) Noframes
     WRITE(ptrack_id) ipoints

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ALLOCATE( out_array ( ipoints, 8 ), STAT=istat )
     IF (istat /= 0) THEN;
       PRINT *, 'ERROR allocating space in output routine for PTRACK'
       STOP
     ENDIF

     ! ... Output tracer concentrations
     kout = 0
     ! ... Assign values to the output array 
     DO j = j1, jm; DO i = i1, im; 
       IF (.NOT. mask2d(i,j)) CYCLE
       l = ij2l(i,j)
       DO k = k1, kmz(i,j)
         kout = kout + 1
         out_array(kout,1) = FLOAT(i)
         out_array(kout,2) = FLOAT(j) 
         out_array(kout,3) = FLOAT(k)
         out_array(kout,4) = hp (k,l)
         out_array(kout,5) = up (k,l)
         out_array(kout,6) = vp (k,l)
         out_array(kout,7) = wp (k,l)
         out_array(kout,8) = Dv (k,l)
       END DO; 
     END DO; END DO
     ! ... Print time stamp followed by the records
     WRITE(ptrack_id) n,year_out,mon_out, day_out,hour_out,  &
     &              ((out_array(m1,m2),m2=1,8),m1=1,ipoints)
     DEALLOCATE (out_array)

 ELSE

     ! ... Time stamp
     year_out = iyr
     mon_out  = imon
     day_out  = iday
     hour_out = ihr

     ALLOCATE( out_array ( ipoints, 5 ), STAT=istat )
     IF (istat /= 0) THEN;
       PRINT *, 'ERROR allocating space in output routine for PTRACK'
       STOP
     ENDIF

     ! ... Output tracer concentrations
     kout = 0
     ! ... Assign values to the output array 
     DO j = j1, jm; DO i = i1, im; 
       IF (.NOT. mask2d(i,j)) CYCLE
       l = ij2l(i,j)
       DO k = k1, kmz(i,j)
         kout = kout + 1
         out_array(kout,1) = hp(k,l)
         out_array(kout,2) = up(k,l)
         out_array(kout,3) = vp(k,l)
         out_array(kout,4) = wp(k,l)
         out_array(kout,5) = Dv(k,l)
       END DO; 
     END DO; END DO
     ! ... Print time stamp followed by the records
     WRITE(ptrack_id) n,year_out,mon_out, day_out,hour_out,  &
     &              ((out_array(m1,m2),m2=1,5),m1=1,ipoints)
     DEALLOCATE (out_array)

 END IF

END SUBROUTINE outp

!***********************************************************************
SUBROUTINE outs
!***********************************************************************
!
!  Purpose: To write ascii output in xml format to a file used for
!           velocity and particle-tracking animations with the Gr
!           application. The file, called 'spacefile.xml', is 
!           essentially a header file for the sequential binary files
!           (spacefile3d.bin  and  spacefile2d.bin) written out in
!           SUB outs_bin. The binary files contain the 2d and 3d data
!           from all the wet spatial nodes in the solution at 
!           snapshots in time. 
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, ios, numdim_2d, numdim_3d, idt2, ihr2, imin2, isec2
   INTEGER :: int_zeta, int_u, int_v, int_uh, int_vh, int_w, int_dz
   INTEGER, SAVE :: itspf1
   REAL(KIND=DPV) :: thrs_dpv, tdays_dpv
   CHARACTER, SAVE :: date*8, time*10, zone*5
   CHARACTER(LEN=10) :: ch_idt, ch_n, ch_int_zeta, ch_int_u, ch_int_v
   CHARACTER(LEN=10) :: ch_int_w, ch_int_dz, ch_numdim_2d, ch_numdim_3d
   CHARACTER(LEN=20) :: ch_thrs_dpv, ch_tdays_dpv
   CHARACTER(LEN=13) :: output_file_xml  = "spacefile.xml"
   CHARACTER(LEN=15) :: output_file2d_bin= "spacefile2d.bin"
   CHARACTER(LEN=15) :: output_file3d_bin= "spacefile3d.bin"
   INTEGER, DIMENSION(8), SAVE :: values
   LOGICAL, SAVE :: first_entry = .TRUE.
     
   !               -----First entry into subroutine-----
   
   IF ( first_entry ) THEN

      ! ... Only allowed if layers are of uniform thickness
      IF ( ibathyf < 0 ) THEN
        PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
        PRINT *, '*******  Please revise your input files   *************'
        STOP
      ENDIF

      !.....Open the xml spacefile.....
      first_entry = .FALSE.
      OPEN ( UNIT=i96, FILE=output_file_xml, IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file_xml, ios )
      
      !.....Open the binary spacefiles.....
      OPEN ( UNIT=i97, FILE=output_file2d_bin, FORM='UNFORMATTED',        &
	      &   ACCESS='SEQUENTIAL', IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file2d_bin, &
         & ios )
      OPEN ( UNIT=i98, FILE=output_file3d_bin, FORM='UNFORMATTED',        &
	      &   ACCESS='SEQUENTIAL', IOSTAT=ios )
      IF(ios /= 0) CALL open_error ( "Error opening "//output_file3d_bin, &
         & ios )
      
      !.....Get date and time of run for insertion into spacefile.xml.....
      CALL date_and_time ( date, time, zone, values )
                  
      !.....Output titles and run parameters into header records of the xml spacefile.....
      WRITE (UNIT=i96, FMT='("<?xml version=""1.0""?>")' )
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("<!-- Si3D Output for Gr -->")')
      WRITE (UNIT=i96, FMT='("<!-- xml header file for binary spacefiles -->")')
      WRITE (UNIT=i96, FMT=1) idt,ddz,im,jm,itrsca,itrap,cd,ismooth,beta,    &
         &                    niter,itrmom,f,wa
      ! idt,ddz,im,jm,iexplt,itrap,cd,ismooth,beta,    &
      !   &                    niter,iextrp,f,tramp  --->
      ! idt,ddz,im,jm,itrsca,itrap,cd,ismooth,beta,    &
      !   &                    niter,itrmom,f,wa
    !1 FORMAT ( "<!-- ", "Run Parameters:  dt =",I5," s", "  dz =",F5.2, " m",& ! idt real
    1 FORMAT("<!-- ", "Run Parameters:  dt =",F5.2," s", "  dz =",F5.2, " m",& ! idt real
         &     "  im =", I5, "  jm =", I5,/"     itrsca =", I2, "  itrap =", &
         &     I2, "  cd = ", F7.4, "  ismooth =", I2, "  beta =", F6.3,     &
         &     "  niter =", I2/"     itrmom =", I2, "  f =", F7.4,           &
         &     "  wa =", F8.1, " -->" )
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("<gov.usgs.gr>")' )
      WRITE (UNIT=i96, FMT='("   <obj class=""gov.usgs.sfhydro.si3d.Si3dOutput""")')
      WRITE (UNIT=i96, FMT='("        runTitle=""",      A, """" )' ) TRIM(title)
      WRITE (UNIT=i96, FMT='("        runNumber=""", A8,A4, """" )' ) date,time(1:4)
      ! Compute the date and time at which data is first written to the spacefiles
      idt2 = ipxml*idt
      IF (itspf == 0 ) THEN
         itspf1 = itspf
         !CALL compute_date ( 0 ) ! idt real
         CALL compute_date (0.0) ! idt real
      ELSE
         itspf1 = INT(itspf/idt2) * idt2
         IF (MOD(itspf,idt2) /= 0) itspf1 = itspf1 + idt2
         !CALL compute_date ( itspf1 ) ! idt real
         CALL compute_date (FLOAT(itspf1)) ! idt real
      END IF
      ihr2  = INT(isec/3600); isec2 = isec - (ihr2*3600)
      imin2 = INT(isec2/60)
      isec2 = isec2 - (imin2*60)
      ch_idt = int_to_char ((idt2),0)
      WRITE (UNIT=i96, FMT='("        startTime=""",I4,"/",I2.2,"/",I2.2," ",I2.2,    &
         &  ":",I2.2,":",I2.2,""""," timeStepInSeconds=",A,">")' ) iyr, imon, iday,   &
         &  ihr2, imin2, isec2, TRIM(ch_idt)
      CALL compute_date ( 0.0 )     ! reset date to start time of run - idt real
      !CALL compute_date ( 0 )     ! reset date to start time of run ! idt real
         
      !.....Output header records for 2-D variables into the xml spacefile.....
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("      <obj class=""gov.usgs.gr.GrObjectContainer""", &
         &                          " title=""2-D x-y Variables"">")' )
      WRITE (UNIT=i96, FMT='("         <dataseries class=""gov.usgs.gr.data.RafData&
         &Series""")' )
      WRITE (UNIT=i96, FMT='("            file=""", A, """")' ) output_file2d_bin
      numdim_2d = 1   ! For the time being, hardwire for one 2D variable only (Zeta)
      ch_numdim_2d = int_to_char (numdim_2d,0)
      WRITE (UNIT=i96, FMT='("            title=""2-D Variables""",                &
         &                                " numDimensions=", A, ">")' )            &
         &                                TRIM(ch_numdim_2d)
      WRITE (UNIT=i96, FMT='("            <dim num=""0"" title=""Zeta"" units=""met&
         &ers"" unitScale=""0.0001""/>")' )      
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("         </dataseries>")' )
      WRITE (UNIT=i96, FMT='("      </obj>")' )    
           
      !.....Output header records for 3-D variables into the xml spacefile.....
      WRITE (UNIT=i96, FMT='(1X)' )
      WRITE (UNIT=i96, FMT='("      <obj class=""gov.usgs.gr.GrObjectContainer""", &
         &                          " title=""3-D Variables"">")' )
      WRITE (UNIT=i96, FMT='("         <dataseries class=""gov.usgs.gr.data.RafData&
         &Series""")' )
      WRITE (UNIT=i96, FMT='("            file=""", A, """")' ) output_file3d_bin
      numdim_3d = 4   ! For the time being, hardwire for four 3D variables only
      ch_numdim_3d = int_to_char (numdim_3d,0)
      WRITE (UNIT=i96, FMT='("            title=""3-D Variables""",                &
         &                                " numDimensions=", A, ">")' )            &
         &                                TRIM(ch_numdim_3d)
      WRITE(UNIT=i96, FMT='("            <dim num=""0"" title=""u"" ",  &
         & "units=""m/s"" unitScale=""0.0001""/>"    )' )
      WRITE(UNIT=i96, FMT='("            <dim num=""1"" title=""v"" ",  &
         & "units=""m/s"" unitScale=""0.0001""/>"    )' )   
      WRITE (UNIT=i96, FMT='("            <dim num=""2"" title=""w"" ",  &
         & "units=""m/s"" unitScale=""0.00001""/>"  )' )
      WRITE (UNIT=i96, FMT='("            <dim num=""3"" title=""dz"" ", &
         & "units=""m**2/s"" unitScale=""0.0001""/>" )' )
      WRITE (UNIT=i96, FMT='(1X)' )
      
      !.....Output closing xml element end tags to spacefile.xml.....
      WRITE (UNIT=i96, FMT='("         </dataseries>")' )
      WRITE (UNIT=i96, FMT='("      </obj>")' )
      WRITE (UNIT=i96, FMT='("   </obj>")' )
      WRITE (UNIT=i96, FMT='("</gov.usgs.gr>")' )
      WRITE (UNIT=i96, FMT='(1X)' )
      
   END IF
         
         
   !              -----Entry into subroutine at time step n----- 
   
   !.....Do not write outout to spacefiles for
   !     first 'itspf' seconds of the simulation.....
   IF (its >= itspf) THEN

      !.....Write out comments to spacefile.xml
      !     only on first and last time steps.....
      IF ((its == itspf1) .OR. (n >= nts)) THEN

         !.....Compute time variables in double precision for output.....
         thrs_dpv  = REAL( its, dpv )/3600.0_dpv     ! time in hours
         tdays_dpv = thrs_dpv/24.0_dpv + begind      ! time in Julian days
   
         !.....Convert certain variables to CHARACTER values for printing.....
         ch_thrs_dpv = double_to_char(thrs_dpv,4,1)
         ch_tdays_dpv = double_to_char(tdays_dpv,6,1)
         ch_n = int_to_char(n,1)
      
         !.....Output xml comment statements with time stamps.....
         WRITE (UNIT=i96, FMT=2) TRIM(ch_thrs_dpv),            &
            &                    TRIM(ch_n),                   &
            &                    TRIM(ch_tdays_dpv)
       2 FORMAT ( "<!-- ", "time = ", A, " hrs  n = ", A, "  Julian day = ", A, &
            &     " -->" )
         WRITE (UNIT=i96, FMT=3) iyr, imon,iday,ihr
       3 FORMAT ( "<!-- ", "Date/time: ", 3(1X,I4), 1X,I4.4, " -->" )
         WRITE (UNIT=i96, FMT='( "<!-- *** 2-D variables *** -->" )' )
         WRITE (UNIT=i96, FMT='( "<!-- Zeta in cm * 100      -->" )' )
         WRITE (UNIT=i96, FMT='( "<!-- *** 3-D variables *** -->" )' )
         WRITE (UNIT=i96, FMT='( "<!-- u,v   in cm/s * 100   -->" )' )
         WRITE (UNIT=i96, FMT='( "<!-- w     in cm/s * 1000  -->" )' )
         WRITE (UNIT=i96, FMT='( "<!-- dz    in cm**2/s      -->" )' )
         WRITE (UNIT=i96, FMT='(" ")')
      
      END IF
   
      !.....Choose between ascii and binary output for raw data.....
      CALL outs_bin  ! Output data for one time step to binary spacefiles

   END IF

END SUBROUTINE outs

!***********************************************************************
SUBROUTINE outs_bin
!***********************************************************************
!
!  Purpose: To write output in binary format to two spacefiles for velocity
!           and particle-tracking animations with the Gr application. The
!           two files are: spacefile3d.bin and spacefile2d.bin.  3-D
!           variables are stored in the '3d' file and 2-d variables in
!           the '2d' file.  The ordering of nodes written to the file
!           must agree with the  si3d_bathy.xml  file written out in 
!           SUB outg.  The velocities written to the file are taken from
!           the faces of the grid cell control volumes.  The two unformatted
!           files written to in this subroutine are sequential. They are
!           written as 'streaming' data (no record breaks in the file) using 
!           the 'little-endian' bit order for the storage layout. The data
!           values themselves are stored as 2-byte integers. Each 
!           call to the subroutine outputs one time step of data to
!           the spacefiles.  The CALL and file OPEN statements for this 
!           subroutine are within SUB outs_xml.
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, k1s, kms, l
   INTEGER(KIND=INT2) :: int2_zeta, int2_u, int2_v, int2_w, int2_dz

   ! ... Only allowed if layers are of uniform thickness
   IF ( ibathyf < 0 ) THEN
     PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
     PRINT *, '*******  Please revise your input files   *************'
     STOP
   ENDIF
   
   !.....Output zeta values at time step  n  to spacefile2d.bin.....
   ! Loop over all nodes
   DO j = j1, jm
      DO i = i1, im
      ! Ignore dry cells
      IF (.NOT. mask2d(i,j)) CYCLE
      ! Convert zeta to 2-byte integer in meters*10000.
      int2_zeta = NINT(s(i,j)*1.E4, int2)  
      WRITE (UNIT=i97) int2_zeta
      END DO
   END DO
      
   !.....Output velocity values at time step  n  to spacefile3d.bin.....
   ! Loop over all nodes
   DO j = j1, jm
      DO i = i1, im
         ! Ignore dry cells
         IF (.NOT. mask2d(i,j)) CYCLE
         l   = ij2l(i,j)
         k1s = k1z (i,j)
         kms = kmz (i,j)
         DO k = k1s, kms
            ! Convert u,v to 2-byte integers in cm/s*100
            int2_u  = NINT(up (k,l)*1.E4, int2)
            int2_v  = NINT(vp (k,l)*1.E4, int2)
            ! Convert w to a 2-byte integer in cm/s*1000
            int2_w  = NINT(wp (k,l)*1.E5, int2)
            ! Convert dz to a 2-byte integer in cm**2/s
            int2_dz = NINT(Dv(k,l)*1.E4, int2)
            WRITE (UNIT=i98) int2_u, int2_v, int2_w, int2_dz
         END DO
      END DO
   END DO
   
END SUBROUTINE outs_bin

!***********************************************************************
SUBROUTINE outg ( h4 )
!***********************************************************************
!
!  Purpose: To process and output bathymetry to a file for use in graphics 
!           post-processing and particle tracking. The bathymetry output by 
!           this subroutine is processed from the original bathymetry
!           read into the program in SUB bathy.  The processed bathymetry
!           is similar to that processed in SUB bathy for the actual
!           model calculations except that depths are defined at the
!           corners of each grid cell rather than the mid-sides. (For
!           graphics and particle tracking it is necessary to have
!           the bathymetry defined at corners.)  The processed bathymetry
!           is stored in the pointer array h44(1:im1,1:jm1,4).  The corner
!           points of each (i,j) cell are stored in h44 with the following
!           numbering system:
!                               4      3
!                                *----* 
!                                |    |
!                                *----*
!                               1      2
!           This subroutine is called from within SUB bathy. 
!
!  Dummy argument:
!  h4 = Bathymetry pointer array to be processed. A single depth is defined
!       at each cell corner. Prior to being passed into this subroutine,
!       the array h4 has already been modified in SUB bathy (from the
!       original model input bathymetry array) by correcting the datum 
!       and adding a fictitious row/column of depths around the grid. The 
!       dimensions of h4 are (0:im1,0:jm1).
!
!-----------------------------------------------------------------------

   !.....Argument.....
   REAL, DIMENSION(:,:), POINTER :: h4

   !.....Local variables.....
   REAL, DIMENSION(:,:,:), POINTER :: h44
   CHARACTER(LEN=14) :: output_file="              "
   CHARACTER :: date*8, time*10, zone*5, ch_itype*5, ch_ikind*3
   CHARACTER(LEN=10) :: ch_idx, ch_idy, ch_idz, ch_nopen
   !CHARACTER(LEN=10):: ch_nbarr
   CHARACTER(LEN=10) :: ch_isbc, ch_jsbc, ch_iebc, ch_jebc
   !CHARACTER(LEN=10):: ch_isbarr, ch_jsbarr, ch_iebarr, ch_jebarr
   CHARACTER(LEN=14) :: ch_xglobal, ch_yglobal, ch_zglobal, ch_rotation
   INTEGER, DIMENSION(8) :: values
   INTEGER :: ios, istat, i, j, k, kb, niters, niterations, nn, ixml
   INTEGER, SAVE :: i90
   REAL :: ztop, hmin, hmax, hmax_e, hmax_n, rotation, xglobal, yglobal, &
         & zglobal

   !.....Allocate space for the processed bathymetry array h44.....
   ALLOCATE ( h44(1:im1, 1:jm1, 4), stat=istat )
   IF (istat /= 0) CALL allocate_error ( istat, 27 )
   
   ! ... Only allowed if layers are of uniform thickness
   IF ( ibathyf < 0 ) THEN
     PRINT *, 'STOP - Output to xml files only allowed for ibathyf > 0'
     PRINT *, '*******  Please revise your input files   *************'
     STOP
   ENDIF

   !.....Choose output file as either .xml or .txt.....
   output_file = "si3d_bathy.xml"

   !.....Define global (UTM) coordinates of the corner
   !     grid cell and the rotation angle of grid.....
   ! Hardwire these variables (for the time being)
   xglobal  = 0.0  ! Corner grid cell location is taken as the center of
   yglobal  = 0.0  !    fictitious cell (1,1)
   zglobal  = 0.0  ! zglobal is just a vertical datum adjustment
   rotation = 0.0  ! measured in degrees counterclockwise from North
   
   !.....Convert variables for .xml header records to characters.....
   ch_xglobal  = real_to_char(xglobal,1,0)
   ch_yglobal  = real_to_char(yglobal,1,0)
   ch_zglobal  = real_to_char(zglobal,1,0)
   ch_rotation = real_to_char(rotation,1,0)
   !ch_idx   = int_to_char(idx,0) ! idt real
   !ch_idy   = int_to_char(idy,0) ! idt real
   ch_idx   = real_to_char(idx,1,0) ! idt real
   ch_idy   = real_to_char(idy,1,0) ! idt real
   ch_idz   = real_to_char(idz,1,0)
   ch_nopen = int_to_char(nopen,0)
   !ch_nbarr= int_to_char(nbarr,0)

   !.....Adjust all four corners of each cell into the bottom layer.....
   DO j = 1, jm1
      DO i = 1, im1
         ! Ignore dry cells
         IF ( .NOT. mask2d(i,j) ) CYCLE
         kb = kmz(i,j)           ! Store the bottom layer number as kb
         ztop = (kb-k1)*ddz      ! depth to top of bottom layer
         hmin = ztop + dzmin     ! minimum depth allowable
         hmax = ztop + ddz       ! maximum depth allowable
         ! Initialize four corner depths of cell in h44 array
         h44(i,j,1) = h4(i-1,j-1); h44(i,j,2) = h4(i,j-1)
         h44(i,j,3) = h4(i  ,j  ); h44(i,j,4) = h4(i-1,j)
         ! Make sure the corner depths are not
         ! less than hmin or greater then hmax
         DO k = 1, 4
            IF (h44(i,j,k) < hmin) h44(i,j,k) = hmin
            IF (h44(i,j,k) > hmax) h44(i,j,k) = hmax
         END DO
      END DO
   END DO

   !.....Adjust corner depths so that vertical faces in
   !     the bottom profile only begin at layer boundaries.....
   niterations = 2
   DO niters = 1, niterations
      ! Loop over all cells in the grid
      DO j = 1, jm1
         DO i = 1, im1
            ! Ignore dry cells
            IF ( .NOT. mask2d(i,j) ) CYCLE
            kb = kmz(i,j)
            ztop = (kb-k1)*ddz
            hmin = ztop + dzmin
            hmax = ztop + ddz
            ! Work on the east side of cell
            IF ( i /= im1 ) THEN
               ! Check if the adjacent cell to the east has
               ! more layers (is deeper) than the present cell
               IF ( kmz(i+1,j) > kb ) THEN
                  ! Change depths in present cell to hmax
                  h44(i,j,2:3) = hmax
               END IF
               ! Check if the adjacent cell to the east has fewer
               ! layers (is shallower) than the present cell
               IF ( kmz(i+1,j) < kb ) THEN
                  hmax_e = (kmz(i+1,j)-k1+1)*ddz
                  ! Change depths in eastern cell to hmax_e
                  h44(i+1,j,1) = hmax_e;  h44(i+1,j,4) = hmax_e
               END IF
               ! Check if the adjacent cell to the east has the
               ! same number of layers as the present cell
               IF ( kmz(i+1,j) == kb ) THEN
                  ! Equate depths along common cell boundary
                  h44(i+1,j,1) = h44(i,j,2);  h44(i+1,j,4) = h44(i,j,3)
               END IF
            END IF
            ! Work on the north side of cell
            IF ( j /= jm1 ) THEN
               ! Check if the adjacent cell to the north has
               ! more layers (is deeper) than the present cell
               IF ( kmz(i,j+1) > kb ) THEN
                  ! Change depths in present cell to hmax
                  h44(i,j,3:4) = hmax
               END IF
               ! Check if the adjacent cell to the north has fewer
               ! layers (is shallower) than the present cell
               IF ( kmz(i,j+1) < kb ) THEN
                  hmax_n = (kmz(i,j+1)-k1+1)*ddz
                  ! Change depths in northern cell to hmax_n
                  h44(i,j+1,1:2) = hmax_n
               END IF
               ! Check if the adjacent cell to the north has the
               ! same number of layers as the present cell
               IF ( kmz(i,j+1) == kb ) THEN
                  ! Equate depths along common cell boundary
                  h44(i,j+1,1) = h44(i,j,4);  h44(i,j+1,2) = h44(i,j,3)
               END IF
            END IF
         END DO
      END DO
   END DO

   !.....Open output bathymetry file.....
   i90 = i6 + 30
   OPEN (UNIT=i90, FILE=output_file, IOSTAT=ios)
   IF(ios /= 0) CALL open_error ( "Error opening "//output_file, ios )

   !.....Get date and time of run for insertion into header.....
   CALL date_and_time ( date, time, zone, values )

   !.....Output header records to .xml bathymetry file.....
   WRITE (UNIT=i90, FMT='("<?xml version=""1.0""?>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!-- Si3D Processed Bathymetry for Gr -->")')
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<gov.usgs.gr>")' )
   WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.sfhydro.si3d.Bathymetry""")')
   WRITE (UNIT=i90, FMT='("       runTitle=""",      A, """" )' ) TRIM(title)
   WRITE (UNIT=i90, FMT='("       runNumber=""", A8,A4, """" )' ) date,time(1:4)
   WRITE (UNIT=i90, FMT='("       xSpacing=",A," ySpacing=",A," zSpacing=",A, &
      &  " rotation=", A, ">")' ) TRIM(ch_idx),                               & 
      &                           TRIM(ch_idy),                               &
      &                           TRIM(ch_idz),                               &
      &                           TRIM(ch_rotation)
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!-- Define global coordinates of grid corner -->")')
   WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.gr.data.Location"""  &
      &                           " x=",A," y=",A," z=",A,"/>" )' )     &
      &                           TRIM(ch_xglobal),                     &
      &                           TRIM(ch_yglobal),                     &
      &                           TRIM(ch_zglobal)
   WRITE (UNIT=i90, FMT='(1X)' )

   ! Output open boundary information into the file if nopen>0
   IF (nopen > 0) THEN
      WRITE (UNIT=i90, FMT='("<!-- Define node numbers (i,j) for start and end&
         & points of open boundaries -->")' )
      WRITE (UNIT=i90, FMT='("  <obj class=""gov.usgs.gr.GrObjectContainer""")')
      WRITE (UNIT=i90, FMT='("       title=""Open Boundaries"" nopen=",A,">")')&
         &                           TRIM(ch_nopen)
      DO nn = 1, nopen
         IF ( itype(nn) == 1) ch_itype = '"wse"'
         IF ( itype(nn) == 2) ch_itype = '"flw"'
         WRITE (UNIT=i90, FMT='("     <dataseries title=""Boundary", I2, """", &
            & " numDimensions=""2"" itype=", A, ">")' ) nn, ch_itype
         ! Convert node numbers from integers to character variables
         ch_isbc = int_to_char(isbc(nn),1); ch_iebc = int_to_char(iebc(nn),1)
         ch_jsbc = int_to_char(jsbc(nn),1); ch_jebc = int_to_char(jebc(nn),1)
         ! Output start and end node numbers of each open boundary
         WRITE (UNIT=i90, FMT='(A,T6,A)' ) TRIM(ch_isbc), TRIM(ch_jsbc)
         WRITE (UNIT=i90, FMT='(A,T6,A)' ) TRIM(ch_iebc), TRIM(ch_jebc)
         WRITE (UNIT=i90, FMT='("     </dataseries>")' )
      END DO
      WRITE (UNIT=i90, FMT='("  </obj>")' )   ! end open boundary object
      WRITE (UNIT=i90, FMT='(1X)' )
   END IF

   ! Output header information for bathymetry dataseries      
   WRITE (UNIT=i90, FMT='("<!-- Bathymetry data series -->")' )
   WRITE (UNIT=i90, FMT='("  <dataseries title=""Bathymetry""",  &
      &                     " numDimensions=""7"">")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""0"" title=""I""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""1"" title=""J""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""2"" title=""Number of Wet Layers""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""3"" title=""SW depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""4"" title=""SE depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""5"" title=""NE depth""/>")' )
   WRITE (UNIT=i90, FMT='("     <dim num=""6"" title=""NW depth""/>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("<!--i    j nwlayers hsw     hse     hne     hnw -->")')

   !.....Output depths at the corners of each grid cell in meters.....
   DO j = 1, jm1
      DO i = 1, im1
         ! Ignore dry cells
         IF ( .NOT. mask2d(i,j) ) CYCLE
         WRITE (UNIT=i90, FMT='(3I5,4F8.1)' ) i,j,(kmz(i,j)-1),(h44(i,j,k), k=1,4)
      END DO
   END DO
   
   !.....Output closing xml element end tags.....
   WRITE (UNIT=i90, FMT='("   </dataseries>")' )
   WRITE (UNIT=i90, FMT='(1X)' )
   WRITE (UNIT=i90, FMT='("  </obj>")' )   ! End bathymetry object
   WRITE (UNIT=i90, FMT='("</gov.usgs.gr>")' )
     
   !.....Deallocate h44 pointer array.....
   DEALLOCATE ( h44 )

END SUBROUTINE outg

!************************************************************************
FUNCTION int_to_char ( ivalue, noquote )
!************************************************************************
!
!  Purpose: To convert an integer value into a 10-character string. Quotes
!           around the number are added if noquote=0.  The returned string
!           has no leading blanks. Typically the string will have trailing
!           blanks.
!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine 
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (10 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  ivalue = An integer number less than 10**8.
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable: 
!  int_to_char = A 10-character string variable with the character
!                representation of the integer number. The non-blank
!                characters in the string are left-justified.
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: ivalue, noquote
   CHARACTER(LEN=10) :: int_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "
   CHARACTER(LEN=9) :: string = "         "

   !.....Write the integer number to an internal file
   !     and add a trailing quotemark if requested.....
   IF (noquote == 1) THEN
      WRITE ( string, FMT='(I8,A)' ) ivalue, blank
   ELSE
      WRITE ( string, FMT='(I8,A)' ) ivalue, q
   END IF
 
   !.....Add a leading quotemark if requested 
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      int_to_char = ADJUSTL(string)//blank
   ELSE 
      int_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION int_to_char


!************************************************************************
FUNCTION real_to_char ( value, nd, noquote )
!************************************************************************
!
!  Purpose: To convert a real value into a 14-character string. Quotes
!           around the number are added if noquote=0. The value is placed
!           in as readable a format as possible considering its range.  An 
!           exponential format is used for very large or very small numbers.
!           Otherwise, the routine will use "nd" decimal places of accuracy.
!           The returned string has no leading blanks. Typically the string
!           will have trailing blanks.  The routine prints out the number
!           according to the following rules:
!              1. value > 9999999.                 ES12.5
!              2. value < -999999.                 ES12.5
!              3. 0.    < ABS(value) < 0.01        ES12.5
!              4. value = 0.0                      F12.nd
!              5. Otherwise                        F12.nd

!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine 
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (14 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  value = A real number.
!  nd = An integer number specifying the number of decimal places to use in
!       the character representation of the real number. (Must be in the 
!       range: 0<nd<10)
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable: 
!  real_to_char = A 14-character string variable with the character
!                 representation of the real number enclosed in
!                 quotemarks. The non-blank characters in the string
!                 are left-justified.
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: nd, noquote
   REAL,INTENT(IN) :: value
   CHARACTER(LEN=14) :: real_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "   
   CHARACTER(LEN=13) :: string = "             "
   CHARACTER(LEN=13) :: fmt    = "             "

   !.....Select proper format and include a trailing quotemark.....
   IF ( noquote == 0 ) THEN 
      IF ( value > 9999999. ) THEN
         fmt = '(ES12.5,"""")'
      ELSE IF ( value < -999999. ) THEN
         fmt = '(ES12.5,"""")'
      ELSE IF ( value == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """""""""", ") " )' ) nd
      ELSE IF ( ABS(value) < 0.01 ) THEN
         fmt = '(ES12.5,"""")'
      ELSE 
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """""""""", ") " )' ) nd
      END IF
   END IF
   
   !.....Select proper format and do not include a trailing quotemark.....
   IF ( noquote /= 0 ) THEN 
      IF ( value > 9999999. ) THEN
         fmt = '(ES12.5," ")'
      ELSE IF ( value < -999999. ) THEN
         fmt = '(ES12.5," ")'
      ELSE IF ( value == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """ """, ") " )' ) nd
      ELSE IF ( ABS(value) < 0.01 ) THEN
         fmt = '(ES12.5," ")'
      ELSE 
         WRITE ( fmt, FMT='( "(" , "F12." , I1 , "," , """ """, ") " )' ) nd
      END IF
   END IF

   !.....Write the real number to an internal file.....
   WRITE ( string, fmt ) value
 
   !.....Add a leading quotemark if requested
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      real_to_char = ADJUSTL(string)//blank
   ELSE
      real_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION real_to_char


!************************************************************************
FUNCTION double_to_char ( value_dpv, nd, noquote )
!************************************************************************
!
!  Purpose: To convert a double precision value into a 20-character string.
!           Quotes around the number are added if noquote=0. The value is
!           placed in as readable a format as possible considering its range.
!           An exponential format is used for very large or very small numbers.
!           Otherwise, the routine will use "nd" decimal places of accuracy.
!           The returned string has no leading blanks. Typically the string
!           will have trailing blanks.  The routine prints out the number
!           according to the following rules:
!              1. value_dpv > 9999999.                 ES18.11
!              2. value_dpv < -999999.                 ES18.11
!              3. 0.    < ABS(value_dpv) < 0.01        ES18.11
!              4. value_dpv = 0.0                      F18.nd
!              5. Otherwise                            F18.nd

!
!  More details:
!     The quotemarks are placed around the value for use in xml formatted
!     output files. The character string that is returned from this routine 
!     can easily be processed with the TRIM function to remove any
!     trailing blanks.  The length of the returned sting (20 characters)
!     was chosen somewhat arbitrarily.
!
!  Dummy argument:
!  value_dpv = A real number of KIND(1.0D0).
!  nd = An integer number specifying the number of decimal places to use in
!       the character representation of the real number. (Must be in the 
!       range: 0<nd<10)
!  noquote = Parameter indicating whether quotes are placed around the
!            number within the character string or not. ( 1--> no quotes,
!            0--> quotes )
!
!  Result variable: 
!  double_to_char = A 20-character string variable with the character
!                   representation of the double precision number enclosed
!                   in quotemarks. The non-blank characters in the string
!                   are left-justified.
!
!
!------------------------------------------------------------------------

   ! .....Argument.....
   INTEGER,INTENT(IN) :: nd, noquote
   REAL(DPV),INTENT(IN) :: value_dpv
   CHARACTER(LEN=20) :: double_to_char

   !.....Local variables.....
   CHARACTER :: q = '"'
   CHARACTER :: blank = " "   
   CHARACTER(LEN=19) :: string = "                   "
   CHARACTER(LEN=14) :: fmt    = "              "

   !.....Select proper format and include a trailing quotemark.....
   IF ( noquote == 0 ) THEN 
      IF ( value_dpv > 9999999. ) THEN
         fmt = '(ES18.11,"""")'
      ELSE IF ( value_dpv < -999999. ) THEN
         fmt = '(ES18.11,"""")'
      ELSE IF ( value_dpv == 0._dpv) THEN
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """""""""", ")  " )' ) nd
      ELSE IF ( ABS(value_dpv) < 0.01 ) THEN
         fmt = '(ES18.11,"""")'
      ELSE 
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """""""""", ")  " )' ) nd
      END IF
   END IF
   
   !.....Select proper format and do not include a trailing quotemark.....
   IF ( noquote /= 0 ) THEN 
      IF ( value_dpv > 9999999. ) THEN
         fmt = '(ES18.11," ")'
      ELSE IF ( value_dpv < -999999. ) THEN
         fmt = '(ES18.11," ")'
      ELSE IF ( value_dpv == 0.) THEN
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """ """, ") "  )' ) nd
      ELSE IF ( ABS(value_dpv) < 0.01 ) THEN
         fmt = '(ES18.11," ")'
      ELSE 
         WRITE ( fmt, FMT='( "(" , "F18." , I1 , "," , """ """, ") "  )' ) nd
      END IF
   END IF

   !.....Write the double precision number to an internal file.....
   WRITE ( string, fmt ) value_dpv
 
   !.....Add a leading quotemark if requested
   !     and left-justify the result variable.....
   IF (noquote == 1) THEN
      double_to_char = ADJUSTL(string)//blank
   ELSE
      double_to_char = q//ADJUSTL(string)
   END IF

END FUNCTION double_to_char

!***********************************************************************
SUBROUTINE nodech ( i, j, nodeno, nchar )
!***********************************************************************
!
!  Purpose: To convert the node number where model time series output is
!           desired to a character variable. 
!
!  Dummy arguments:
!  i,j    = x and y node mumbers expressed as integers (maximum allowable
!           value for i or j is 9999)
!  nodeno = node number expressed as a character variable
!  nchar  = number of characters in nodeno (excluding trailing blanks)
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!  6/18/01           P.E. Smith        Original f90 code
!
!-----------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN)  :: i, j
   INTEGER, INTENT(OUT) :: nchar
   CHARACTER(LEN=9), INTENT(OUT) :: nodeno

   !.....Convert node numbers to character variables.....
   IF ( i < 10) THEN
      WRITE ( nodeno(1:2), FMT='(I1,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(3:9), FMT='(I1,"      ")' ) j
         nchar = 3
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(3:9), FMT='(I2,"     " )' ) j
         nchar = 4
      END IF
      IF ((j >= 100) .AND. ( j < 1000)) THEN
         WRITE ( nodeno(3:9), FMT='(I3,"    "  )' ) j
         nchar = 5
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(3:9), FMT='(I4,"   "   )' ) j
         nchar = 6
      END IF
   END IF
   IF (( i >= 10 ) .AND. ( i < 100 )) THEN
      WRITE ( nodeno(1:3), FMT='(I2,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(4:9), FMT='(I1,"     ")' ) j
         nchar = 4
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(4:9), FMT='(I2,"    " )' ) j
         nchar = 5
      END IF
      IF ((j >= 100) .AND. ( j < 1000)) THEN
         WRITE ( nodeno(4:9), FMT='(I3,"   "  )' ) j
         nchar = 6
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(4:9), FMT='(I4,"  "   )' ) j
         nchar = 7
      END IF
   END IF
   IF (( i >= 100 ) .AND. ( i < 1000)) THEN
      WRITE ( nodeno(1:4), FMT='(I3,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(5:9), FMT='(I1,"    ")' ) j
         nchar = 5
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(5:9), FMT='(I2,"   " )' ) j
         nchar = 6
      END IF
      IF (( j >= 100) .AND. (j < 1000)) THEN
         WRITE ( nodeno(5:9), FMT='(I3,"  "  )' ) j
         nchar = 7
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(5:9), FMT='(I4," "   )' ) j
         nchar = 8
      END IF
   END IF
   IF ( i >= 1000 ) THEN
      WRITE ( nodeno(1:5), FMT='(I4,"_")' ) i
      IF ( j < 10 ) THEN
         WRITE ( nodeno(6:9), FMT='(I1,"   ")' ) j
         nchar = 6
      END IF
      IF(( j >= 10 ) .AND. ( j < 100 )) THEN
         WRITE ( nodeno(6:9), FMT='(I2,"  " )' ) j
         nchar = 7
      END IF
      IF (( j >= 100) .AND. (j < 1000)) THEN
         WRITE ( nodeno(6:9), FMT='(I3," "  )' ) j
         nchar = 8
      END IF
      IF ( j >= 1000) THEN
         WRITE ( nodeno(6:9), FMT='(I4      )' ) j
         nchar = 9
      END IF
   END IF

END SUBROUTINE nodech

!***********************************************************************
SUBROUTINE cputimes
!***********************************************************************
!
!  Purpose: To output CPU times for various subroutines in the model.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Calculate total time in subroutines.....
   tot_subs=t_exmom+t_matmom+t_matcon+t_solver+t_vel+t_exsal+t_salin &
      &    +t_outt+t_turb+t_settrap+t_save

   !.....Print CPU times.....
   PRINT '(A, I10)', " Number of calls to turb  =", n_turb
   PRINT '(A, I10)', " Number of calls to exmom =", n_exmom
   PRINT '(A)' , " "
   PRINT '(A)', "        --CPU Times--"
   PRINT '(A, F10.3, A)', " t_turb   =", t_turb ,  " seconds"
   PRINT '(A, F10.3, A)', " t_exmom  =", t_exmom,  " seconds"
   PRINT '(A, F10.3, A)', " t_matmom =", t_matmom, " seconds"
  !PRINT '(A, F10.3, A)', "   t_trid =", t_trid,   " seconds"
   PRINT '(A, F10.3, A)', " t_matcon =", t_matcon, " seconds" 
   PRINT '(A, F10.3, A)', " t_solver =", t_solver, " seconds"
   PRINT '(A, F10.3, A)', " t_vel    =", t_vel,    " seconds"
   PRINT '(A, F10.3, A)', " t_exsal  =", t_exsal,  " seconds"
   PRINT '(A, F10.3, A)', " t_salin  =", t_salin,  " seconds"
   PRINT '(A, F10.3, A)', " t_settrap=", t_settrap," seconds"
   PRINT '(A, F10.3, A)', " t_outt   =", t_outt ,  " seconds"
   PRINT '(A, F10.3, A)', " t_save   =", t_save ,  " seconds"
   PRINT '(A, F10.3, A)', " t_setmask=", t_setmask," seconds"
   PRINT '(A, F10.3, A)', " tot_subs =", tot_subs, " seconds" 
   PRINT '(A)' , " "

END SUBROUTINE cputimes

!***********************************************************************
SUBROUTINE outlog
!***********************************************************************
!
!  Purpose: To write the current simulation time and time step number
!           to a log file that can be used to externally monitor the
!           progress of the run.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: ios

   !.....Open the log file.....
   OPEN ( UNIT=i82, FILE="si3d_log.txt", IOSTAT = ios )
   IF ( ios /= 0 ) CALL open_error ( "Error opening si3d_log.txt", ios )

   !.....Write time in hours and time step number to file.....
   WRITE (UNIT=i82, FMT='(" time=", F10.4, " hours", "    n=", I6)' ) thrs, n

   !.....Close the log file during each subroutine call to cause
   !     any output held in a buffer to be written to the file.....
   IF ( n < nts ) THEN
      CLOSE (UNIT=i82, STATUS="KEEP")
   ELSE
      CLOSE (UNIT=i82, STATUS="DELETE")
   END IF

END SUBROUTINE outlog

!***********************************************************************
SUBROUTINE check_stopfile
!***********************************************************************
!
!  Purpose: To check whether the execution of the program is to be
!           stopped, and then, if requested, to stop the program.  The
!           routine reads a variable 'istop' that can be set manually
!           during the execution of the program.  The value of 'istop'
!           is initialized to zero.  If during the execution of the
!           program the value of 'istop' is manually reset to 1, the
!           program will close all files and immediately terminate 
!           execution of the simulation. 
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!-----------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: ios, istop
   INTEGER, SAVE :: nfirst = 0

   !.....Open stop file.....
   OPEN (UNIT=i83,FILE="si3d_stop.txt",ACTION="READWRITE",IOSTAT=ios)
   IF ( ios /= 0 ) CALL open_error ( "Error opening si3d_stop.txt", ios )

   !.....Write istop=0 on the first entry into the subroutine.....
   IF ( nfirst == 0 ) THEN
      nfirst = 1
      WRITE (UNIT=i83, FMT='(" istop = 0")')
      ENDFILE i83
      CLOSE (UNIT=i83)
      RETURN
   END IF

   !.....Read the value of istop.....
   READ (UNIT=i83, FMT='(8X, I2)' ) istop
   
   !.....Stop execution of program if istop=1.....
   IF ( istop == 1 ) THEN
      CLOSE (UNIT=i83, STATUS="DELETE")  ! Close and delete stop file
      PRINT *, " ****STOPPING si3d after istop was manually set to 1"
      WRITE (UNIT=i6, FMT='(" ")')
      WRITE (UNIT=i6, FMT='(" ****STOPPING si3d after istop was manually &
         &set to 1")')
      WRITE (UNIT=i6, FMT='(" n = ", I7)') n
      STOP
   END IF

   !.....Close the stop file..... 
   IF ( n <  nts ) THEN
      CLOSE (UNIT=i83)   ! Close stop file every time step for PC 
   ELSE
      CLOSE (UNIT=i83, STATUS="DELETE")   ! Delete file when n=nts
   END IF

END SUBROUTINE check_stopfile

!***********************************************************************
SUBROUTINE compute_date ( itime_sec )
!***********************************************************************
!
!  Purpose: To compute the date and time given the time in seconds from
!           a previous date and time.  The previous date and time are
!           saved from the last call to the subroutine.
!
!  Variables:
!     iyr0,imon0,iday0,ihr0,isec0 -- Start date and time of the run defined
!                                    in the input data file. (Used only
!                                    in the first call to the subroutine.)
!     iyrp,imonp,idayp,ihrp,isecp -- Previous date and time.
!                                    (Saved during the previous call to
!                                    the subroutine.)
!     iyr ,imon ,iday ,ihr ,isec  -- New date and time.
!                                    (Computed in the subroutine and 
!                                    returned to the calling program.)
!     doy, doyp                   -- New & previous times (fractional julian day)
!
!  More details:
!     iyr  = year  (4 digit integer)
!     imon = month (2 digit integer)
!     iday = day   (2 digit integer)
!     hrs  = time in decimal hours from beginning of iday (real)
!     ihr  = time in decimal hours multiplied by 100 and rounded to the
!            nearest integer (4 digit integer) (ihr=NINT(hrs*100.))
!     isec = time in seconds from beginning of iday (5 digit integer)
!
!  Dummy argument:
!  itime_sec = time in seconds from the previous date and time.
!
!-----------------------------------------------------------------------

   !.....Argument.....
   REAL, INTENT(IN) :: itime_sec ! idt real

   !.....Local variables.....
   INTEGER, DIMENSION(12) :: month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   INTEGER, SAVE :: iyrp, imonp, idayp, ihrp  ! idt real
   REAL   , SAVE :: isecp                     ! idt real
   INTEGER :: max_possible_months, jmonths
   REAL :: hrs

   !.....Initialize date on first entry into subroutine
   !     with the starting date of the simulation.....
   IF (itime_sec == 0) THEN

      ! . Year, month, day, hour, and seconds
      iyrp=iyr0; imonp=imon0; idayp=iday0; ihrp=ihr0; isecp=isec0
      iyr =iyr0; imon =imon0; iday =iday0; ihr =ihr0; isec =isec0

      ! . Julian day
      IF (leap_year(iyrp)) month(2) = 29
      IF (imon0 == 1) THEN
        doyp = iday0+FLOAT(ihr)/2400.
        doy  = iday0+FLOAT(ihr)/2400.
      ELSE
        doyp = SUM(month(1:imon0-1))+iday0+FLOAT(ihr)/2400.
        doy  = SUM(month(1:imon0-1))+iday0+FLOAT(ihr)/2400.
      ENDIF
      RETURN
 
   END IF

   ! ... Store previous values for day of year doy
   doyp = doy;

   !.....Compute new date.....
   isecp = isecp + itime_sec
   IF (isecp < 86400) THEN
      isec = isecp; 
      hrs  = isec/3600.; 
      ihrp = NINT(hrs*100.); 
      ihr  = ihrp
      !doy  = doyp + FLOAT(itime_sec)/86400. ! idt real
      doy  = doyp +       (itime_sec)/86400. ! idt real
      RETURN
   END IF
   idayp = idayp + isecp/86400
   !isecp = MOD(isecp,86400); isec = isecp  ! idt real
   isecp = MOD(isecp,86400.); isec = isecp  ! idt real
   hrs   = isec/3600.
   ihrp  = NINT(hrs*100.); 
   ihr = ihrp
   ! Take care of leap year
   IF (imonp == 2) THEN
      month(2) = 28
      IF (leap_year(iyrp)) month(2) = 29
   END IF
   IF (idayp <= month(imonp)) THEN
      iday = idayp
      IF (imon0 == 1) THEN
        !doy = iday+FLOAT(isec)/86400. ! idt real
        doy = iday+      (isec)/86400. ! idt real
      ELSE
        !doy = SUM(month(1:imon0-1))+iday0+FLOAT(isec)/86400. ! idt real
        doy = SUM(month(1:imon0-1))+iday0+      (isec)/86400. ! idt real
      ENDIF
      RETURN
   END IF
   ! Take care of case where 'itime_sec' is very large (> 1 month)
   max_possible_months = (itime_sec/2500000) + 1
   DO jmonths = 1, max_possible_months
      idayp = idayp - month(imonp)
      imonp = imonp + 1
      IF (imonp == 13) THEN
         imonp = 1; iyrp = iyrp + 1
         IF (leap_year(iyrp)) month(2) = 29
         IF (.NOT. leap_year(iyrp)) month(2) = 28
      END IF
      IF (idayp <= month(imonp)) THEN
         iday = idayp
         EXIT
      END IF
   END DO
   IF (imon0 == 1) THEN
     !doy  = iday+FLOAT(isec)/86400. ! idt real
     doy  = iday+      (isec)/86400. ! idt real
   ELSE
     !doy  = SUM(month(1:imon0-1))+iday0+FLOAT(isec)/86400.  ! idt real
     doy  = SUM(month(1:imon0-1))+iday0+      (isec)/86400.  ! idt real
   ENDIF
   imon = imonp
   iyr  = iyrp

END SUBROUTINE compute_date

!***********************************************************************
LOGICAL FUNCTION leap_year ( iyear )
!***********************************************************************
!
!  Purpose: To identify whether 'iyear' is a leap year or not. The
!           function returns .TRUE. if 'iyear' is a leap year.
!
!  (Note: Formally, a leap year must be divisible by 4 but not by 100,
!         or it must be divisible by 400.  For this subroutine, I only
!         check if a year is divisible by 4.) 
! 
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------

   !.....Argument.....
   INTEGER, INTENT(IN) :: iyear

   leap_year  = MOD(iyear,4) == 0

END FUNCTION leap_year

!************************************************************************
  SUBROUTINE OutScalarEnergyBalance
!************************************************************************
!
!  Purpose: To output mass & energy existing at the basin 
!           scale.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!---------------------------------------------------------------

   ! ... Local Variables
   CHARACTER(LEN=25):: flux_file="ScalarBalance.txt"
   INTEGER, SAVE	:: i899 = 899
   INTEGER              :: ios, i,j,k,l,k1s, kms
   REAL                 :: elev, uijk, vijk, wijk, rijk
   REAL(real_G1)	:: HeatB, OxygB, PotE, KinE
   REAL, PARAMETER	:: SpecificHeat = 4181.6

   ! ... Open output file on first call
   IF ( n == 0 ) THEN
      OPEN ( UNIT=i899, FILE=flux_file, IOSTAT=ios)
      IF(ios /= 0) PRINT *, "Error opening "//flux_file
      ShearCum = 0.0
      BuocyCum = 0.0
      DisspCum = 0.0      
   END IF

   ! ... Heat & Oxygen in the domain
   HeatB = 0.0E0
   OxygB = 0.0E0
   DO l = 1, lm;  
     i = l2i(l); j = l2j(l);
     kms = kmz(i,j)
     k1s = k1z(i,j)
     DO k = k1s, kms
       HeatB = HeatB+sal   (k,l    )*h(k,l)  
       OxygB = OxygB+h(k,l)*dx*dy     
       IF (ntr > 0) & ! Change the condition iopss to ntr > 0
       OxygB = OxygB+tracer(k,l,ntr)*h(k,l)       
      END DO
   END DO; 


   ! ... Basin scale potential energy
   PotE = 0.0
   DO l = 1, lm;  
     i = l2i(l); j = l2j(l);
     kms = kmz(i,j)
     k1s = k1z(i,j)
     elev = hhs(i,j)-h(kms,l)/2.
     rijk = densty_s(salp(kms,l),0.0)+1000.
     PotE = PotE + rijk*g*elev*h(kms,l) 
     IF (k1s == kms) CYCLE
     DO k = kms-1, k1s, -1
       elev = elev + (hp(k+1,l)+hp(k,l))/2.
       rijk = densty_s(salp(k,l),0.0)+1000.
       PotE = PotE + rijk*g*elev*hp(k,l)             
     END DO
   END DO;
   PotE = PotE*dx*dy

   ! ... Basin scale kinetic energy
   KinE = 0.0
   DO l = 1, lm;  
     i = l2i(l); j = l2j(l);
     kms = kmz(i,j)
     k1s = k1z(i,j)
     ! ... Basin scale potential energy
     uijk = (up(kms,l) + up(kms  ,lWC(l)))/2.
     vijk = (vp(kms,l) + vp(kms  ,lSC(l)))/2.
     wijk = (wp(kms,l) + wp(kms+1,    l ))/2.
     rijk = densty_s(salp(kms,l),0.0)+1000.
     KinE = KinE + 0.5*rijk*(uijk**2.+vijk**2.+wijk**2.)*h(kms,l) 
     IF (k1s == kms) CYCLE
     DO k = kms-1, k1s, -1
       uijk = (up(k,l) + up(k  ,lWC(l)))/2.
       vijk = (vp(k,l) + vp(k  ,lSC(l)))/2.
       wijk = (wp(k,l) + wp(k+1,    l ))/2.
       rijk = densty_s(salp(k,l),0.0)+1000.
       KinE = KinE + 0.5*rijk*(uijk**2.+vijk**2.+wijk**2.)*h(k,l) 
     END DO
   END DO;
   KinE = KinE*dx*dy;

   WRITE (UNIT=i899, FMT='(9E20.11)') thrs, HeatB, OxygB, &
                                      PotE, KinE, TKinE,  &
                                      ShearCum, BuocyCum, DisspCum 

   ! ... Set to zero integral variables
   ShearCum = 0.0
   BuocyCum = 0.0
   DisspCum = 0.0

  END SUBROUTINE OutScalarEnergyBalance

!***********************************************************************
SUBROUTINE exTracer  (nt)
!***********************************************************************
!
!  Purpose: To solve transport equation for tracer, using Flux-limiters.
!           nt denotes tracer number to be solved
!
!-----------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: nt

   ! ... Local variables
   INTEGER :: i, j, k, l, k1s, kms, gamma1, istat
   REAL    :: vel, ratio, C_f, delz, twodt1, hd
   REAL, DIMENSION (4          ) :: ss
   REAL, DIMENSION (km1) :: vsettle
   REAL, DIMENSION (lm1) :: fluxdep
   !REAL, DIMENSION(:,:), POINTER :: hdxpp, hdypp !FJR25032011

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   !.....Allocate hdxpp, hdypp, th, and th1 arrays.....
   !ALLOCATE (hdxpp(km1,lm1), hdypp(km1,lm1), STAT=istat)
   !IF ( istat /= 0) CALL allocate_error ( istat, 26 )

   ! ... Constants used in solution
   twodt1 = twodt*tz

   ! ... Calculate hdxpp & hdypp arrays for diffusion terms & 
   hdxpp = 0.0; 
   hdypp = 0.0;
   DO l = 1, lm

      ! ... 3D-(i,j) indexes for l
      i = l2i(l); j = l2j(l);

      ! ... Retrieve top & bottom wet sal-pts .................
      kms = kmz(i,j)
      k1s = k1z(i,j)

      ! ... Calculate hdxpp & hdypp array at u-&v- pts ........
      DO k = k1, kms
        hdxpp(k,l) = Ax0*hupp(k,l)
        hdypp(k,l) = Ay0*hvpp(k,l)
      ENDDO
                        
   END DO

   vsettle = 0;
   ! ... Calculate vertical settling velocity
   SELECT CASE (ecomod)
   CASE (0,1)
       IF (isinjected == 1) THEN
           vsettle = -wsettle(nt);
       ELSE
           vsettle = 0.0E0
       ENDIF
   CASE (2:3)
      vsettle = -vdep(nt) ! use negative as for settling velocity
   CASE (4)
      vsettle(1:km1) = 0 ! FABM vels are specified in the loop later - Chris rectified 21/11/2016
   END SELECT

   ! ... Initialize ex & flux arrays to zeros
   ex    = 0.0; 
   fluxX = 0.0; 
   fluxY = 0.0; 
   fluxZ = 0.0;

   DO l = 1, lm; 

     ! ... Map l- into (i,j)-indexes .........................
     i = l2i(l); j = l2j(l);
	 
     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(i,j)
     k1s = k1z(i,j)

     IF (ecomod == 4) THEN
         vsettle(k1s:kms) = - wVELtr(i,j,k1s:kms,nt) ! FABM output is positive and should come with a negative sign - Chris rectified 21/11/2016
     ENDIF
 
     DO k = k1s, kms;
	    
       ! ... EW fluxes .......................................
       IF (hup(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2  
         vel  = (uhpp(k,l) + uh(k,l))/2.

         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,    l ,nt); 
         ss(3)  = tracerpp(k,lEC(l),nt); 
         IF (hpp(k,lWC(l))<=ZERO)THEN; ss(1)=ss(2); 
           ELSE; ss(1)=tracerpp(k,lWC(l),nt); ENDIF           
         IF (hpp(k,lEC(lEC(l)))<=ZERO)THEN; ss(4)=ss(3);
           ELSE; ss(4)=tracerpp(k,lEC(lEC(l)),nt); ENDIF

         ! ... Calculate Cf factor to use in flux calculations
         gamma1 = -SIGN (1., vel)
         C_f    = 0.0; 
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes at x-faces
         fluxX(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.

       ENDIF

       ! ... NS fluxes .......................................
       IF (hvp(k,l)> ZERO) THEN

         ! ... Velocities at time n+1/2  
         vel  = (vhpp(k,l) + vh(k,l))/2.
          
         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,        l  ,nt ); 
         ss(3)  = tracerpp(k,    lNC(l) ,nt ); 
         IF (hpp(k,lSC(l))<=ZERO)THEN; ss(1)=ss(2); 
           ELSE; ss(1)=tracerpp(k,lSC(l),nt); ENDIF           
         IF (hpp(k,lNC(lNC(l)))<=ZERO)THEN; ss(4)=ss(3); 
           ELSE; ss(4)=tracerpp(k,lNC(lNC(l)),nt); ENDIF

         ! ... Calculate Cf factor to use in flux calculations
         C_f    = 0.0; ! Default value is for upwinding 
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))
           ! MC flux limiter (VanLeer, 1977)
           !C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF
 
         ! ... Calculate fluxes at y-faces
         fluxY(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/dx*C_f)*(ss(3)-ss(2))/2.

       ENDIF

       ! ... UD fluxes .......................................
       IF (hp(k-1,l)> ZERO) THEN         

         !IF (idbg == 1) PRINT *, " settling vel",vsettle(k), i, j, k, wVELtr(i,j,k,nt)
         ! ... Velocities at time n+1/2 (include settling velocity) Chris: a negative value means going down
         vel  = wp(k,l) + vsettle(k); IF (k == k1s .AND. vel > 0.0) vel = 0.0;

         ! ... Define stencil for scalar transport
         ss(2)  = tracerpp(k,l,nt); 
         ss(3)  = tracerpp(k-1,l,nt)  
         IF(hpp(k-2,l)<=ZERO)THEN;ss(4)=ss(3);
           ELSE;ss(4)=tracerpp(k-2,l,nt);ENDIF      
         IF(hpp(k+1,l)<=ZERO)THEN;ss(1)=ss(2);
           ELSE;ss(1)=tracerpp(k+1,l,nt);ENDIF;  
         ! ... Calculate ratio of slope of solution across interfaces & 
         !     estimate flux limiter
         C_f = 1.0 ! Default method is Lax-Wendroff
         gamma1 = -SIGN (1., vel)
         IF (ss(3) - ss(2) /= 0 ) THEN 
           ratio =(ss(3+gamma1)-ss(2+gamma1))/(ss(3)-ss(2))  
           ! MC flux limiter (VanLeer, 1977)   
           ! C_f   = MAX(0., MIN( 2*ratio, (1+ratio)/2., 2. ))
           ! ... Roe's Superbee Limiter
           C_f = MAX(0., MIN(1.,2.*ratio),MIN(2.,ratio))
         ENDIF

         ! ... Calculate fluxes at z-faces
         delz = (hp(k,l) + hp(k-1,l))/2.
         fluxZ(k,l) = vel/2.*(ss(3)+ss(2))- &
         & ((1.-C_f)*ABS(vel)+vel**2.*twodt1/delz*C_f)*(ss(3)-ss(2))/2.

       ENDIF

     ENDDO; 

     ! ... account for deposition in the sediment
     IF ((isseddep == 1 .AND. vsettle(kms) /= 0) .AND. hp(kms-1,l) > ZERO) THEN

         ! ... Velocities at time n+1/2 (include settling velocity) Chris: a negative value means going down
         vel  = wp(k,l) + vsettle(kms); IF (vel > 0.0) vel = 0.0;
         IF (idbg == 1) PRINT *, " deposition params ",hp(kms,l), i, j
         fluxdep (l) = vel * tracerpp(kms,l,nt);

     ELSE
         fluxdep (l) = 0;
     ENDIF

   ENDDO; 

   ! ... Update ex array with divergence of advective fluxes & diffusion
   DO l = 1, lm; 

     ! ... Map l- into (i,j)-indexes .........................
     i = l2i(l); j = l2j(l);
	 
     ! ... Retrieve top & bottom wet sal-pts .................
     kms = kmz(i,j)
     k1s = k1z(i,j)
 
     DO k = k1s, kms;	   

       !.....Horizontal diffusion.....
       hd= (hdxpp(k,    l )*(tracerpp(k,lEC(l),nt) - tracerpp(k,    l ,nt))       &
         & -hdxpp(k,lWC(l))*(tracerpp(k,    l ,nt) - tracerpp(k,lWC(l),nt)))/dxdx &
         &+(hdypp(k,    l )*(tracerpp(k,lNC(l),nt) - tracerpp(k,    l ,nt))       &
         & -hdypp(k,lSC(l))*(tracerpp(k,    l ,nt) - tracerpp(k,lSC(l),nt)))/dydy

       ex(k,l) =   hpp(k,l)*tracerpp(k,l,nt)/twodt1   & 
                - (fluxX(k,l) - fluxX(k,lWC(l))) / dx &
                - (fluxY(k,l) - fluxY(k,lSC(l))) / dy &
                - (fluxZ(k,l) - fluxZ(k+1,l   )) + hd * ihd

     ENDDO; 

     ! calculate removal rate if deposition exists
     ex(kms,l)=ex(kms,l) + fluxdep(l);

   ENDDO



   ! ... Modify explicit term to account for flow boundary conditions
   CALL MODextracer4openbc (nt)

   !.....Deallocate pointers.....
   !DEALLOCATE ( hdypp, hdxpp )

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_exsal = t_exsal + (etime - btime)

END SUBROUTINE ExTracer

!***********************************************************************
SUBROUTINE ImTracer (nt)
!***********************************************************************
!
!  Purpose: To solve for active scalar concentration.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-----------------------------------------------------------------------


   ! ... Arguments
   INTEGER, INTENT (IN) :: nt

   !.....Local variables.....
   REAL :: twodt1, Osource, Qsource
   INTEGER :: i, j, k, l, k1s, kms, kt, nwlayers, inn, kk, noc, nn
   REAL, DIMENSION (1:km1) :: hn
   REAL, DIMENSION (3,1:km1) :: aa
   REAL, DIMENSION (1:km1) :: dsourcesink

   !.....Timing.....
   REAL, EXTERNAL :: TIMER
   REAL :: btime, etime
   btime = TIMER(0.0)

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over interior sal-pts to solve for
   !     matrix from the active scalar equation.....
   DO l = 1, lm
            
      ! ... 3D-(i,j) indexes for l
      i = l2i(l); j = l2j(l);
            
      !.....Compute top & bottom layer numbers & No. of layers ....
      kms = kmz(i,j)
      k1s = k1z(i,j)
      nwlayers = (kms-k1s)+1

      ! ... Define layer thikness at time n - The corrections for
      !     surface and recently submerged cells are needed to 
      !     keep mass conservation
      hn(k1s+1:kms) = h(k1s+1:kms,l)
      hn(k1s      ) = twodt1*wp(k1s,l)+hpp(k1s,l)
      IF (hpp(k1s,l)<= ZERO) THEN
        hn(k1s+1) = hpp(k1s+1,l)
      ENDIF

      ! allocation different types of sources or sinks DISABLED
      IF (nt == -1) THEN
        ! allocate sourcesink value to dsourcesink
        dsourcesink = sourcesink(:,l,nt);
        ! set sourcesink to zero
        sourcesink(:,l,nt) = 0;
      ELSE
        dsourcesink = 0;
      ENDIF


      !.....Calculate active scalar for case of a single layer.....
      SELECT CASE (nwlayers)
      CASE (1)   

        aa( 2,k1s) = hn(k1s)/twodt1
        ds(   k1s) = ex(k1s,l) + sourcesink(k1s,l,nt)
        ! ... Include the effect of Sediment Oxygen Demand 
        !     on bottom cell (Note: dims. MT^-1L^-2, same as ds)          
        IF ( iopss > 0 .AND. nt == ntr) THEN
            ds(k1s) = ds(k1s)
        ENDIF
        tracer(k1s,l,nt) = ds(k1s  )/aa(2,k1s) + dsourcesink(k1s)
        ! Check for Nans
        IF (isnan(tracer(k1s,l,nt))) THEN
            tracer(k1s,l,nt) = 0
        ENDIF

      !.....Calculate active scalar for case of two or more layers.....
      CASE (2:)

         !.....Form coefficient matrix [aa]
         ! Define upper diagonal terms
         aa(3,k1s:kms-1) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(3,kms)       =  0.0
         ! Define lower diagonal terms
         aa(1,k1s+1:kms) = -Dv(k1s+1:kms,l)/(hn(k1s:kms-1)+hn(k1s+1:kms))*2.
         aa(1,k1s)       =  0.0
         ! Define center diagonal terms
         aa(2,k1s:kms)   =  hn(k1s:kms)/twodt1-aa(1,k1s:kms)-aa(3,k1s:kms)
     
         !.....form r.h.s. matrix [ds]..... 
         DO k = k1s, kms
            ds(k) = ex(k,l) + sourcesink(k,l,nt)
         ENDDO

          ! ... Modify transport eqs. to accont for sources & sinks. 30Mar2016 not limited to LDO Chris & Jet model inclusion
         IF (iopss > 0) THEN

             DO inn = 1, iopss

                 ! ... tracer modif for adjacent cell to the jet cell - Chris
                 IF (ptype(iodev(inn)) == 0 .AND. inn /= iopss) THEN                                 !  add source to jet only when it's jet model
                     nn  = iodev(inn);
                     IF  ( (l == lEC(ij2l(ipss(inn),jpss(inn))) .AND. uEpss(nn) > 0) .OR.  &         !  east exit of the jet cell
                         & (l == lWC(ij2l(ipss(inn),jpss(inn))) .AND. uWpss(nn) > 0) .OR.  &
                         & (l == lNC(ij2l(ipss(inn),jpss(inn))) .AND. vNpss(nn) > 0) .OR.  &
                         & (l == lSC(ij2l(ipss(inn),jpss(inn))) .AND. vSpss(nn) > 0)       ) THEN

                         DO k = k1s, kms
                             ds(k)=ds(k) + Qpss(k,inn)*Rpss(k,inn,nt)/(dx*dy)*(1-SQRT(1-(sineagl(nn))**2))/2     !  diffusion angle correction same coeff applied on wp
                         ENDDO

                     ENDIF
                 ENDIF

                 IF ( j /= jpss(inn) .OR. i /=ipss(inn) ) CYCLE
                 SELECT CASE (ptype(iodev(inn)))
                     ! ... account the mixing of jet model
                     CASE (0 )
                         DO k = k1s, kms
                             IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                             Qsource  = Qpss(k,inn)/(dx*dy)  &! Inflow (velocity) (m/s)
                             & *(uWpss(nn)+uEpss(nn)+vNpss(nn)+vSpss(nn))
                             Osource  = Rpss(k,inn,nt)       ! Concentration (kg/m3)
                             ds(k)=ds(k)+Qsource*Osource*SQRT(1-(sineagl(nn))**2)     ! kg/m2/s = conc.* thickness / time - Corrected for flow angles by the jet model
                         ENDDO
                     CASE (1:)
                         DO k = k1s, kms
                             IF (ABS(Qpss(k,inn))<1.E-10) CYCLE
                             Qsource  = Qpss(k,inn)/(dx*dy)  ! Inflow (velocity) (m/s)
                             Osource  = Rpss(k,inn,nt)       ! Concentration (kg/m3)
                             ds(k)=ds(k)+Qsource*Osource     ! kg/m2/s = conc.* thickness / time
                         ENDDO
                 END SELECT
             ENDDO

         ENDIF

         !.....Solve tridiagonal system for the
         !     vertical distribution of active scalar.....
         CALL trid1 (aa, ds, sal1, k1s, kms, km1, nwlayers)
 
         !.....Define scalars at new time step....
         tracer(k1s:kms  ,l,nt) = sal1(1:nwlayers) + dsourcesink(1:nwlayers);
         tracer(k1 :k1s-1,l,nt) = sal1(1         ) + dsourcesink(1);


         DO k = k1s, kms
             IF (isnan(tracer(k,l,nt)) .OR. tracer(k,l,nt) < 0) THEN
                 tracer(k,l,nt) = 0
             ENDIF
         ENDDO

      END SELECT

   !.....End loop over scalar-pts.....
   END DO

   !.....Compute CPU time spent in subroutine.....
   etime = TIMER(0.0)
   t_salin = t_salin + (etime - btime)

END SUBROUTINE ImTracer

!************************************************************************
SUBROUTINE PointSourceSinkInput
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open any files with data for point sources/sinks - it will
!           read the time series data and assign the initial values at time t=0.  
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, k, l, nn, itr, ios, istat, nptspss, ncdev, psstr
   CHARACTER(LEN=14) :: pssfmt, pssfile

   ! ... Allocate space for device characteristics           
   ALLOCATE (  ptype    (npssdev), &               ! Type of device simulated
               dfL      (npssdev), &               ! Length of diffuser (not allways used)
               JtL      (npssdev), &               ! Length of water jet (in jet model only) Chris
               pdt      (npssdev), &               ! Update frequency of forcing variables

               diammb   (npssdev), &               ! Initial diameter of bubbles
               lambda   (npssdev), &               ! Half-width of the plume
               idetr    (npssdev), &               ! Portion of momentum of the sources/sinks
               eCoeff   (npssdev), &               ! Entrainment coefficient near jet exit
               JetArea  (npssdev), &               ! Area of the jet exit
               isrelcell(npssdev), &               ! is relative cell or not
               plumcell (npssdev), &               ! layer of the diffuser
               STAT = istat)   
   IF (istat /= 0) CALL allocate_error ( istat, 22 )

   ! ... Allocate space for input information on forcing variables pss - 
   ALLOCATE (  flpss (npssdev), scpss (npssdev), &
               uEpss (npssdev), uWpss (npssdev), & 
			   vNpss (npssdev), vSpss (npssdev), &
			   wUmtn (npssdev), wLmtn (npssdev), &
			   sineagl (npssdev),                &
               qthrs (npssdev), STAT = istat)
   IF (istat /= 0) CALL allocate_error ( istat, 22 )
   IF (ntr > 0) THEN
      ALLOCATE (trpss (npssdev,ntr), STAT=istat) 
      IF (istat /= 0) CALL allocate_error ( istat, 23 )
      ! initialize trpss mod by Chris 05/10/2016
      trpss = 0;
   ENDIF

   ! ... Allocate space for variables holding column information
   ALLOCATE (  kdetr (iopss), STAT = istat)
   IF (istat /= 0) CALL allocate_error ( istat, 24 )
    
   !               -----Read files with pss data-----

   ! ----Loop over npssdev to Open & read files with pss data-----  
   DO nn = 1, npssdev

      ! ... Construct file name (beware that nopen cannot be > 99)
      pssfile = "pss0 .txt"
      IF ( nn < 10  ) WRITE ( pssfile(5:5), FMT='(I1)' ) nn
      IF ( nn >= 10 ) WRITE ( pssfile(4:5), FMT='(I2)' ) nn

      ! ... Open IO unit
      OPEN (UNIT=i52, FILE=pssfile, STATUS="old", IOSTAT=ios)
      IF (ios /= 0) CALL open_error ( "Error opening "//pssfile, ios )

      ! Skip over first five header records in open boundary condition file
      READ (UNIT=i52, FMT='(//////)', IOSTAT=ios)
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Read type of source-sink simulated (plumes, boundary conditions, pumped inflows)
      READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) ptype(nn)
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Read number of points in file (it has to be equal in all arrays)
      READ (UNIT=i52, FMT='(10X,I11)', IOSTAT=ios) nptspss
      IF (ios /= 0) CALL input_error ( ios, 48 )

      ! Allocate space for the array of data first time in loop
      IF ( nn == 1) THEN
        ALLOCATE ( varspss(npssdev, ntr+2, nptspss), STAT=istat )
        IF (istat /= 0) CALL allocate_error ( istat, 24 )
        varspss = 0;
      ENDIF


      SELECT CASE (ptype(nn))

      ! ************** Section boundary inflows **************
      CASE (-2)  
        

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )
 
        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 49 )

        ! Write the format of the data records into an internal file
        WRITE (UNIT=pssfmt, FMT='("(10X,",I3,"G11.2)")') ntr+2

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
         READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
              (varspss(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 49 )
        END DO

        ! ... Assign initial values to forcing variables
        flpss (nn) = varspss(nn,1,1) ! Flow rate 
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.)
        IF (ntr > 0) THEN
          DO itr = 1, ntr
            trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load 
          ENDDO
        ENDIF

        ! ... Set idetr to 1 - (default value) 
        idetr(nn) = 1.; 

      ! *************** Cell inflows ***********************
      CASE (-1)   
        
        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )
 
        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 50 )

        ! Write the format of the data records into an internal file
        WRITE (UNIT=pssfmt, FMT='("(10X,",I3,"G11.2)")') ntr+2

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
         READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
              (varspss(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 50 )
        END DO

        ! ... Assign initial values to forcing variables
        flpss (nn) = varspss(nn,1,1) ! Flow rate 
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.)
        IF (ntr > 0) THEN
          DO itr = 1, ntr
            trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load 
          ENDDO
        ENDIF

        ! ... Set idetr to 1 - (default value) 
        idetr(nn) = 1.; 

      ! ************* Water Pumped inflows *****************
      CASE (0) 
!
!        PRINT *, '***************** ERROR *****************'
!        PRINT *, 'Water pumped inflow still NOT incorporated'
!        PRINT *, '***************** ERROR *****************'
!	    STOP
!

      !  ------------------- Chris constructing ---------------
	    ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read how water is from upper cell face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) wUmtn(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read how water is from lower cell face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) wLmtn(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read water exiting angle
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) sineagl(nn)
        IF (ios /= 0) CALL input_error ( ios, 52 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water velocity at detrainment cell is calculated
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) idetr(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read entrainment coefficient
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) eCoeff(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read Area of the jet exit (m2)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) JetArea(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read ambient salinity (constant, uS/cm) -- not used
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) salamb
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read No. of cells from the bottom where jets locate Chris
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) jetcell
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read No. of cells from the jet layer
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) jettopc
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read number of point source columns in file
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) psstr
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read frequency of update
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) pdt(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Determine the validity of the input psstr
        IF (psstr > ntr) THEN
            PRINT *, '***************** ERROR *******************'
            PRINT *, '"psstr" must be less than "ntr" in pss.txt '
            PRINT *, '***************** ERROR *******************'
            STOP
        ENDIF

        ! Write the format of the data records into an internal file
        WRITE (UNIT=pssfmt, FMT='("(10X,",I3,"G11.2)")') psstr+2

        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
           READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
                (varspss(nn,i,j), i=1,psstr+2)
           IF (ios /= 0) CALL input_error ( ios, 51 )
        END DO


        ! assign values to variables
        flpss (nn) = varspss(nn,1,1) ! Flow rate
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.) - Not used here
        DO itr = 1, psstr            ! Tracer loaded
          trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load
        ENDDO

        ! ... Calculate water jet length
        ncdev = 0
        DO j = 1, iopss
          IF (iodev(j) == nn) THEN
            ncdev    = ncdev + 1
          ENDIF
        ENDDO
        JtL(nn) = ncdev * idx ! calculate diffuser length


      ! ******** Plume oxygenation models *******************
      CASE (1:) 
	    
        ! Make sure that, at least, one tracer is simulated 
        ! which corresponds to oxygen 
        IF (ntr < 1) THEN
          PRINT *, '***************** ERROR *****************'
          PRINT *, 'No. of tracers should be > 1 to simulate '
          PRINT *, 'the addition of oxygen through the plumes'
          PRINT *, '***************** ERROR *****************'
          STOP         
        ENDIF

        ! Read how water is detrained at west face
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uEpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water is detrained at north face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) uWpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read how water is detrained at east face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vNpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )
 
        ! Read how water is detrained at south face
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) vSpss(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read flow threshold - to determine when it works and when not
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) qthrs(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )
      
        ! Read how water velocity at detrainment cell is calculated
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) idetr(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read half diffuser length (m)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) lambda(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read initial bubble diameter (mm)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) diammb(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read ambient salinity (constant, uS/cm)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) salamb
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read atmospheric pressure (Pascals)
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) patm
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read the way to locate the diffuser in the layers 0 absolute 1 relative to bottom
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) isrelcell(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Read layers above the bottom where the diffusers locate
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) plumcell(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )
  
        ! Read constant number of the tracer columns
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) psstr
        IF (ios /= 0) CALL input_error ( ios, 51 )

        ! Determine the validity of the input psstr
        IF (psstr > ntr) THEN
            PRINT *, '***************** ERROR *******************'
            PRINT *, '"psstr" must be less than "ntr" in pss.txt '
            PRINT *, '***************** ERROR *******************'
            STOP
        ENDIF

        ! Write the format of the data records into an internal file
        WRITE (UNIT=pssfmt, FMT='("(10X,",I3,"G11.2)")') psstr+2

        ! Read frequency of update 
        READ (UNIT=i52, FMT='(10X,G11.2)', IOSTAT=ios) pdt(nn)
        IF (ios /= 0) CALL input_error ( ios, 51 )
        
        ! Read data array
        READ (UNIT=i52, FMT='(A)', IOSTAT=ios) commentline
        DO j = 1, nptspss
           READ (UNIT=i52, FMT=pssfmt, IOSTAT=ios) &
                (varspss(nn,i,j), i=1,psstr+2)
           IF (ios /= 0) CALL input_error ( ios, 51 )
        END DO

        ! ... Assign values to forcing variables 
        !     (psstr should be >= 1 otherwise the model issues an error message)
        flpss (nn) = varspss(nn,1,1) ! Flow rate 
        scpss (nn) = varspss(nn,2,1) ! Active scalar concentration (temp.) - Not used here
        DO itr = 1, psstr            ! Tracer loads -
          trpss(nn,itr) = varspss(nn,2+itr,1) ! Tracer load 
        ENDDO

        ! ... Calculate diffuser length
        ncdev = 0
        DO j = 1, iopss
          IF (iodev(j) == nn) THEN
            ncdev    = ncdev + 1
          ENDIF
        ENDDO
        dfL(nn) = ncdev * idx 

     END SELECT

     ! ... Close IO unit
     CLOSE (i52)

   ENDDO

END SUBROUTINE PointSourceSinkInput


!************************************************************************
SUBROUTINE PointSourceSinkForcing (nn)
!************************************************************************
!
!  Purpose: To assign values of flow, temperature and tracer loads 
!           in sources and sinks devices
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER, INTENT(IN) :: nn 
   REAL    :: dthrs_pss
   INTEGER :: i, j, k, l, itr
 
   ! ... Time (hrs) between consecutive records in files
   dthrs_pss = dtsecpss/3600.

   ! ... Get new values (NewOpenBC)
   flpss(nn)  = nearestn(0.,thrs,varspss(nn,1,:),dthrs_pss)
   scpss(nn)  = nearestn(0.,thrs,varspss(nn,2,:),dthrs_pss)
   IF (ntr > 0) THEN
     DO itr = 1, ntr
       trpss(nn,itr) = nearestn(0.,thrs,varspss(nn,2+itr,:),dthrs_pss)
      ENDDO
   ENDIF

END SUBROUTINE PointSourceSinkForcing

!***********************************************************************
SUBROUTINE PointSourceSinkSolver
!***********************************************************************
!
!  Purpose: Interface between si3d and other models requiring 
!           sources and sinks (including VT-plume model).
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!   29-may-09        FJRUeda           Modifies call to rwps to be consisten
!                                      with latest version of VTPlume
!   21-jul-10        FJRueda           Version that includes boundary conditions
!                                      & sources sinks of WATER
!
!-----------------------------------------------------------------------

   INTEGER :: nn, inn, i, j, k, kk, l, k1s, kms, nwl, ktop, ksrc, plmdim, itr, njets
   REAL    :: areatot, Tsource, Rsource, jetdiam, injettemp !, ot ! 12-10-13 KAB | Chris - jetdiam
   REAL(8) :: wselev, dfelev, dflgth, hcell, rjulday, lambnot, diamm
   REAL(8) :: elevt, qwt, tplumt, comgpt, tplum0, comgp0, qscfm, frconot
   REAL(8), DIMENSION (1:km1) :: zamb, Tamb, DOamb ! B.C. for plume model
   REAL(8), DIMENSION (1:km1) :: qwd               ! Outflow rate for plume
   REAL(8), DIMENSION (1:km1) :: bwd               ! Radius of circular plume
   REAL   , DIMENSION (1,1,ntr) :: Rsourcetrs
   LOGICAL, SAVE :: DiffON 
   REAL*8 :: oteff

   ! ... initialize the device layer
   ksrc = 0;

   ! ... Return if no points sources/sinks are specified
   IF (iopss <= 0) RETURN 

   ! ... initialise counter
   njets = 0;

   ! ... Define Qpss, Tpss, & Rpss for columns in each device
   DO nn = 1, npssdev

     SELECT CASE (ptype(nn))

     ! ************ Section Boundary Conditions ***********************
     CASE (-2) 

       ! ... Only do computations when there is flow
       IF (ABS(flpss(nn)) < qthrs(nn)) CYCLE

       ! ... Only determine boundary conditions on leapfrog iterations
       IF( (n > 1) .AND. (istep > 1)) CYCLE

       ! ... Interpolate forcing variables (flpss, scpss, trpss) 
       !     to present time using time series input - 
       CALL PointSourceSinkForcing (nn)

       ! ... Find total volume of cells holding boundary conditions
       areatot = 0.0
       DO inn = 1, iopss  
         IF (iodev(inn) .NE. nn) CYCLE  
         ! ... Define i,j,l indexes
         i = ipss(inn); 
         j = jpss(inn); 
         l = ij2l(i,j);
         ! ... Define top and bottom cells & no. of layers 	
         !     The diffuser is set to be one cell above the bottom
         k1s = k1z(i,j) ;
         kms = kmz(i,j) ;
         DO k = k1s, kms
           areatot = areatot + hp(k,l) 
         ENDDO
       ENDDO
     
	   ! ... Determine flow rate, temp. and tracer for each cell
	   !     Inflow rate is pressumed uniform in space
       DO inn = 1, iopss  

         IF (iodev(inn) .NE. nn) CYCLE 

         ! ... Define i,j,l indexes
         i = ipss(inn); 
         j = jpss(inn); 
         l = ij2l(i,j);

         ! ... Define k- indexes
         k1s = k1z(i,j) ;
         kms = kmz(i,j) ;

         ! ... Initialize (also done in subroutine init)
		 Qpss(:,inn  ) = 0.0; 
		 Tpss(:,inn  ) = 0.0; 
		 Rpss(:,inn,:) = 0.0; 
         ! ... Loop over cells in the water column
         DO k = k1s, kms
           Qpss(k,inn) = flpss(nn) * hp(k,l) / areatot 
           IF (scpss(nn)<0.0 .OR. flpss(nn)<=0.0) THEN
             Tpss(k,inn) = salp (k,l)
           ELSE 
             Tpss(k,inn) = scpss(nn )
           ENDIF
           IF (ntr > 0) THEN
               itr = LDO;
               IF (trpss(nn,itr)<0.0 .OR. flpss(nn)<=0.0) THEN 
                   Rpss(k,inn,itr) = tracerpp(k,l,itr)
               ELSE
                   Rpss(k,inn,itr) = trpss(nn,itr)
               ENDIF
           ENDIF
         ENDDO
       ENDDO


     ! ************ Bottom cell Boundary Conditions  ********************
     CASE (-1) 

       ! ... Only do computations when there is flow
       IF (ABS(flpss(nn)) < qthrs(nn)) CYCLE

       ! ... Only determine boundary conditions on leapfrog iterations
       IF( (n > 1) .AND. (istep > 1)) CYCLE

       ! ... Interpolate forcing variables (flpss, scpss, trpss) for each device
       !     to present time,  using time series input - 
       CALL PointSourceSinkForcing (nn)
      
       DO inn = 1, iopss  ! Loop over water columns in device 

         IF (iodev(inn) .NE. nn) CYCLE 

         ! ... Define i,j,l indexes
         i = ipss(inn); 
         j = jpss(inn); 
         l = ij2l(i,j);

         ! ... Define k- indexes
         kms = kmz(i,j) ;

         ! ... Initialize (also done in subroutine init)
         Qpss(:,inn  ) = 0.0; 
         Tpss(:,inn  ) = 0.0; 
         Rpss(:,inn,:) = 0.0; 
         ! ... Set values of vars. at bottom cell
         Qpss(kms,inn  ) = flpss(nn)
         IF (scpss(nn)<0.0 .OR. flpss(nn)<=0.0) THEN
           Tpss(kms,inn) = salp (kms,l)
         ELSE 
           Tpss(kms,inn) = scpss(nn   )
         ENDIF
         IF (ntr > 0) THEN
             itr = LDO;
             IF (trpss(nn,itr)<0.0 .OR. flpss(nn) <= 0.0) THEN 
                 Rpss(kms,inn,itr) = tracerpp(kms,l,itr)
             ELSE
                 Rpss(kms,inn,itr) = trpss(nn,itr)
             ENDIF
         ENDIF
       ENDDO ! Loop over columns in device
   

     ! ************ Water jet systems ****************************
     CASE (0)
         !
         !        PRINT *, '***************** ERROR *****************'
         !        PRINT *, 'Water pumped inflow still NOT incorporated'
         !        PRINT *, '***************** ERROR *****************'
         !	    STOP
         ! ............................ constructing ( Chris )..................
         ! ... ksrc - the layer from the surface where devices locate
         ! ... trpss - oxygen flow rate in water

         ! ... Interpolate forcing variables (flpss, scpss, trpss) for each device
	     !     to present time,  using time series input -


         ! ... Only do computations when there is flow
         IF (ABS(flpss(nn)) < qthrs(nn)) CYCLE

         CALL PointSourceSinkForcing (nn)

         ! ... Calculate water jet releasing induced by diffuser
         IF (ABS(flpss(nn)) < qthrs(nn)) THEN  ! jet OFF
             DiffON = .FALSE.
             DO inn = 1, iopss
                 IF (iodev(inn) .NE. nn) CYCLE
                 Qpss (:,inn) = 0.0
             ENDDO
         ELSE                                  ! Update jet FLOWS
             IF &
                 ( (DiffON == .FALSE.)       .OR.  & ! Jet is TURNED ON
                 (  istep  ==  1            .AND.  & ! Update on first iterations
                 (MOD(n,MAX(pdt(nn),1))==0))) THEN   ! Update every pdt time steps

                 DiffON = .TRUE.                     ! Jet remains ON

                 DO inn = iopss, 1, -1               ! from iopss to 1, iopss should be the place where the pump locates
                     IF (iodev(inn) .NE. nn) CYCLE

                     ! ... set counter for no. of jets
                     njets = njets+1;

                     ! ... Define i,j,l indexes
                     i = ipss(inn);
                     j = jpss(inn);
                     l = ij2l(i,j);

                     ! ... Define k- indexes
                     k1s = k1z(i,j) ;
                     kms = kmz(i,j) ;

                     ! ... Initialise jet/pump location and diameter
                     ksrc    = kms - jetcell;
                     jetdiam = sqrt(JetArea(nn)/3.1416)*2;  ! unit m

                     ! ... Initialize (also done in subroutine init)
                     Qpss(:,    inn  ) = 0.0;
                     Tpss(:,    inn  ) = 0.0;
                     Rpss(:,    inn,:) = 0.0;

                     ! ... Overwrite eCoeff



                     ! ... Overwrite idetr and eCoeff to adjust the velocity at the exit of the cell
                     IF (uWpss (nn) > 0 .AND. uEpss (nn) > 0) THEN  ! EW velocity
                         ! ... alpha entrainment coefficient
                         eCoeff(nn) = ( dy- jetdiam )/6/dx  ! Account for different aspect ratios - Normal distribution
                         ! ... u velocity
                         idetr(nn) = dy*(ddz*jettopc)/JetArea(nn)/(1+6*eCoeff(nn)*dx/2/(jetdiam/2)) ! corrected with the actual vertical length from multicells
                     ELSE IF (vNpss (nn) > 0 .AND. vNpss (nn) > 0) THEN  ! NS velocity
                         ! ... alpha entrainment coefficient
                         eCoeff(nn) = ( dx- jetdiam )/6/dy  ! Account for different aspect ratios - Normal distribution
                         ! ... v velocity
                         idetr(nn) = dx*(ddz*jettopc)/JetArea(nn)/(1+6*eCoeff(nn)*dy/2/(jetdiam/2))
                     ENDIF

                     ! otherwise idetr use default value specified in files


                     IF (inn == iopss) THEN                                                  ! The last cell is the submersible pump for jet model only
                         ! ... for submersible pump
                         Qpss(ksrc-jettopc+1:ksrc, inn  ) = - flpss (nn) * (njets-1)/jettopc !&
                                                            !& * idetr(nn)
                                                                                             ! as a sink use '-'
                         Tpss(:,    inn  ) = salpp(:,l)
                         injettemp         = salpp(ksrc,l)                                   ! define output temperature of the jet
                     ELSE
                         ! ... Assign temperature profile into source/sink
                         Qpss(ksrc-jettopc+1:ksrc, inn  ) = flpss (nn)/jettopc               !* idetr(nn)
                                                                                             ! JtL (nn)?? unit m3/s
                         Tpss(:,    inn  ) = salpp (:   ,     l)                             ! not use temperature in input files
                         Tpss(ksrc-jettopc+1:ksrc, inn  ) = injettemp                        ! overwrite the jet core temperature to ambient temperature
                     ENDIF


                     ! ... Define tracer concentrations in entrained & detrained water

                     IF (ntr > 0) THEN
                         ! initialise all tracers in source/sink as the ambient values.
                         Rpss(ksrc-jettopc+1:ksrc,inn,:) = tracerpp(ksrc-jettopc+1:ksrc,l,:);
                         itr = LDO;
                         IF (trpss(nn,itr)<0.0 .OR. flpss(nn) <= 0.0) THEN
                             Rpss(ksrc-jettopc+1:ksrc,inn,itr) = tracerpp(ksrc-jettopc+1:ksrc,l,itr);
                         ELSE
                             IF (inn == iopss) THEN                   ! assume no oxygen is taken away from the water so no existing oxygen considered later.
                                 Rpss(ksrc-jettopc+1:ksrc,inn,itr) = 0;
                             ELSE                                     ! oxygen concentration(g/m3 or mg/L) take
                                 Rpss(ksrc-jettopc+1:ksrc,inn,itr) = tracerpp(ksrc,l,itr) + trpss(nn,itr)/jettopc;
                             ENDIF
                         ENDIF
                         IF ((MOD(n,MAX(500,1))==0 .OR. n == 1) .OR. idbg == 1) THEN ! update to screen every 500 steps
                             PRINT *, '****************************************************'
                             PRINT *, '-------OUTPUT FROM JET MODEL ROUTINES---------------'
                             PRINT *, '****************************************************'
                             PRINT *, 'LAYER AND FLOW RATE:', inn, njets, eCoeff(nn), injettemp
                             PRINT *, '****************************************************'
                         ENDIF
                     ENDIF
                 ENDDO





             ENDIF
         ENDIF


      
     ! ************ Oxygen-gas diffuser ****************************
     CASE (1:)

       ! ... Interpolate forcing variables (flpss, scpss, trpss) for each device
       !     to present time,  using time series input - 
       CALL PointSourceSinkForcing (nn)
       IF (idbg == 1) PRINT *, "debug point - pss1", kdetr(inn),km1, varspss(nn,1,1:12),varspss(nn,2,1:12)

       ! ... Calculate en-detrainment flows induced by diffuser
       IF (ABS(flpss(nn)) < qthrs(nn)) THEN  ! Diffuser OFF
         DiffON = .FALSE.  
         DO inn = 1, iopss  
           IF (iodev(inn) .NE. nn) CYCLE 
           Qpss (:,inn) = 0.0E0
           kdetr(inn  ) = km1
           IF (idbg == 1) PRINT *, "debug point - pss2", kdetr(inn),km1
         ENDDO
       ELSE                                  ! Update diffuser FLOWS          
         IF & 
         ( (DiffON == .FALSE.)       .OR.  & ! Diffuser is TURNED ON
         (  istep  ==  1            .AND.  & ! Update on first iterations
         (MOD(n,MAX(pdt(nn),1))==0))) THEN   ! Update every pdt time steps

           DiffON = .TRUE.                     ! Diffuser remains ON
           DO inn = 1, iopss  ! Loop over columns in device         
             IF (iodev(inn) .NE. nn) CYCLE 
          
             ! ... Define i,j,l indexes
             i = ipss(inn); 
             j = jpss(inn); 
             l = ij2l(i,j);

             ! ... Define k- indexes
             k1s = k1z(i,j) ;
             kms = kmz(i,j) ;
             nwl = kms-k1s+1;


         
             ! ... Define ambient temperatures
             Tamb(k1s:kms) = salpp(k1s:kms,l)
             Tamb(kms+1  ) = Tamb(kms);
             Tamb(1      ) = Tamb(k1s);
          
             ! ... Define DO concentrations
             DOamb(k1s:kms) = tracerpp(k1s:kms,l,LDO) * DO_mmolpm3_to_mgpl;! unit mg/l
             DOamb(kms+1  ) = DOamb(kms); ! unit mg/l
             DOamb(1      ) = DOamb(k1s); ! unit mg/l
          
             ! ... Depths for cells in plume column from datum
             zamb(k1s  ) = hp(k1s,l)/2.
             DO k = k1s+1, kms
               zamb(k) = zamb(k-1) + (hp(k-1,l)+hp(k,l))/2.
             END DO
             zamb(kms+1) =  zamb(kms)+hp(kms,l)
             zamb(1    ) = -zamb(k1s) 

             IF (isrelcell (nn) == 1) THEN
                 ksrc    = kms-plumcell(nn);   ! Layer No. where diffuser is located from bottom
             ELSE
                 ksrc    = plumcell(nn)        ! Layer No. where diffuser is located from surf
             ENDIF

             IF (idbg == 1) PRINT *, "entering bubble plume module", DOamb(ksrc), flpss(nn)

             ! ... Inputs for plume model
             dfLgth  = dfL(nn)         ;             ! Length of diffuser
             rjulday = doy             ;             ! Julian day (arbitrary)
             wselev  = 0.0000          ;       ! Elevation of free surface
             !ksrc    = plumcell(nn)    ;      ! Layer No. where diffuser is located
             ksrc = kms -1 ;
             dfelev  = -zamb(ksrc)     ;       ! Elevation of diffuser  
             hcell   = ddz             ;       ! Pressumed constant - thickess of cells
             qwd     = 0.0E0           ;       ! Initialize qwd
             bwd     = 0.0E0           ;       ! Initialize perimeter (FJRplume)
             qscfm   = flpss(nn)       ;       ! Air flow rate 
             frconot = 0.965           ;       ! Fraction of O2 in air (not used?); Changed to 0.965 KAB 12-16-13
             lambnot = lambda(nn)      ;       ! Half-width 
             diamm   = diammb(nn)      ;       ! Initial bubble diameter
             IF (ptype(nn) < 2) THEN
               plmdim = 1 ! Linear Plume
             ELSE 
               plmdim = 2 ! Circular Plume
             ENDIF 

             ! ... Run plume model
             CALL lineplu_v1(iyr,rjulday,wselev,dfelev,kms,dfLgth, &
                               lambnot,salamb,patm, diamm, plmdim, &
                               qscfm, frconot,                     & ! boundary conditions
                               ksrc, hcell,                        &
                               zamb (1:kms),                       &
                               Tamb (1:kms),                       &
                               DOamb(1:kms),                       &
                               elevt, qwt , tplumt, comgpt,        &
                               qwd(1:ksrc), bwd(1:ksrc), ktop, oteff) 

             ! ....Detrainment cell - save it for later 
             kdetr(inn) = ktop

             ! ... Define flow at cells
             Qpss(:,inn)   = 0.0; 

             ! ... Define flow at entrainment cells 
             DO k = ktop+1,ksrc  
                Qpss(k,inn) = -qwd(k)*dy/dfLgth 
             ENDDO
    
             ! ... Define flow at detrainment cell to force volume conservation
             Qpss(ktop,inn) = -SUM(Qpss(ktop+1:ksrc,inn)); 

	         ! ... Oxygen transfer efficiency - save for later - KAB 12-10-13
	         oteff = oteff / 100
             !PRINT *, 'O2 Trans Eff (%), nn, inn:', oteff, nn, inn
	         IF (oteff <= 1) THEN
                 oteffnn(nn,inn) = oteff
               ELSE
                 oteffnn(nn,inn) = 1
             ENDIF

             IF ((MOD(n,MAX(100,1))==0 .OR. n == 1) .OR. idbg == 1) THEN ! update to screen every 100 steps
                 PRINT *, '****************************************************'
                 PRINT *, '----OUTPUT FROM PLUME ROUTINES----------------------'
                 PRINT *, '****************************************************'
                 PRINT *, 'Plume number and column number:', nn, inn
                 PRINT *, 'Results of Plume Model elev', elevt, dfelev
                 PRINT *, 'Resutls of Plume Model qwt ', qwt, ktop
                 PRINT *, 'Results of Plume Model T&O ', tplumt, comgpt
                 PRINT *, 'Results of Plume Model pwd ', ksrc, bwd(ksrc), bwd(ktop+1)
                 PRINT *, 'Results of Plume Model O2ef', oteff
                 PRINT *, '****************************************************'
                 PRINT *, 'Plume width from 1-ktop', bwd(ktop+1:ksrc)
             ENDIF
           ENDDO
         ENDIF 
       ENDIF


       IF (idbg == 1) PRINT *, "assigning source terms"
       ! ... Define temperature of entrained and detrained water in plume 
       !     based on existing values of Qpss & temperatures  
       DO inn = 1, iopss  ! Loop over columns in device         
         IF (iodev(inn) .NE. nn) CYCLE 
         
         ! ... Define i,j,l indexes
         i = ipss(inn); 
         j = jpss(inn); 
         l = ij2l(i,j);

         ! ... Define k- indexes
         k1s = k1z(i,j) ;
         kms = kmz(i,j) ;
         nwl = kms-k1s+1;
         Tpss(:,inn) = salpp(:,l)
         k = kdetr(inn); 
         !ot = oteffnn(nn,inn); ! KAB 12-10-13
         IF (k < kms) THEN
             Tsource = 0.0 ! FJRPlumes
             DO kk = k+1,kms
                 Tsource  = Tsource + salpp(kk,l)*Qpss(kk,inn)
             ENDDO
             Tsource = Tsource / SUM(Qpss(k+1:kms,inn))
         ELSE
             Tsource = salpp(k,l)
         ENDIF
         Tpss(k,inn) = Tsource
         !IF (idbg == 1) PRINT *, "debug point 4", kdetr(inn),Tsource,Qpss(kdetr(inn),inn),LDO

         IF (idbg == 1) PRINT *, "assigning DO source...", LDO
         ! ... Define tracer concentrations in entrained & detrained water 
         !     in the plume (assumes dx = dy) based on existing values of 
         !     Qpss, tracer concs. and location of detrainment cell
         IF (ntr > 0) THEN
             DO itr = 1, ntr
                 !itr=LDO; % removed by Chris mod 29Sep2016
                 ! Initialize the tracer sources - Chris mod 15Mar2016 -MOD 28Mar2016
                 ! Rpss(:,    inn,:) = 0.0;
                 ! assign actual concentration to all sources. origin: Rpss(:,inn,itr) = tracerpp(:,l,itr)
                 Rpss(:,inn,:) = tracerpp(:,l,:);
                 k = kdetr(inn);
                 IF (k < kms) THEN
                     Rsource = 0.0
                     Rsourcetrs = 0.0
                     DO kk = k+1,kms ! comment: integrate along the water depth. As for nonsource/sink cell, Qpss=0 bypasses them
                         ! modified by Chris integrate all sources. old version: Rsource  = Rsource + tracerpp(kk,l,itr) * Qpss(kk,inn);
                         Rsourcetrs(1,1,:) = Rsourcetrs(1,1,:) + tracerpp(kk,l,:) * Qpss(kk,inn);
                     ENDDO
                     !IF (idbg == 1) PRINT *, "debug point 5", ksrc, k
                     Rsourcetrs(1,1,:) = Rsourcetrs(1,1,:) / SUM(Qpss(k+1:kms,inn));
                     Rsource = Rsourcetrs(1,1,itr) +  &
                         (trpss(nn,itr)*oteffnn(nn,inn))*dy/dfL(nn)/Qpss(k,inn) ! added oteffnn 12-10-13 KAB
                 ELSE
                     !IF (idbg == 1) PRINT *, "debug point 6", Rsource, kdetr(inn)
                     Rsourcetrs(1,1,:) = tracerpp(k,l,:)
                     Rsource = (trpss(nn,itr)*oteffnn(nn,inn))*dy/dfL(nn)/Qpss(k,inn) ! added oteffnn 12-10-13 KAB
                 ENDIF
                 !IF (idbg == 1) PRINT *, "debug point 7", Rsourcetrs
                 Rpss(k,inn,:) = Rsourcetrs(1,1,:);                 ! added Chris Apr 07 2016
                 Rpss(k,inn,itr) = Rsource;             ! overwrite the oxygen source added by the bubble-plume model
                 !IF (idbg == 1) PRINT *, "debug point 8", Rsource,k,Rpss(k,inn,itr)
             ENDDO
         ENDIF
       ENDDO
          
     END SELECT

   ENDDO

END SUBROUTINE PointSourceSinkSolver

END MODULE si3d_procedures
