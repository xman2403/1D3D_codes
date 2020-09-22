!************************************************************************
MODULE si3d_ecomod
    !************************************************************************
    !
    !  Purpose: Procedures that implement routines related to the modelling
    !           of ecological processees
    !
    !-------------------------------------------------------------------------

    USE si3d_boundaryconditions
    USE si3d_types

    IMPLICIT NONE
    SAVE
  
CONTAINS


    !***********************************************************************
    SUBROUTINE srcsnk00
        !***********************************************************************
        !
        !  Purpose: Set sources and sinks to zero - used when conservative
        !           tracers (ecomod == 0) are modelled
        !
        !-----------------------------------------------------------------------
 
        sourcesink = 0.0E0

    END SUBROUTINE srcsnk00

    !***********************************************************************
    SUBROUTINE TRCinput
        !***********************************************************************
        !
        !  Purpose: To read in values of peak concentrations, location of peaks
        !           and dispersion of clouds -
        !
        !-----------------------------------------------------------------------

        !. . . .Local Variables
        CHARACTER(LEN=12) :: trciofile="si3d_tr.txt"
        CHARACTER(LEN=18) :: iotrcfmt
        INTEGER:: ios, itr ,idep, istat
        INTEGER, PARAMETER :: niotrc = 9

        !.....Open input parameter file.....
        OPEN (UNIT=i99, FILE=trciofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//trciofile, ios )

        ! ... Allocate space for arrays holding information on size groups
        ALLOCATE ( trct0(ntr), trcpk(ntr), trctn(ntr), &
            trcx0(ntr), trcy0(ntr), trcz0(ntr), &
            trcsx(ntr), trcsy(ntr), trcsz(ntr), STAT=istat)
        IF (istat /= 0) CALL allocate_error ( istat, 121 )

        ! Skip over first five header records in SIZE file
        READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! Write the format of the data records into an internal file
        WRITE (UNIT=iotrcfmt, FMT='("(10X,",I3,"G11.2)")') niotrc

        ! Read information for each size group (diameter, growth rate & deposition)
        DO itr = 1, ntr
            READ (UNIT=i99, FMT=iotrcfmt, IOSTAT=ios) &
                trct0(itr), trctn(itr), trcpk(itr),  &
                trcx0(itr), trcy0(itr), trcz0(itr),  &
                trcsx(itr), trcsy(itr), trcsz(itr)
            PRINT *, itr, trct0(itr), trctn(itr)
            IF (ios /= 0) CALL input_error ( ios, 91 )
        END DO

        CLOSE (i99)

    END SUBROUTINE TRCinput

    !***********************************************************************
    SUBROUTINE InitTracerCloud
        !***********************************************************************

        INTEGER:: itr, i, j, k, l
        REAL   :: x, y, z

        DO itr = 1, ntr

            DO l = 1, lm
                i = l2i(l);
                j = l2j(l);
                DO k = k1,kmz(i,j)
                    x = ( REAL ( i - i1 ) + 0.5 ) * dx - trcx0 (itr)
                    y = ( REAL ( j - j1 ) + 0.5 ) * dy - trcy0 (itr)
                    z = ( REAL ( k - k1 ) + 0.5 ) * ddz- trcz0 (itr)
                    tracer  (k,l,itr) = trcpk (itr) *                        &
                        & EXP( - x ** 2. / ( 2 * trcsx(itr) ** 2. ) ) *        &
                        & EXP( - y ** 2. / ( 2 * trcsy(itr) ** 2. ) ) *        &
                        & EXP( - z ** 2. / ( 2 * trcsz(itr) ** 2. ) )
                ENDDO
            ENDDO

        ENDDO

    END SUBROUTINE InitTracerCloud

    !***********************************************************************
    SUBROUTINE SDinput
        !***********************************************************************
        !
        !  Purpose: To read in values of parameters for sediment transport
        !           model E = alpha*(tau-taucr/taucr)^m & vdep
        !
        !-----------------------------------------------------------------------

        !. . . .Local Variables
        CHARACTER(LEN=12) :: sdiofile="si3d_sed.txt"
        CHARACTER(LEN=18) :: iosdfmt
        INTEGER:: ios, itr , idep, istat
        INTEGER, PARAMETER :: niosd = 3

        !.....Open input parameter file.....
        OPEN (UNIT=i99, FILE=sdiofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//sdiofile, ios )

        ! ... Allocate space for arrays holding information for sediments
        ALLOCATE ( sdenSED(ntr), diamSED(ntr), vdep(ntr), STAT=istat)
        IF (istat /= 0) CALL allocate_error ( istat, 99 )

        ! ....Skip over first five header records in SEDIMENT file
        READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! ....Read hours between consecutive WIBSS fields
        READ (UNIT=i99, FMT='(10X,G11.2)', IOSTAT=ios) dthrsWIBSS
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! ....Read hours from start to first WIBSS field
        READ (UNIT=i99, FMT='(10X,G11.2)', IOSTAT=ios) thrsWIBSS0
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! ....Skip over one header records in SEDIMENT file
        READ (UNIT=i99, FMT='(//)', IOSTAT=ios)
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! ....Write the format of the data records into an internal file
        WRITE (UNIT=iosdfmt, FMT='("(10X,",I3,"G11.2)")') niosd

        ! ....Read information for each sediment category
        DO itr = 1, ntr
            READ (UNIT=i99, FMT=iosdfmt, IOSTAT=ios) &
                diamSED(itr), sdenSED(itr), vdep(itr)
            PRINT *, itr, sdenSED(itr), diamSED(itr), vdep(itr)
            IF (ios /= 0) CALL input_error ( ios, 91 )
        END DO

        CLOSE (i99)

        ! ... Set thrsWIBSS to zero
        thrsWIBSSp = 0.0
        thrsWIBSS  = 0.0


    END SUBROUTINE SDinput

    !***********************************************************************
    SUBROUTINE srcsnkSD
        !***********************************************************************
        !
        !  Purpose: To evaluate the source-sink terms in the reactive-transport
        !           equations for sediments - the only sources considered
        !           occur at the bottom cell due to resuspension E
        !
        !-----------------------------------------------------------------------
 
        ! ... Local variables
        INTEGER :: i, j, k, l, k1s, kms, itr, kmyp, kmym, kmxp, kmxm
        REAL    :: Rep, fRep, Zu5, ubott, vbott, ustar, taubx, tauby, taub

        ! ... Initialize to zero on first call
        IF (n .LE. 1) sourcesink = 0.0E0

        ! ... Calculate Wave-Induced Bottom Shear Stress WIBSS fields
        !     (computed with stwave in pre-process mode)
        CALL CalculateWIBSS

        ! ... Loop over tracers
        DO itr = 1, ntr

            ! ... Calculate particle Reynolds number
            Rep  = sqrt(g*sdenSED(itr)*(diamSED(itr)**3.))/1.E-6;
            IF      ( Rep .GT. 1. ) THEN
                fRep = 0.586*(Rep**1.23);
            ELSE IF ( Rep .LE. 1. ) THEN
                fRep = Rep**3.75;
            ENDIF
            IF      ( Rep .GT. 3. ) THEN
                PRINT *, '!!!!!!!!!!!! WARNING !!!!!!!!!!!'
                PRINT *, 'Rep is not within range of      '
                PRINT *, 'validity of sediment model      '
                PRINT *, 'Size Class = ', itr
                PRINT *, 'Rep = ', Rep
            ENDIF

            ! ... Loop over wet columns
            DO l = 1, lm
           
                ! ... 3D-(i,j) indexes for l ...........
                i = l2i(l); j = l2j(l);
            
                !.....Compute bottom layer numbers  ....
                kms = kmz(i,j);
                kmyp= MIN(kmz(i,j+1), kmz(i,j))
                kmym= MIN(kmz(i,j-1), kmz(i,j))
                kmxp= MIN(kmz(i+1,j), kmz(i,j))
                kmxm= MIN(kmz(i-1,j), kmz(i,j))

                ! ....Compute Currend-Induced Bottom Shear Stress CIBSS (cd*rho*U^2)
                ubott = (uhp(kmxp,l)+uhp(kmxm,lWC(l)))/2./hp(kms,l)
                vbott = (vhp(kmyp,l)+vhp(kmym,lSC(l)))/2./hp(kms,l)
                taubx = cd * (rhop(kms,l)+1000.) * ubott**2.
                tauby = cd * (rhop(kms,l)+1000.) * vbott**2.

                ! ... Add Wave-Induced Bottom Shear Stress & calculate friction velocity
                taubx = taubx + wibssx(i,j)
                tauby = tauby + wibssy(i,j)
                taub  = sqrt(taubx**2.+tauby**2.)
                ustar = sqrt(taub/1000.)

                ! ... Compute Erosion flux E (kg/m2/s) at bottom cell
                Zu5  = ( (ustar/vdep(itr)) * fRep )**5.;
                sourcesink(kms,l,itr)= Ased*(Zu5) / ( 1. + (Ased/0.3)*(Zu5))

            ENDDO ! End loop over wet columns

        ENDDO ! End loop over tracers

    END SUBROUTINE srcsnkSD

    !***********************************************************************
    SUBROUTINE CalculateWIBSS
        !***********************************************************************

        !.....Local variables.....
        INTEGER           :: i, j, k, l, ios, istat
        INTEGER, SAVE     :: nn = 0
        REAL              :: weight
        CHARACTER(LEN=12) :: wibssfmt, wibssfile

        ! ... Allocate space & initialize on first call
        IF ( nn .EQ. 0) THEN
            ALLOCATE ( wibssIOx(im1,jm1), wibssIOxp(im1,jm1),         &
                & wibssIOy(im1,jm1), wibssIOyp(im1,jm1),         &
                & wibssx  (im1,jm1), wibssy   (im1,jm1), STAT=istat)
            IF (istat /= 0) CALL allocate_error ( istat, 90 )
            wibssIOx = 0.0E0; wibssIOxp = 0.0E0
            wibssIOy = 0.0E0; wibssIOyp = 0.0E0
            wibssx   = 0.0E0; wibssx    = 0.0E0
        ENDIF

        ! Save and read new frame if thrs > thrsWIBSS.......
        IF (thrs > thrsWIBSS) THEN

            ! ... Update times
            nn = nn + 1
            thrsWIBSSp = thrsWIBSS
            thrsWIBSS  = thrsWIBSS + dthrsWIBSS

            ! ... Set initial time to thrsWIBSS0 for nn = 1
            IF (nn .EQ. 1) thrsWIBSS = thrsWIBSS0

            ! ... Update wibssIOxp & wibssIOyp
            wibssIOxp = wibssIOx
            wibssIOyp = wibssIOy

            ! ... Input file name
            wibssfile = "wibss00 .txt"
            IF (                nn <  10) WRITE ( wibssfile(8:8), FMT='(I1)' ) nn
            IF ( nn >= 10 .AND. nn < 100) WRITE ( wibssfile(7:8), FMT='(I2)' ) nn
            IF ( nn >= 100              ) WRITE ( wibssfile(6:8), FMT='(I3)' ) nn

            !.....Open wibss file.....
            OPEN (UNIT=i99, FILE=wibssfile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
            IF (ios /= 0) CALL open_error ( "Error opening "//wibssfile, ios )

            ! ... Read wibssIOx & wibssIOy
            WRITE (UNIT=wibssfmt, FMT='("(", I4, "G9.0)")') im-i1+1
            DO j = jm, j1, -1
                READ (UNIT=i99, FMT=wibssfmt, IOSTAT=ios) (wibssIOx(i,j), i = i1, im)
                IF (ios /= 0) CALL input_error ( ios, 92 )
            END DO
            DO j = jm, j1, -1
                READ (UNIT=i99, FMT=wibssfmt, IOSTAT=ios) (wibssIOy(i,j), i = i1, im)
                IF (ios /= 0) CALL input_error ( ios, 92 )
            END DO
            CLOSE(i99)

        ENDIF

        ! ... Estimate wibss at present time step by interpolation
        weight  = (thrs - thrsWIBSSp)/(thrsWIBSS-thrsWIBSSp)
        DO l = 1, lm
            i = l2i(l);
            j = l2j(l);
            wibssx (i,j) = wibssIOx (i,j)*(   weight) + &
                wibssIOxp(i,j)*(1.-weight)
            wibssy (i,j) = wibssIOy (i,j)*(   weight) + &
                wibssIOyp(i,j)*(1.-weight)
        ENDDO

    END SUBROUTINE CalculateWIBSS

    !***********************************************************************
    SUBROUTINE SZinput
        !***********************************************************************
        !
        !  Purpose: To read in values for diameter, max growth rate & settling
        !           velocity for each size group (up to ntr)
        !
        !-----------------------------------------------------------------------

        !. . . .Local Variables
        CHARACTER(LEN=12) :: sziofile="si3d_sz.txt"
        CHARACTER(LEN=18) :: ioszfmt
        INTEGER:: ios, itr ,idep, istat
        INTEGER, PARAMETER :: niosize = 4

        !.....Open input parameter file.....
        OPEN (UNIT=i99, FILE=sziofile, STATUS="OLD",  FORM="FORMATTED", IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//sziofile, ios )

        ! ... Allocate space for arrays holding information on size groups
        ALLOCATE ( diamsz(ntr), rmaxsz(ntr), vdep(ntr), STAT=istat)
        IF (istat /= 0) CALL allocate_error ( istat, 99 )

        ! Skip over first five header records in SIZE file
        READ (UNIT=i99, FMT='(/////)', IOSTAT=ios)
        IF (ios /= 0) CALL input_error ( ios, 90 )

        ! Write the format of the data records into an internal file
        WRITE (UNIT=ioszfmt, FMT='("(10X,",I3,"G11.2)")') niosize

        ! Read information for each size group (diameter, growth rate & deposition)
        DO itr = 1, ntr
            READ (UNIT=i99, FMT=ioszfmt, IOSTAT=ios) &
                diamsz(itr), rmaxsz(itr), idep, vdep(itr)
            PRINT *, itr, diamsz(itr), rmaxsz(itr), idep, vdep(itr)
            IF (ios /= 0) CALL input_error ( ios, 91 )
            IF (idep > 0) THEN ! Stokes Law
                vdep(itr) = 10. ** (1.112 * LOG(diamsz(itr))/LOG(10.) - 1.647);                ! Speed in m/day
                vdep(itr) = vdep(itr) / 86400. ! Speed in m/s
            ENDIF
        END DO

        CLOSE (i99)

    END SUBROUTINE SZinput

    !***********************************************************************
    SUBROUTINE srcsnkSZ
        !***********************************************************************
        !
        !  Purpose: To evaluate the source-sink terms in the reactive-transport
        !           equations, used in the analysis of size structures
        !
        !-----------------------------------------------------------------------
 
        ! ... Local variables
        INTEGER :: i, j, k, l, k1s, kms, itr
        REAL    :: r, zt

        ! ... Loop over tracers
        DO itr = 1, ntr

            ! ... Loop over wet columns
            DO l = 1, lm
            
                ! ... 3D-(i,j) indexes for l
                i = l2i(l); j = l2j(l);
            
                !.....Compute top & bottom layer numbers & No. of layers ....
                kms = kmz(i,j);
                k1s = k1z(i,j);

                ! ... Growth rate (express it in s^-1)
                k = k1s;
                zt = hpp(k,l)/2.
                sourcesink(k,l,itr)= rmaxsz(itr)*tracerpp(k,l,itr)* &
                    EXP(-0.3*zt)*hpp(k,l) / 86400.
                DO k = k1s+1, kms
                    zt = zt + (hpp(k-1,l)+hpp(k,l))/2.
                    sourcesink(k,l,itr)= rmaxsz(itr)*tracerpp(k,l,itr)* &
                        EXP(-0.3*zt)*hpp(k,l) / 86400.
                ENDDO

            ENDDO ! End loop over wet columns

        ENDDO ! End loop over tracers

    END SUBROUTINE srcsnkSZ

    !************************************************************
    SUBROUTINE WQinput
        !***********************************************************
        !
        !            Purpose: To read wq_inp.txt
        !
        !------------------------------------------------------------


        !. . . .Local Variables
        CHARACTER(LEN=12):: wq_input_file="si3d_wq.txt"
        INTEGER::		ios

        !.....Open input parameter file.....
        OPEN (UNIT=i99, FILE=wq_input_file, STATUS="OLD", IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//wq_input_file, ios )

        !.....Read header record containing comments about run................
        READ (UNIT=i99, FMT='(/(A))', IOSTAT=ios) title
        IF (ios /= 0) CALL input_error ( ios, 91)

        !. . Read list of tracerpps modeled
        READ (UNIT=i99,FMT='(///(14X,I20))',IOSTAT=ios) iDO, iARB, iPON, iDON, &
            &    iNH4, iNO3, iPOP, iDOP, iPO4, iALG, iDOM, iPOM, iSOD
        IF (ios /= 0) CALL input_error ( ios, 92)

        !. . . Read model control switches
        READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) Vsod, sedkmod, ktmod
        IF (ios /= 0) CALL input_error ( ios, 101)

        !. . . Read model stochiometeric constants
        READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) anc, acc, apc, roc, &
            &    ron, KDOM, KNIT, KSN, KSP, FNH4
        IF (ios /= 0) CALL input_error ( ios, 93)

        !. . . Read model rates
        READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) sedk, wodk, Gmax, k_r, k_m,   &
            &    k_z, k_ZZ, Z_phy, k_hc, k_hn, k_hp, k_mn, k_mor, k_mp, k_n ,   &
            &    k_ra, k_rs, k_set, k_setarb, k_vn, mu_max
        IF (ios /= 0) CALL input_error ( ios, 94)
  
!        print*, k_a, k_arb, k_dn, k_DOM, k_ex
!        print*, k_gr, k_hc, k_hn, k_hp, k_mn
!        print*, k_mor, k_mp, k_n, k_ra, k_rs
!        print*, k_set, k_setarb, k_vn, mu_max
  
        !... Convert model rates to time of [1/sec]
!        k_a   = k_a/3600.0
!        k_arb = k_arb/3600.0
!        k_dn  = k_dn/3600.0
!        k_DOM = k_DOM/3600.0
!        k_ex  = k_ex/3600.0
!        k_gr  = k_gr/3600.0
!        k_hc  = k_hc/3600.0
!        k_hn  = k_hn/3600.0
!        k_hp  = k_hp/3600.0
!        k_mn  = k_mn/3600.0
!        k_mor = k_mor/3600.0
!        k_mp  = k_mp/3600.0
!        k_n   = k_n/3600.0
!        k_ra  = k_ra/3600.0
!        k_rs  = k_rs/3600.0
!        k_set = k_set/3600.0
!        k_setarb = k_setarb/3600.0
!        k_vn   = k_vn/3600.0
!        mu_max = mu_max/3600.0
  

        !. . . Read model temperature rates
        READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) Theta_a, Theta_dn , &
            &    Theta_DOM, Theta_gr, Theta_hc , Theta_hn, Theta_hp , Theta_mn, Theta_mor, &
            &    Theta_mp , Theta_mu, Theta_n  , Theta_PON, Theta_ra, Theta_SOD, &
            &    Theta_vn
        IF (ios /= 0) CALL input_error ( ios, 95)

        !. . . Read miscillaneous rates
        READ (UNIT=i99,FMT='(///(14X,G20.2))',IOSTAT=ios) ATM_DON, ATM_NH4,    &
            &    ATM_NO3, ATM_PO4, ATM_POP, GW_NH4, GW_NO3, GW_PO4, J_NH4, J_NO3, &
            &    J_PO4
        IF (ios /= 0) CALL input_error ( ios, 96)

        IF (idbg == 1) THEN
            PRINT*, "iARB = ", iARB, "iDO  = ", iDO , "iPON = ", iPON
            PRINT*, "iDON = ", iDON, "iNH4 = ", iNH4, "iNO3 = ", iNO3
            PRINT*, "iPOP = ", iPOP, "iDOP = ", iDOP, "iALG = ", iALG
        END IF

    END SUBROUTINE WQinput

    !************************************************************************
    SUBROUTINE WQinit
        !************************************************************************
        !
        ! Purpose: initializes water quality variables
        !
        !-----------------------------------------------------------------------

        !. . . Local Variables
        INTEGER, DIMENSION(ntrmax):: tracerpplocal
        INTEGER:: i,j, sumtr

        ! ... Initialize tracerpp local
        tracerpplocal = 0

        !. . Initialize constituent locations
        LDO =0; LARB=0; LPON=0; LDON=0;
        LNH4=0; LNO3=0; LPOP=0; LDOP=0; LPO4=0
        LALG=0; LDOM=0; LPOM=0; LSOD= 0

        !. . Assign Lxx to each constituent modeled
        !. . .. first need to define intermediate array tracerpplocal
        i=1
        IF (iDO==1) THEN
            tracerpplocal(i) = 1
            i = i+1
        END IF

        IF (iARB == 1) THEN         ! always for arbitrary constituent
            tracerpplocal(i) = 2
            i = i+1
        END IF

        IF (iPON == 1) THEN
            tracerpplocal(i) = 3
            i = i+1
        END IF

        IF (iDON == 1) THEN
            tracerpplocal(i) = 4
            i = i+1
        END IF

        IF (iNH4 == 1) THEN
            tracerpplocal(i) = 5
            i = i+1
        END IF

        IF (iNO3 == 1) THEN
            tracerpplocal(i) = 6
            i = i+1
        END IF

        IF (iPOP == 1) THEN
            tracerpplocal(i) = 7
            i = i+1
        END IF

        IF (iDOP == 1) THEN
            tracerpplocal(i) = 8
            i = i+1
        END IF

        IF (iPO4 == 1) THEN
            tracerpplocal(i) = 9
            i = i+1
        END IF

        IF (iALG == 1) THEN
            tracerpplocal(i) = 10
            i = i+1
        END IF

        IF (iDOM == 1) THEN
            tracerpplocal(i) = 11
            i = i+1
        END IF

        IF (iPOM == 1) THEN
            tracerpplocal(i) = 12
            i = i+1
        END IF

        IF (iSOD ==1) THEN
            tracerpplocal(i) = 13
        END IF

        sumtr = 0
        DO j = 1, ntrmax
            IF (tracerpplocal(j)>0) sumtr = sumtr + 1
        END DO

        IF (sumtr .ne. ntr) THEN
            print*,('****************************************')
            print*,('               WARNING                  ')
            print*,('    NTR AND NO. tracerpp DO NOT MATCH     ')
            print*,('           EXITING MODEL                ')
            print*,('****************************************')
            STOP
        END IF

        IF (idbg == 1) THEN
            PRINT*, "ntr = ", ntr
            PRINT*, "sumtr = ", sumtr
        END IF
 
        !. . Next define Lxx
        DO i = 1, ntr
            IF (tracerpplocal(i) == 1) THEN
                LDO  = i
            ELSEIF (tracerpplocal(i) == 2) THEN
                LARB = i
            ELSEIF (tracerpplocal(i) == 3) THEN
                LPON = i
            ELSEIF (tracerpplocal(i) == 4) THEN
                LDON = i
            ELSEIF (tracerpplocal(i) == 5) THEN
                LNH4 = i
            ELSEIF (tracerpplocal(i) == 6) THEN
                LNO3 = i
            ELSEIF (tracerpplocal(i) == 7) THEN
                LPOP = i
            ELSEIF (tracerpplocal(i) == 8) THEN
                LDOP = i
            ELSEIF (tracerpplocal(i) == 9) THEN
                LPO4 = i
            ELSEIF (tracerpplocal(i) == 10) THEN
                LALG = i
            ELSEIF (tracerpplocal(i) == 11) THEN
                LDOM = i
            ELSEIF (tracerpplocal(i) == 12) THEN
                LPOM = i
            ELSEIF (tracerpplocal(i) == 13) THEN
                LSOD = i
            END IF
        END DO

        IF (idbg == 1) THEN
            PRINT*, "LDO  = ", LDO , "LARB = ", LARB, "LPON = ", LPON
            PRINT*, "LDON = ", LDON, "LNH4 = ", LNH4, "LNO3 = ", LNO3
            PRINT*, "LPOP = ", LPOP, "LDOP = ", LDOP, "LALG = ", LALG
        END IF
    END SUBROUTINE WQinit
   
    !************************************************************
    SUBROUTINE srcsnkWQ
        !***********************************************************
        !
        !   Purpose: to call all source subroutines for all cells
        !
        !------------------------------------------------------------

        !... Local variables
        INTEGER:: i,j,k,l, k1s, kms,itr

        ! reset soursesink = 0
        sourcesink = 0;

        DO l = 1, lm;

            ! ... Map l- into (i,j)-indexes .........................
            i = l2i(l); j = l2j(l);
	 
            ! ... Retrieve top & bottom wet sal-pts .................
            kms = kmz(i,j)
            k1s = k1z(i,j)
 
            DO k = k1s, kms;

                IF (iARB == 1) THEN
                    CALL sourceARB(k,l)
                END IF
                IF (iDO == 1) THEN
                    CALL sourceDO(k,l,kms,k1s)
                END IF
                IF (iPON == 1) THEN
                    CALL sourcePON(k,l)
                END IF
                IF (iDON == 1) THEN
                    CALL sourceDON(k,l)
                END IF
                IF (iNH4 == 1) THEN
                    CALL sourceNH4(k,l)
                END IF
                IF (iNO3 == 1) THEN
                    CALL sourceNO3(k,l)
                END IF
                IF (iPOP == 1) THEN
                    CALL sourcePOP(k,l)
                END IF
                IF (iDOP == 1) THEN
                    CALL sourceDOP(k,l)
                END IF
                IF (iPO4 == 1) THEN
                    CALL sourcePO4(k,l)
                END IF
                IF (iALG == 1) THEN
                    CALL sourceALG(k,l)
                END IF
                IF (iDOM == 1) THEN
                    CALL sourceDOM(k,l)
                END IF
                IF (iPOM == 1) THEN
                    CALL sourcePOM(k,l)
                END IF
      
            END DO
        END DO

    END SUBROUTINE srcsnkWQ

    !*********************************************************************
    SUBROUTINE sourceARB(kwq,lwq)
        !********************************************************************
        !
        ! Purpose: if arbitrary constituent is modeled, this subroutine
        !  calculates source and sink terms that depend on arbitrary
        !  constituent concentrations
        !
        !  for arbitrary constitucent, sink terms include decay and settling
        !  Chris - chla modelling 13/01/2016 (chla in mg/m3
        !-----------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        REAL    :: theta_T      ! Temperature multiplier
        REAL    :: pf           ! Preference factor for grazing of zooplankton on phytoplankton
        REAL    :: I_s          ! Saturation light intensity of phytoplankton (uE/m2/s)
        REAL    :: I_i          ! (uE/m2/s)
        REAL    :: IP_i         ! average internal P of algae
        REAL    :: IP_min       ! minimum internal P
        REAL    :: IN_i         ! average internal N of algae
        REAL    :: IN_min       ! minimum internal N
        REAL    :: fI           ! light limit
        REAL    :: fIP          ! Phosphorous limit
        REAL    :: fIN          ! N limit
        REAL    :: LNlimit      ! light and nutrient limitations



        I_s     = 125.0;
        I_i     = Qsw*QswFr(kwq,lwq)/0.219; ! sunlight conversion 0.219 photons to W/m2
        IP_min  = 0.3; ! mg P / mg Chla Averaged from Jorgensen et al. 1986
        IN_min  = 2.6; ! mg N / mg Chla Averaged from Jorgensen et al. 1986
        theta_T = 1.02;
        pf      = 0.1;

        IP_i   = 2;  ! mg P / mg Chla
        IN_i   = 11.5;

        !Z_phy  = (tracer(kwq,lwq,LARB)-1.58)/4.97*1000; ! mg/m3 Desortova (2007)


        ! ... estimate light and nutrient limitation
        fI  = I_i / I_s * EXP (1 - I_i / I_s);
        fIP = (IP_i - IP_min) / IP_i;
        fIN = (IN_i - IN_min) / IN_i;

        ! set limit = light limit ! Chris mod 21 MAR 17
        !LNlimit = MIN (fI, fIP, fIN);
        LNlimit = fI;



        sourcesink(kwq,lwq,LARB) = Gmax * theta_T * tracer(kwq,lwq,LARB) * LNlimit &
            & - k_r * theta_T * tracer(kwq,lwq,LARB) - k_m * theta_T * tracer(kwq,lwq,LARB) &
            & - k_z * Z_phy * theta_T * pf * tracer(kwq,lwq,LARB)/(k_ZZ + tracer(kwq,lwq,LARB));

        sourcesink(kwq,lwq,LARB)=sourcesink(kwq,lwq,LARB)/86400; ! Growth rate - per day - Chris 22/03/17

        ! return 0 when nan
        IF (sourcesink(kwq,lwq,LARB) .NE. sourcesink(kwq,lwq,LARB)) THEN
           sourcesink(kwq,lwq,LARB) = 0;
        ENDIF

        ! ... print calculation
        IF (idbg == 1) PRINT *, 'kms = ', kwq, '; ', 'chla_src = ', sourcesink(kwq,lwq,LARB)



    END SUBROUTINE sourceARB

    !*********************************************************************
    SUBROUTINE sourceDO(kwq,lwq,kms,k1s)
        !********************************************************************
        !
        !  Purpose: if dissolved oxygen is modeled, this subroutine
        !  calculates source and sink terms that depend on dissolved
        !  oxygen concentrations
        !
        !-----------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq,kms,k1s

        ! Variables used for surface oxygen tranfer added by KEVIN
        REAL    :: Jatm       ! Oxygen flux at air-water interface (g/m2/s)
        REAL    :: ksurf      ! Oxygen transfer velocity at air-water interface (m/s)
        REAL    :: eps_surf   ! Dissipation rate in surface layer
        REAL    :: Sc_surf    ! Schmidt number in surface layer
        REAL    :: DOsat      ! Saturation concentration for DO in water (mg/L)
        REAL    :: lnDOsat    ! ln of saturation concentration for DO in water
        REAL    :: D_surf     ! Diffusivity of DO in water at surface temperature (cm2/s)
        REAL    :: epsu       ! Dissipation rate in surface layer due to wind shear (W/kg)
        REAL    :: epsw       ! Dissipation rate in surface layer due to convection (W/kg)
        REAL    :: tauo       ! Total shear stress (N/m2)
        REAL    :: taut       ! Tangental shear stress (N/m2)
        REAL    :: Ke         ! Keulegan number
        REAL    :: ustar_surf ! friction velocity in water in surface cell (m/s)
        REAL    :: Uten       ! Wind velocity at 10m above water surface (m/s)
        REAL    :: Jbo        ! Surface buoyancy flux due to heat flux @ water surface, W/kg (Physics & Chemistry of Lakes, 1995)
        REAL    :: atherm     ! thermal expansion coefficient of water (K^-1)
        REAL    :: sedkt      ! Rate constant for sed kinetics adjusted for temperature
        REAL    :: wodkt      ! Water Oxygen Demand constant adjusted for temperature
        REAL    :: drhosurf, DO2t, T_Kb, mu3
        INTEGER :: nwlayers, i, j

        nwlayers = (kms-k1s)+1
          ! ... 3D-(i,j) indexes for l
        i = l2i(lwq); j = l2j(lwq);

        T_K   = sal(2,lwq)+273.15    ! Convert temperature at surface from Celsius to Kelvin

        ! Calculate saturation concentration of DO in water
        lnDOsat = -173.4292+249.6339*(100/T_K)+143.3483*log(T_K/100)-21.8492*(T_K/100) ! Weiss 1970
        DOsat = 1.427*exp(lnDOsat)                             ! Weiss 1970, mg/L, g/m3

        ! Calculate Schmidt number
        mu2 = 0.00002414*(10.0**(247.8/(T_K-140.)))        ! Dynamic Viscosity (kg/m/s)
        rho = 1000*(1.0-0.000019549*abs(sal(2,lwq)-3.84)**1.68)  ! Density of water (kg/m3)
        kin_visc= mu2/rho                  ! Kinematic viscosity, m2/s
        D_surf  = (Dt1*T_K*mu1)/(298.15*mu2)           ! Diffusivity of DO in water, cm2/s (Stokes-Einstein equation)
        Sc_surf = kin_visc/(D_surf/(100*100))          ! Schmidt number

        ! Calculate turbulence in surface cell due to wind shear
        Uten = sqrt(uair(i,j)**2+vair(i,j)**2)     ! Wind velocity at 10m above air-water interface, m/s
        tauo = cdw(i,j)*Uten**2*rhoair             ! Total wind shear stress, N/m2
        ustar_surf = (tauo/rho)**0.5               ! Friction velocity, m/s (Read 2010)
        Ke   = ustar_surf**3/(g*kin_visc)          ! Keulegan number (Soloviev 2007)
        taut = tauo/(1+(Ke/Kecr))              ! Tangental shear stress, N/m2 (Soloviev 2007)
        epsu = (taut/rho)**2/kin_visc              ! Turbulence due to wind shear, N/m2 (Read 2010; Soloviev 2007)

        ! Calculate turbulence in surface cell due to convection
        ! rhoH20 = 10^3*(1-1.9549*10^-5*abs(T-3.84)^1.68   ! Verburg 2010
        ! drhoH20/dT = -3.2842*10^-2*(T-3.84)^0.68     ! Derivative of above equation, KAB
        drhosurf = (-0.032842)*(sal(2,lwq)-3.84)**0.68       ! Drho/DTemp at surface
        atherm= -(1/rho)*drhosurf              ! thermal expansion coefficient
        Jbo  = -(atherm*g)/(Cp*rho)*Hnet(i,j)          ! Surface buoyancy flux

        IF (Jbo > 0) THEN
            epsw = 0        ! Force epsw to zero, per Read 2010 if positive buoyancy flux
        ELSE
            epsw = -Jbo
        ENDIF

        ! Calculate total turbulence in surface cell, transfer velocity, oxygen flux
        eps_surf=epsu + epsw                       ! Total dissipation rate in surface cell
        ksurf = N_eta*(eps_surf*kin_visc)**0.25*Sc_surf**(-n_surf)     ! Transfer velocity, m/s (Read 2010)
        Jatm  = ksurf*(DOsat-tracer(2,lwq,LDO))               ! Oxygen flux at water surface, g/m2/s


        T_Kb  = sal(kms,lwq)+273.15             ! Convert temperature at SWI from Celsius to Kelvin
        sedkt = sedk*theta_sedk**(sal(kms,lwq)-20.0)    ! Adjust sedk for temperature in bottom cell
        mu3   = 0.00002414*(10.0**(247.8/(T_Kb-140.)))    ! Adjust mu for temperature at SWI
        DO2t  = ((Dt1*T_Kb*mu1)/(298.15*mu3))/10000.0 ! Calculate DO2 at Tswi, m2/s

        ! sediment oxygen demand
        SELECT CASE (Vsod)
            ! zero order
            CASE (0)
                ds(k1s) = ds(k1s) - k4sod + Jatm
            ! first order
            CASE (1:)

                SELECT CASE (sedkmod)
                    CASE (0)
                        k4sod=(-(DO2t*sedkt/86400.0)+SQRT((DO2t*sedkt/86400.0)**2.0+2.0*DO2t*sedkt/86400.0*tracer(kms,lwq,LDO)*kt(i,j)**2.0))/kt(i,j)
                    !k4sod=ABS(k4sod)
                    CASE (1:)
                        k4sod=(kt(i,j)*tracer(kms,lwq,LDO))/(1.0+(kt(i,j)/(DO2t*sedkt/86400.0)**(1./2.)))
                END SELECT

                IF (kt(i,j) == 0) THEN
                    k4sod = 0
                ENDIF
        END SELECT



        SELECT CASE (nwlayers)
            ! for single layer column
            CASE (1)
                IF (kwq==k1s) THEN
                    sourcesink(kwq,lwq,LDO)=-k4sod+Jatm;
                ENDIF

            ! mutilayer column
            CASE (2:)

                IF (kwq==k1s) THEN
                    sourcesink(kwq,lwq,LDO)=Jatm;
                ELSEIF (kwq==kms) THEN
                    sourcesink(kwq,lwq,LDO)=-k4sod;
                ENDIF


        END SELECT



    ! ... older version
    !
    !
    !  !. . . Local Variables
    !  REAL		::	 Tk, lnOS, OS, ln_Pwv, Pwv, theta2, Patm
    !
    !  ! Calculate DO saturation
    !  Tk = salp(kwq,lwq) + 273
    !  lnos = -139.34410 + 1.575701*1E5 /(Tk    ) &
    !  &                 - 6.642308*1E7 /(Tk**2.) &
    !  &	                + 1.243800*1E10/(Tk**3.) &
    !  &                 - 8.621949*1E11/(Tk**4.)
    !  os = EXP(lnos)
    !
    !  ! Correct for Patmospheric (Pa - declared in si3d_types and defined in surfbc0)
    !  Patm   = Pa * 0.00000986923; ! Transform atmospheric pressure from Pa to atm
    !  ln_Pwv = 11.8751 - 3840.70/Tk - 216961/Tk
    !  Pwv    = EXP(ln_Pwv)
    !  theta2 = 0.000975 - 1.426*1E-5 * salp(kwq,lwq) + &
    !                      6.436*1E-8 * salp(kwq,lwq)**2
    !  os = os*Pa*((1-Pwv/Pa) *(1-theta2*Pa))&
    !  &           /((1-Pwv)*(1-theta2) )
    !
    !  ! Calculate reaertaion
    !  ! for now using constant rearation defined in wq_inp, but in future, can have
    !  ! alternatives for calcualting reaeration.
    !
    !  IF (kwq .eq. k1z(l2i(lwq),l2j(lwq))) THEN
    !  sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDON)               +		&
    !  &                         k_a*(tracerpp(kwq,lwq,LDO) - OS) 		!reaeration
    !					! + photosynthesis	- only if IALG = 1; calculated in sourceALG
    !					! - Respiration		- only if IALG = 1; calculated in sourceALG
    !					! - Nitrification	- only if INH4 = 1; calculated in sourceNH4
    !					! - SOD				- only if ISOD = 1; calculated in sourceSOD
    !					! - Oxidation of OM	- only if IDOM = 1; calculated in sourceDOm
    ! END IF
    ! ... end

    END SUBROUTINE sourceDO

    !************************************************************************
    SUBROUTINE sourceDON(kwq,lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on DON
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . . Local variables
        REAL:: minrl

        !. . mineralization
        minrl = k_mn * theta_mn**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDON)

        sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)	-				&
            &                          minrl
        ! + hydrolysis	- only if IPON = 1; caluclated in sourcePON
        ! + Resp/excr	- only if IAlG = 1; calculated in sourceALG

        !. . .Add contribution from atmospheric deposition to top layer
        IF (kwq .eq. k1z(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq,lwq, LDON) = sourcesink(kwq, lwq, LDON) + &
                & hpp(kwq,lwq)*ATM_DON		! Atmoshperic deposition
        END IF
	
        IF (INH4 == 1) THEN		! add mineralization of DON to NH4
            sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) + minrl
        END IF

    END SUBROUTINE sourceDON

    !************************************************************************
    SUBROUTINE sourcePON(kwq,lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on PON
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !... Local variables
        REAL:: hydrol


        !. Calculate hydrolysis
        hydrol = k_hn *theta_hn**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPON)
        sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON)       -	&
            &                          hydrol					        -	&	! hydrolysis
            &                          k_set*tracerpp(kwq,lwq,LPON)	-	&	! settling
            &                          k_rs * tracerpp(kwq,lwq,LPON)			! resusupension
						                                            ! + mortality	- only if IALG = 1; calcualted in sourceALG

        ! Add contribution of mineralization to DON concentration
        IF (IDON == 1) THEN
            sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)	+  hydrol
        END IF

    END SUBROUTINE sourcePON

    !************************************************************************
    SUBROUTINE sourceNH4(kwq,lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on NH4
        !
        !--------------------------------------------------------------------------
        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . . Local Variables
        REAL:: nitrif, f_DO

        !. . . Calculate nitrification
          ! Calculate DO inhibition of nitrification
        IF (IDO == 1) THEN
            f_DO = tracerpp(kwq,lwq,LDO) /(KNIT + tracerpp(kwq,lwq,LDO) )
        ELSE
            f_DO = 1.0
        END IF

        nitrif = k_n*theta_n**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LNH4) * f_DO

        sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4)		-	&
            &  nitrif			-	&	! nitrificatin
            &  k_vn * theta_vn**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LNH4) ! volitilization
						                                            ! - Algal Uptake		- if IALG = 1; calculated in sourceALG
						                                            ! + mineralization		- if IDON = 1; calculated in sourceDON

        !. . Add contribution from sediment release and GW flux into bottom cells
        IF (kwq .eq. kmz(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) +   &
                & hpp(kwq,lwq)*J_NH4     +	 &	! sediment release
                & hpp(kwq,lwq)*GW_NH4	     	! GW flux
        END IF

        !. . Add contribution from atmoshperic deposition into top cells
        IF (kwq .eq. k1z(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq, lwq, LNH4) = sourcesink(kwq, lwq, LNH4) + &
                & hpp(kwq, lwq) * ATM_NH4	 ! ATM deposition
        END IF


        !. . Add contribution from nitrification to nitrate sourcesink

        IF (INO3 == 1) THEN
            sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) + nitrif
        END IF



    END SUBROUTINE sourceNH4

    !************************************************************************
    SUBROUTINE sourceNO3(kwq,lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on NO3
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq


        !sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3)		-	&
        !    &  k_dn*theta_dn**(salp(kwq,lwq) - 20)* tracerpp(kwq,lwq,LNO3)	! denitrification
						                                            ! + nitrification		- if INH4 = 1; calculated in sourceNH4
						                                            ! - algal uptake		- if IALG = 1; calculated in sourceALG

        !. . Add contribution from sediment release and GW flux into bottom cells
        IF (kwq .eq. kmz(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) +   &
                & hpp(kwq,lwq)*J_NO3     +	 &	! sediment release
                & hpp(kwq,lwq)*GW_NO3	     	! GW flux
        END IF

        !. . Add contribution from atmoshperic deposition into top cells
        IF (kwq .eq. k1z(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq, lwq, LNO3) = sourcesink(kwq, lwq, LNO3) + &
                & hpp(kwq, lwq) * ATM_NO3	 ! ATM deposition
        END IF

    END SUBROUTINE sourceNO3

    !************************************************************************
    SUBROUTINE sourcePOP(kwq,lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on POP
        !
        !--------------------------------------------------------------------------
        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . . Local Variables
        REAL:: hydrol

        ! Calculate hydrolysis

        hydrol = k_hp*theta_hp**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPOP)

        sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)				+	&
            &  k_rs*tracerpp(kwq,lwq,LPOP)	-	&	! Resuspension
            &  hydrol						! hydrolysis

        !... Add atmoshperic deposition to top layer
        IF (kwq == k1z(l2i(lwq), l2j(kwq))) THEN
            sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP)				+	&
                &  h(kwq,lwq)*ATM_POP					! Atmospheric deposition
        END IF

        ! Caclulate hydrolysis contribution to DOP
        IF (IDOP == 1) THEN
            sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP) + hydrol
        END IF
						
    END SUBROUTINE sourcePOP

    !************************************************************************
    SUBROUTINE sourceDOP(kwq, lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on DOP
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . .Local Variables
        REAL:: minrl

        ! Calculate mineralization

        minrl = k_mp*theta_mp**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDOP)

        sourcesink(kwq,lwq,LDOP) = sourcesink(kwq,lwq,LDOP)		-	&
            &  minrl				! hydrolysis
						                                            ! + resp/excretion	- if IAlG = 1; calculated in sourceALG
						                                            ! + hydrolysis		- if iPOP = 1; calculated in sourcePOP

        ! Caclulate mineralization contribution to PO4
        IF (IPO4 == 1) THEN
            sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + minrl
        END IF
						
    END SUBROUTINE sourceDOP

    !************************************************************************
    SUBROUTINE sourcePO4(kwq, lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on PO4
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq
 
        sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4)		+	&
            &                          h(kwq,lwq)*J_PO4			+	&	! Sediment Release
            &                          h(kwq,lwq)*GW_PO4			+	&	! GW flux
            &                          h(kwq,lwq)*ATM_PO4				! Atmospheric deposition
           ! + mineralization		- if IDOP = 1; calculated in sourceDOP
           ! + algal resp			- if IALG = 1; calculated in sourceALG
           ! - algal uptake		- if IALG = 1; calculated in sourceALG

        !. . Add contribution from sediment release and GW flux into bottom cells
        IF (kwq .eq. kmz(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) +   &
                & hpp(kwq,lwq)*J_PO4     +	 &	! sediment release
                & hpp(kwq,lwq)*GW_PO4	     	! GW flux
        END IF

        !. . Add contribution from atmoshperic deposition into top cells
        IF (kwq .eq. k1z(l2i(lwq), l2j(lwq))) THEN
            sourcesink(kwq, lwq, LPO4) = sourcesink(kwq, lwq, LPO4) + &
                & hpp(kwq, lwq) * ATM_PO4	 ! ATM deposition
        END IF

						
    END SUBROUTINE sourcePO4

    !************************************************************************
    SUBROUTINE sourcePOM (kwq, lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on POM
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . Local Variables
        REAL:: hydrol

        ! Calculate hydrolysis
        hydrol = k_hc*theta_hc**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LPOM)
        sourcesink(kwq,lwq,LPOM) = sourcesink(kwq,lwq,LPOM)				-	&
            &  k_set * tracerpp(kwq,lwq,LPOM) 	-	&	! settling
            &  hydrol							! hydrolysis
						                                            ! + mortality		- if IALG = 1; calculated in sourceALG

        ! Caclulate hydrolysis contribution to DOM
        IF (IDOM == 1) THEN
            sourcesink(kwq,lwq,LDOM) = sourcesink(kwq,lwq,LDOM) + hydrol
        END IF
						
    END SUBROUTINE sourcePOM

    !************************************************************************
    SUBROUTINE sourceDOM(kwq, lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on DOM
        !
        !--------------------------------------------------------------------------

        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . .Local Variables
        REAL:: F_DO, oxid

        !. . Calculate how DO concn impedes oxidation of DOM
        IF (IDO ==1) THEN
            F_DO = (tracerpp(kwq,lwq,LDO))/(kDOM + tracerpp(kwq,lwq,LDO) )
        END IF

        ! Calculate oxidation
        !oxid = K_DOM*theta_DOM**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LDOM)*F_DO

        sourcesink(kwq,lwq,LDOM) = sourcesink(kwq,lwq,LDOM) - oxid
						                                            ! + hydrolysis	-if IPOM = 1; calculated in sourcePOM

        IF (IDO == 1) THEN
            sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) - roc*oxid
        END IF
											
    END SUBROUTINE sourceDOM

    !************************************************************************
    SUBROUTINE sourceALG(kwq, lwq)
        !*********************************************************************
        !
        ! Purpose: To calculate sourcesink terms that depend on ALG
        !
        !--------------------------------------------------------------------------
        ! ... Arguments
        INTEGER, INTENT (IN) :: kwq,lwq

        !. . .Local Variables
        REAL::	mu, f_L, f_T, f_N, f_P, N_conc
        REAL::	resp, excr, mort, growth

        ! Calculate mu, growth rate
	
        !. .  Calculate growth limiting factors
        ! light limitation
        !		f_L = SolarFR(depth)*0.47/light_sat *e**(-SolarFR(depth)*0.47/light_sat +1)
        ! f_L construct needs to be fixed, need more information.  for now, use f_L = 1
        f_L = 1.0

        ! temperature limitaton
        f_T = theta_mu**(salp(kwq,lwq -20))

        ! nutrient limitation - but only if the nutrients are modeled
        IF ((INH4 ==1) .AND. (INO3 ==1)) THEN
            N_conc = tracerpp(kwq,lwq,LNH4) + tracerpp(kwq,lwq,LNO3)
            f_N = N_conc/(KSN + N_conc)
        ELSE IF (INH4 == 1 .AND. INO3 ==0) THEN
            N_conc = tracerpp(kwq,lwq,LNH4)
            f_N = N_conc/(KSN + N_conc)
        ELSE IF (INH4 == 0 .AND. INO3 ==1) THEN
            N_conc = tracerpp(kwq,lwq,LNO3)
            f_N = N_conc/(KSN + N_conc)
        ELSE
            f_N = 1.0
        END IF

        IF (IPO4 ==1) THEN
            f_P = tracerpp(kwq,lwq,LPO4) /(KSP + tracerpp(kwq,lwq,LPO4) )
        ELSE
            f_P = 1.0
        END IF

        mu = mu_max*f_L*f_T*f_N*f_P

        !. . Calculate growth
        growth = mu_max * tracerpp(kwq,lwq,LALG)
        !. . Calculate respiration
        resp   = k_ra*theta_ra**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG)
        !. . Calculate excretion
        !excr   = k_ex*theta_ra**(salp(kwq,lwq) - 20) * tracerpp(kwq,lwq,LALG)
        !. . Calculate mortality
        mort   = k_mor*theta_mor**(salp(kwq,lwq) - 20)*tracerpp(kwq,lwq,LALG)


        !sourcesink(kwq,lwq,LALG) = sourcesink(kwq,lwq,LALG)	+									&
        !    &  growth		-										& ! growth
        !    &  resp		-										& ! respiration
        !    &  excr		-										& ! excretion
        !    &  k_gr*theta_gr**(salp(kwq,lwq)-20) * tracerpp(kwq,lwq,LALG) -	& ! grazing
        !    &  k_set*tracerpp(kwq,lwq,LALG)								 ! settling

        ! If dissolved oxygen is modeled, alter sourcesink(kwq,lwq,LDO) to reflect growth and resp
        ! of algae population
        IF (IDO ==1) THEN
            sourcesink(kwq,lwq,LDO) = sourcesink(kwq,lwq,LDO) + acc*roc*(growth-resp)
        END IF

        ! IF PON is modeled, alter sourcesink(kwq,lwq,LPON) to reflect mortality of algae
        IF (IPON == 1) THEN
            sourcesink(kwq,lwq,LPON) = sourcesink(kwq,lwq,LPON) + anc*mort
        END IF

        ! If DON is modeled, alter soucesink(kwq,lwq,LDON) to reflect excretion of algae
        IF (IDON ==1) THEN
            sourcesink(kwq,lwq,LDON) = sourcesink(kwq,lwq,LDON)  + anc*(resp + excr)
        END IF

        ! IF NH4 is modeled, alter sourcesink(kwq,lwq,LNH4) to include uptake of NH4 by algae
        IF (INH4 == 1) THEN
            sourcesink(kwq,lwq,LNH4) = sourcesink(kwq,lwq,LNH4) - anc*FNH4*growth
        END IF

        ! IF NO3 is modeled, alter sourcesink(kwq,lwq,NO3) to include uptake of NO3 by algae
        IF (INO3 == 1) THEN
            sourcesink(kwq,lwq,LNO3) = sourcesink(kwq,lwq,LNO3) - anc*(1-FNH4) * growth
        END IF

        ! If POP is modeled, alter sourcesink(kwq,lwq,POP) to include mortality of algae
        IF (IPOP == 1) THEN
            sourcesink(kwq,lwq,LPOP) = sourcesink(kwq,lwq,LPOP) + apc*mort
        END IF

        ! If PO4 is modeled, alter sourcesink(kwq,lwq,LPO4) to include resp and uptake
        IF (IPO4 == 1) THEN
            sourcesink(kwq,lwq,LPO4) = sourcesink(kwq,lwq,LPO4) + apc*(resp - growth)
        END IF

        ! If DOM is modeled, alter sourcesink(kwq,lwq,LDOM) to include respiration
        IF (IDOM == 1) THEN
            sourcesink(kwq,lwq,LDOM) = sourcesink(kwq,lwq,LDOM) + resp + excr
        END IF

    END SUBROUTINE sourceALG

END MODULE si3d_ecomod

