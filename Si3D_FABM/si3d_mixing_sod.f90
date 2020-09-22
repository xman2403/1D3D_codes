!************************************************************************
MODULE si3d_mixing
    !************************************************************************
    !
    !  Purpose: Define vertical mixing coefficients for momentum and mass
    !           Includes the schemes originally coded by Peter Smith and
    !           the 2 Eqs. model proposed by Mellor&Yamada (level 2.5) and later
    !           modified by Kantha & Clayson.
    !
    !-------------------------------------------------------------------------

    USE si3d_types
    USE turbulence,  only: turb_method				! si3dgotm
    USE turbulence,  only: init_turbulence, do_turbulence		! si3dgotm
    USE turbulence,  only: num,nuh,nus				! si3dgotm
    USE turbulence,  only: eps, L, tke, tkeo			! si3dgotm
    USE turbulence,  only: const_num,const_nuh			! si3dgotm
    USE turbulence,  only: gamu,gamv,gamh,gams			! si3dgotm
    USE turbulence,  only: kappa					! si3dgotm
    USE kpp,         only: init_kpp,do_kpp			! si3dgotm
    USE mtridiagonal,only: init_tridiagonal,clean_tridiagonal	! si3dgotm
    USE eqstate,     only: init_eqstate				! si3dgotm

    IMPLICIT NONE
    SAVE

CONTAINS


    !***********************************************************************
    SUBROUTINE UpdateMixingCoefficients
        !***********************************************************************
        !
        !  Purpose: To compute the vertical distribution of the vertical
        !           eddy coefficients at all horizontal nodal points. The
        !           choice of a turbulence model depends on the parameter
        !           `iturb'.
        !
        !             if iturb = 0 --> eddy coefficients are constants
        !             if iturb = 1 --> eddy viscosities & diffusivities calculated
        !                              from TKE & a length scale using the
        !                              quasi-equilibrium turbulent energy model of
        !                              Galperin et. al 81988) (Mellor Yamada 2.5)
        !             it iturb = -1--> GOTM libraries
        !
        !           Note: The vertical eddy coefficients are defined at the
        !                 interfaces located vertically between pressure-pts.
        !                 Later in subroutine  matmom  they are interpolated to
        !                 the interfaces located vertically between the u-pts
        !                 and the v-pts where they are used in the computations.
        !                 The index  k  for a particular layer references the
        !                 eddy coefficients at the upper interface of the layer.
        !
        !-----------------------------------------------------------------------

        !.....Local variables.....
        INTEGER :: i, j, k, ll, kms, k1s, nwlayers, m
        REAL (real_G1), DIMENSION(1 :km1) :: drhodz, rhoavg, N2, xl, deltaz
        REAL (real_G1) :: Sm, Sh, Gh, Rt
        !REAL, DIMENSION(1:km1) :: dudz, dvdz, ssu, ssv
        !REAL, DIMENSION(0:km ) :: NN, MM, hl
        !REAL  :: depth, taub, ustarb, taus, ustars, z0s, z0b, usurf, vsurf, &
        !         uup, uupo, udn, udno, vup, vupo, vdn, vdno
        double precision :: depth, ustars, ustarb, z0s, z0b, copydt  ! si3dgotm
        REAL, DIMENSION(1:km1) :: dudz, dvdz, ssu, ssv               ! si3dgotm
        double precision, DIMENSION(0:km ) :: NN, MM, hl             ! si3dgotm
        REAL  :: taub, taus, usurf, vsurf, &                         ! si3dgotm
            uup, uupo, udn, udno, vup, vupo, vdn, vdno          ! si3dgotm
        REAL, PARAMETER :: charnock_val = 1400.0, z0s_min = 0.02, h0b = 0.05
        REAL  :: tint, kab ! Integral time scale (seconds)
        !.....Timing.....
        REAL, EXTERNAL :: TIMER
        REAL :: btime, etime
        REAL    :: Dff_wtr     ! Diffusivity of DO in water (cm2/s)

        ! ... Time counter & No. of calls to turbulence submodel
        btime = TIMER(0.0)
        n_turb = n_turb + 1
        copydt = dt ! si3dgotm

        SELECT CASE (iturb)

            !          ----- START Inteface si3dgotm -------
            CASE(-1)

                ! ... Solve only on the last iteration in each time step
                IF ( lastiter < 1) RETURN

                !.....Sweep over wet pressure points
                DO ll = 1, lm

                    ! .... Map 3D-(i,j) and l- index
                    i = l2i(ll); j = l2j(ll);

                    !.....Compute the number of wet layers.....
                    kms = kmz(i,j); k1s = k1z(i,j)
                    nwlayers = (kms-k1s) + 1

                    ! ... Ignore 1-layer colums
                    IF (nwlayers <= 1) CYCLE

                    ! ... Initialize & calculate MM & NN frequencies
                    ! Note that MM, NN, and hl are ordered from bottom to top
                    ! i.e. bottom layer (km) = 1 and top layer = nwlayers
                    ! following Gotm convention
                    MM = 0.0; NN = 0.0; m = 0; depth = 0; hl = 0.0;
                    DO k = kms, k1s+1,-1

                        ! ... Update counter
                        m = m + 1

                        ! ... Depth & layer thikness
                        depth = depth + hp(k,ll)
                        hl(m) = hp(k,ll)

                        ! ... Calculate buoyancy frequency MM
                        deltaz(k) = 0.5 * (hp  (k-1,ll) + hp  (k,ll))
                        rhoavg(k) = 0.5 * (rhop(k-1,ll) + rhop(k,ll)) + 1000.
                        drhodz(k) = (rhop(k-1,ll) - rhop(k,ll))/deltaz(k)
                        NN    (m) = - g * drhodz(k) / rhoavg(k)
         
                        ! ... Evaluate shear frequency MM
                        !dudz  (k) = 0.5*(up(k-1,ll)+up(k-1,lWC(ll))-   &
                        !&                up(k  ,ll)-up(k  ,lWC(ll)))/deltaz(k)
                        !dvdz  (k) = 0.5*(vp(k-1,ll)+vp(k-1,lSC(ll))-   &
                        !&                vp(k  ,ll)-vp(k  ,lSC(ll)))/deltaz(k)
                        !MM    (m) =  dudz(k) **2. + dvdz(k)**2.
         
                        ! ... The following lines follow GOTM meanflow routines
                        uup  = 0.5*(u  (k-1,ll)+u  (k-1,lWC(ll)))
                        udn  = 0.5*(u  (k  ,ll)+u  (k  ,lWC(ll)))
                        uupo = 0.5*(upp(k-1,ll)+upp(k-1,lWC(ll)))
                        udno = 0.5*(upp(k  ,ll)+upp(k  ,lWC(ll)))
                        ssu(k)= 0.5*(                                      &
                            (theta    *(uup -udn )*(uup -udno)+    &
                            (1.-theta)*(uupo-udno)*(uupo-udn ))    &
                            /deltaz(k) / hp(k,ll)                  &
                            +(theta    *(uup -udn )*(uupo-udn )+    &
                            (1.-theta)*(uupo-udno)*(uup -udno))    &
                            /deltaz(k) / hp(k,ll)                  )

                        vup  = 0.5*(v  (k-1,ll)+v  (k-1,lSC(ll)))
                        vdn  = 0.5*(v  (k  ,ll)+v  (k  ,lSC(ll)))
                        vupo = 0.5*(vpp(k-1,ll)+vpp(k-1,lSC(ll)))
                        vdno = 0.5*(vpp(k  ,ll)+vpp(k  ,lSC(ll)))
                        ssv(k)= 0.5*(                                      &
                            (theta    *(vup -vdn )*(vup -vdno)+    &
                            (1.-theta)*(vupo-vdno)*(vupo-vdn ))    &
                            /deltaz(k) / hp(k,ll)                  &
                            +(theta    *(vup -vdn )*(vupo-vdn )+    &
                            (1.-theta)*(vupo-vdno)*(vup -vdno))    &
                            /deltaz(k) / hp(k,ll)                  )
                        MM(m) = ssu(k) + ssv(k)
                    END DO

                    ! ... Add conttribution of the top-most layer
                    hl(m+1) = hp(k1s,ll)
                    depth   = depth + hp(k1s,ll)

                    ! ... Bottom friction velocity .......................
                    taub  = cd  * ( (0.5*(u(kms,ll) + u(kms,lWC(ll))))**2.&
                        &            +  (0.5*(v(kms,ll) + v(kms,lSC(ll))))**2.)
                    ustarb = sqrt(taub)

                    ! ... Surface friction velocity ......................
                    usurf = 0.5*(up(k1u(i,j),ll) + up(k1u(i-1,j),lWC(ll)) )
                    vsurf = 0.5*(vp(k1v(i,j),ll) + vp(k1v(i,j-1),lSC(ll)) )
                    taus  = cdw(i,j)*rhoair * ((uair(i,j)-usurf)**2.+     &
                        &                          (vair(i,j)-vsurf)**2.)
                    ustars = sqrt(taus/1000.)

                    ! ... Calculate bottom roughness z0b
                    z0b = 0.03*h0b

                    ! ... Calculate surface roughness z0s
                    z0s = MAX(charnock_val * (ustars**2.)/g, z0s_min)

                    ! ... Set tke, eps, L & diff. to stored values in si3d
                    m = 0
                    DO k = kms, k1s+1,-1
                        m = m + 1
                        L   (m) = si3dlen(k,ll)
                        eps (m) = si3deps(k,ll)
                        tke (m) = si3dtke(k,ll)
                        tkeo(m) = si3dtke(k,ll)
                        num (m) = Av     (k,ll) - AvMolecular
                        nuh (m) = Dv     (k,ll) - DvMolecular
                    ENDDO

                    ! ... Udate tke, eps, L & diff.
                    CALL do_turbulence(nwlayers,copydt,depth,ustars,ustarb,z0s,z0b,  &
                        hl(0:nwlayers),NN(0:nwlayers),MM(0:nwlayers))
 
                    ! ... Save updated tke, eps, L & diff. in si3d variables
                    m = 0
                    DO k = kms, k1s+1,-1
                        m = m + 1
                        si3dlen(k,ll) = L   (m)
                        si3deps(k,ll) = eps (m)
                        si3dtke(k,ll) = tke (m)
                        Av     (k,ll) = num (m) + AvMolecular
                        Dv     (k,ll) = nuh (m) + DvMolecular
                    ENDDO

                    ! ... Assign transfer coefficients at top & bottom surfaces
                    Av(1:k1s,ll) = 0.0; Av(kms+1:km1,ll) = 0.0
                    Dv(1:k1s,ll) = 0.0; Dv(kms+1:km1,ll) = 0.0

                ENDDO

            !          ----- END Inteface si3dgotm -------

            !          ----- Method 0: Constant eddy coefficients -------
            CASE (0)

                ! ... Eddy viscosity & diffusivity already initialized
                RETURN

            !          ----- Method 1: 2-equation level coefficients -----
            CASE(1:)

                ! ... Solve for q2 and q2l for next time step
                CALL SolveTKEandLengthScale

                !.....Sweep over wet pressure points
                DO ll = 1, lm

                    ! .... Map 3D-(i,j) from 2L-l index
                    i = l2i(ll); j = l2j(ll);

                    !.....Compute the number of wet layers.....
                    kms = kmz(i,j); k1s = k1z(i,j)
                    nwlayers = (kms-k1s) + 1

                    ! ... Ignore 1-layer colums
                    IF (nwlayers <= 1) CYCLE

                    DO k = k1s+1, kms

                        ! ... Calculate density gradients @ interfaces
                        deltaz(k) = 0.5 * (hp  (k-1,ll) + hp  (k,ll))
                        rhoavg(k) = 0.5 * (rhop(k-1,ll) + rhop(k,ll)) + 1000.
                        drhodz(k) = (rhop(k-1,ll) - rhop(k,ll))/deltaz(k)
                        N2(k)     = - g * drhodz(k) / rhoavg(k)

                        ! ... Calculate length scale at time n
                        xl(k)    = q2lp (k,ll) / q2p (k,ll)
  
                        ! ... Length scale limitation from Galperin et al (1988)
                        IF(N2(k)>0.0) xl(k)=MIN(xl(k),0.53*DSQRT(q2p(k,ll)/N2(k)))
         
                        ! ... Evaluate transfer coefficients at time n
                        Gh = - N2(k) * xl(k)**2 / q2p(k,ll)
                        Gh = MAX( Gh_min_Kc,MIN(Gh_max_Kc, Gh) )
                        Sh = t_1/(1.- Gh * t_2)
                        Sm = (t_3 + Sh * Gh * t_4 )/(1.-t_5*Gh)
                        Rt = DSQRT (q2p(k,ll)) * xl(k)
                        Av(k,ll) =  Rt * Sm +  AvMolecular
                        Dv(k,ll) =  Rt * Sh +  DvMolecular

                        ! ... Evaluate kt ...
                        IF (k == kms) THEN
                            T_K   = sal(k,ll)+273.15    ! Convert temperature at surface from Celsius to Kelvin
                            mu2 = 0.00002414*(10.0**(247.8/(T_K-140.)))     ! Dynamic Viscosity (kg/m/s)
                            rho = 1000*(1.0-0.000019549*abs(sal(k,ll)-3.84)**1.68)  ! Density of water (kg/m3)
                            kin_visc= (mu2/rho)*(100.**2)           ! Kinematic viscosity, cm2/s
                            Dff_wtr  = (Dt1*T_K*mu1)/(298.15*mu2)            ! Diffusivity of DO in water, cm2/s (Stokes-Einstein equation)

                            SELECT CASE (ktmod)
                                !Lorke & Peeters
                                CASE (1)
                                    kt(i,j) = (1/(2*pi))*(kin_visc/Dff_wtr)**((-1.)/2.)*(epskt(i,j)*kin_visc/10000.)**(1./4.) ! kt in m/s, assumes epskt in m2/s3, kin_visc in cm2/s
                                    !                  PRINT *, ' Av, Dv, epskt', Av(k,ll), Dv(k,ll), epskt(i,j)
                                    ! Large Eddy
                                CASE (2:)
                                    !                  kab=u(k,ll)**2. + v(k,ll)**2.
                                    !                  tint = xl(k)/(SQRT(kab)) ! xl in m, u&v in m/s, tint in s
                                    !kt(i,j) = DSQRT((4*Dv(k,ll))/(xl(k)*pi))
                                    !                  kt(i,j) = SQRT((4*Dff_wtr)/(tint*pi))/100. ! kt in m/s
                                    tint = q2(k,ll)/epskt(i,j)                ! Integral time scale (s)
                                    kt(i,j) = SQRT((4*(Dff_wtr/10000.)*(1./tint))/pi)  ! Tranfer velocity (m/s)
                            END SELECT

                        END IF
                    END DO

                    IF (hp(k1s,ll)<1.E-2) THEN
                        Av(k1s+1,ll) = MAX(Av(k1s+1,ll),1.E-4);
                        Dv(k1s+1,ll) = MAX(Dv(k1s+1,ll),1.E-4);
                    ENDIF

                    ! ... Assign transfer coefficients at top & bottom surfaces
                    Av(1:k1s,ll) = 0.0; Av(kms+1:km1,ll) = 0.0
                    Dv(1:k1s,ll) = 0.0; Dv(kms+1:km1,ll) = 0.0

                END DO

        END SELECT

        !.....Compute CPU time spent in subroutine.....
        etime = TIMER(0.0)
        t_turb = t_turb + (etime - btime)

    END SUBROUTINE UpdateMixingCoefficients

    !***********************************************************************
    SUBROUTINE SolveTKEandLengthScale
        !***********************************************************************
        !
        !  Purpose: Solve a 2 equation model for turbulence quantities (2*TKE &
        !           a generic length scale variable - in this case 2*TKE*l).
        !
        !  Replace l (index) by m so that there ins not any overlap with variable
        !  L (length scale) used in gotm
        !
        !-----------------------------------------------------------------------

        !.....Local variables.....
        REAL (real_G1), DIMENSION(1:km1) ::  zfromb, zfromt, deltaz, deltazp
        REAL (real_G1), DIMENSION(1:km1) ::  drhodz, rhoavg, xl, dvdz, dudz, &
            &                                     ShearP, BuoyP, dcoeff, N2, N2p, rho_c
        REAL (real_G1) :: Sm, Sh, Gh, dum_1, taurb, taurs, zfromb0, zfromt0, Dissipijk
        INTEGER        :: i, j, k, m, istat, nwlayers, kms, k1s, nn, in, jn
        REAL           :: twodt1, htot, vsurf, usurf, ustar, rhoxh

        ! .. Initialize twodt1
        twodt1 = twodt/istep

        ! ... Initialize variables for basin scale energy balances
        IF ( iobal > 0 .AND. istep > 1) THEN
            ShearProduction    = 0.E0
            BuoyancyProduction = 0.E0
            Dissipation        = 0.E0
            TKinE              = 0.E0
        ENDIF

        !.....Sweep over wet pressure points
        DO m = 1, lm

            ! ... 3D-(i,j) indexes for l (= m)
            i = l2i(m); j = l2j(m);

            !.....Compute the number of wet layers.....
            kms = kmz(i,j); k1s = k1z(i,j);
            nwlayers = (kms-k1s) + 1

            !.....Ignore 1-layer columns
            IF (nwlayers <= 1) CYCLE

            !.....Calculate total depth
            htot = SUM(hp(k1s:kms,m));

            ! Compute the array of vertical distances from the bottom
            ! of the water column to the top of each layer
            zfromb0 = 0.0
            DO k = kms, k1s+1, -1
                zfromb0 = zfromb0 + hp(k,m)
                zfromb(k) = zfromb0
            END DO
            zfromb(k1s  ) = htot
            zfromb(kms+1) = 0.0E0

            ! Compute the array of vertical distances from the top
            ! of the water column to the top of each layer
            zfromt0 = 0.0
            DO k = k1s+1, kms
                zfromt0   = zfromt0 + hp(k-1,m)
                zfromt(k) = zfromt0
            END DO
            zfromt(k1s  ) = 0.0E0
            zfromt(kms+1) = htot

            ! ... Calculate length scale at time n
            xl(k1s+1:kms) = q2lp(k1s+1:kms,m)/ q2p(k1s+1:kms,m)
            xl(k1s)   = 0.0
            xl(kms+1) = 0.0


            ! ----- TURBULENT KINETIC ENERGY ( q2 ) equations -----------------

            ! ... Calculate shear & buoyancy source terms
            DO k = k1s+1, kms

                ! ... Calculate density  gradients @ time n
                deltaz(k) = 0.5 * (hp  (k-1,m) + hp  (k,m))
                rhoavg(k) = 0.5 * (rhop(k-1,m) + rhop(k,m)) + 1000.
                drhodz(k) =       (rhop(k-1,m) - rhop(k,m))/deltaz(k)
                N2(k)     = - g / rhoavg(k) * drhodz (k)

                ! ... Evaluate shear
                dudz(k)  = 0.5*( up(k-1,m)+up(k-1,lWC(m)) -   &
                    &                up(k  ,m)-up(k  ,lWC(m)) ) / deltaz(k)
                dvdz(k)  = 0.5*( vp(k-1,m)+vp(k-1,lSC(m)) -   &
                    &                vp(k  ,m)-vp(k  ,lSC(m)) ) / deltaz(k)

                ! ... Length scale limitation from Galperin et al (1988)
                IF (N2(k)>0.0E0) xl(k)=MIN(xl(k),0.53*DSQRT(q2p(k,m)/N2(k)))
 
                ! ... Calculate production of TKE @ interfaces
                ShearP(k)=  Av(k,m)*( dudz(k) **2. + dvdz(k)**2. )
                BuoyP (k)= -Dv(k,m)*( N2 (k)                     )

            END DO

            ! ... Form diffusion coefficients for q2 & q2l & time n
            dcoeff (k1s:kms) = Sq * 0.50* (xl(k1s:kms  )+ xl(k1s+1:kms+1  )) *  &
                &               DSQRT ( 0.50*(q2p(k1s:kms,m)+q2p(k1s+1:kms+1,m)) )
    
            ! ... Form tridiagonal matrix [aa] ...............................

            ! a. Define upper diagonal terms
            aaT(3,k1s     ) = 0.0
            aaT(3,k1s+1:kms) = - dcoeff(k1s+1:kms) /                            &
                &  (hp(k1s+1:kms,m) * 0.50 * (hp(k1s+1:kms,m)+hp(k1s:kms-1,m)))
            aaT(3,kms+1   ) = 0.0

            ! b. Define lower diagonal terms
            aaT(1,kms+1   ) = 0.0
            aaT(1,k1s+1:kms) = - dcoeff(k1s :kms-1) /                           &
                &  (hp(k1s:kms-1,m)  * 0.50 * (hp(k1s+1:kms,m)+hp(k1s:kms-1,m)))
            aaT(1,k1s     ) = 0.0

            ! c. Define center diagonal terms
            aaT(2,k1s     )  = 3./twodt1
            aaT(2,k1s+1:kms) = 3./twodt1 - aaT(1,k1s+1:kms) - aaT(3,k1s+1:kms)  &
                &          + 2. * ( DSQRT(q2p(k1s+1:kms,m))/ (B_1*xl(k1s+1:kms)) )
            aaT(2,kms+1   )  = 3./twodt1

            !.....Form r.h.s. matrix [ds_t]................................
            dsT(k1s+1:kms)  = 4. * q2p (k1s+1:kms,m)/twodt1                     &
                -      q2pp(k1s+1:kms,m)/twodt1
            DO k = k1s+1, kms
                IF ( BuoyP(k) >= 0.0 ) THEN
                    dsT(k) = dsT(k) + 2. * (ShearP(k) + BuoyP(k))
                ELSE
                    dsT(k) = dsT(k) + 2. * (ShearP(k)           )
                    aaT(2,k) = aaT(2,k)               - BuoyP(k) * 2. / q2p(k,m)
                ENDIF
            ENDDO

            ! ... Bottom boundary condition ...............................
            ustar = cd  * ( (0.5*(u(kms,m) + u(kms,lWC(m))))**2                 &
                &            +  (0.5*(v(kms,m) + v(kms,lSC(m))))**2 )
            dsT(kms+1) = 3.*B_1**(2./3.)*ustar/twodt1
 
            ! ... Surface boundary condition ..............................
            usurf = 0.5*(up(k1u(i,j),m) + up(k1u(i-1,j),lWC(m)) )
            vsurf = 0.5*(vp(k1v(i,j),m) + vp(k1v(i,j-1),lSC(m)) )
            ustar = cdw(i,j)*rhoair * ((uair(i,j)-usurf)**2.+     &
                (vair(i,j)-vsurf)**2.)/1000.
            dsT(k1s  ) = 3.*B_1**(2./3.)*ustar/twodt1
       
            !.....Solve tridiagonal system for vertical distribution of TKE
            CALL tridT (aaT, dsT, sal1T, k1s, kms+1, kms+2, nwlayers+1)

            !.....Define TKE at next time step....
            q2(k1s:kms+1,m) = MAX(q2_min,sal1T(1:nwlayers+1))
            q2(k1 :k1s-1,m) = q2(k1s,m)

            ! ... Define variables for basin scale energy balances
            IF (iobal > 0 .AND. istep > 1) THEN
                DO k = k1s+1, kms
                    rhoxh = (rhop(k,m)+1000.)*hp(k,m)
                    Dissipijk  = q2(k,m) * DSQRT(q2p(k,m))/ (B_1*xl(k))
                    ShearProduction    = ShearProduction    + ShearP(k)*rhoxh
                    BuoyancyProduction = BuoyancyProduction + BuoyP (k)*rhoxh
                    Dissipation        = Dissipation        + Dissipijk*rhoxh
                    TKinE              = TKinE              + q2(k,m)  *rhoxh
                    IF (k == kms) THEN
                        epskt(i,j) = Dissipijk ! m2/s3
                    !              PRINT *, 'Column number and Cell number:', m, k
                    !              PRINT *, 'Index i, j:', i, j
                    !              PRINT *, 'q2, xl:', q2(k,m), xl(k)
                    !              PRINT *, 'Dissipation rate:', epskt(i,j)
                    ENDIF
                ENDDO
            ENDIF

            ! ----- EQUATIONS for TURBULENT LENGTH SCALE ( q2l ) ----------


            ! ... Form tridiagonal matrix [aa] ............................

            ! a. Define upper diagonal terms (same as for q2)

            ! b. Define lower diagonal terms (same as for q2)

            ! c. Define center diagonal terms - Wall function as proposed by M&Y
            DO k = k1s+1, kms
                aaT(2,k) = 3./twodt1 - aaT(1,k) - aaT(3,k)                      &
                    + DSQRT(q2p(k,m))/ (B_1*xl(k))  *                      &
                    &   ( 1. + E_2 * ( xl(k)/(kappaS*zfromb(k)))** 2. +              &
                    &          E_3 * ( xl(k)/(kappaS*zfromt(k)))** 2. )
            ENDDO
            aaT(2,k1s  ) = 3./twodt1
            aaT(2,kms+1) = 3./twodt1

            !.....Form r.h.s. matrix [ds_t].....................................
            SELECT CASE (iturb)
                CASE (1) ! Original MY2.5 formulation
                    DO k = k1s+1, kms
                        IF ( BuoyP(k) >= 0 ) THEN
                            dsT(k  ) = xl(k   )*(E_1*ShearP(k)+E_1*BuoyP(k))
                        ELSE
                            dsT(k  ) = xl(k   )*(E_1*ShearP(k)             )
                            aaT(2,k) = aaT(2,k)              - E_1*BuoyP(k)/q2p(k,m)
                        ENDIF
                    ENDDO
                CASE (2:) !k-kl implementation of Burchard et. al 1999
                    DO k = k1s+1, kms
                        IF ( BuoyP(k) >= 0 ) THEN
                            dsT(k  ) = xl(k   )*(E_1*ShearP(k)+2.00*BuoyP(k))
                        ELSE
                            dsT(k  ) = xl(k   )*(E_1*ShearP(k)              )
                            aaT(2,k) = aaT(2,k)              - 5.06*BuoyP(k)/q2p(k,m)
                        ENDIF
                    ENDDO
            END SELECT

            dsT(k1s+1:kms) = dsT(k1s+1:kms)                       &
                &             + 4 * q2lp (k1s+1:kms,m)/twodt1        &
                &             -     q2lpp(k1s+1:kms,m)/twodt1

            ! ... Enforce Dirichlet bottom boundary condition ................
            dsT(kms+1) = 0.0
 
            ! ... Enforce Dirichlet surface boundary condition ...............
            dsT(k1s  ) = 0.0
         
            !.....Solve tridiagonal system for vertical distribution of TKE .
            CALL tridT (aaT, dsT, sal1T, k1s, kms+1, kms+2, nwlayers+1)

            !.....Define TKE at next time step....
            !     Make sure that Lc > alpha * Lk -Ivey&Imberger, 1991)
            q2l(k1s:kms+1,m) = MAX(q2l_min,sal1T(1:nwlayers+1))
            q2l(k1 :k1s-1,m) = q2l(k1s,m)

        END DO

        IF (iobal > 0 .AND. istep == 1) THEN
            ShearProduction    = ShearProduction    * dx * dy * dt
            BuoyancyProduction = BuoyancyProduction * dx * dy * dt
            Dissipation        = Dissipation        * dx * dy * dt
            TKinE              = 0.5 * TKinE        * dx * dy * dt
            ShearCum = ShearCum + ShearProduction
            BuocyCum = BuocyCum + BuoyancyProduction
            DisspCum = DisspCum + Dissipation
        ENDIF

    END SUBROUTINE SolveTKEandLengthScale

    !***********************************************************************
    SUBROUTINE tridT ( acoef, ds, sal, k1, km, km1, n )
        !***********************************************************************
        !
        !  Purpose: Tridiagonal matrix solver for scalar equation using
        !           the double-sweep method
        !
        !  Revisions:
        !    Date            Programmer        Description of revision
        !    ----            ----------        -----------------------
        !   3/7/00           F.F. Rueda        Original f90 code
        !
        !-----------------------------------------------------------------------

        !.....Dimensioning parameter.....
        INTEGER, PARAMETER :: kmax = 500

        !.....Arguments.....
        INTEGER, INTENT(IN) :: k1, km, km1, n
        REAL(real_G1), DIMENSION(km), INTENT(INOUT)  :: ds
        REAL(real_G1), DIMENSION(n),  INTENT(INOUT)  :: sal
        REAL(real_G1), DIMENSION(3,km1), INTENT(IN)  :: acoef

        !.....Local variables.....
        INTEGER :: k, kk
        REAL(real_G1), DIMENSION(kmax) :: a, b, c, d, e, f

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

          !.....Compute scalar in bottom layer.....
1       sal(n) = (d(n)-a(n)*f(n-1))/(b(n)+a(n)*e(n-1))

        !.....Backward sweep (solve for scalar vector).....
        DO k = n-1, 1, -1
            sal(k) = e(k) * sal(k+1) + f(k)
        END DO

    END SUBROUTINE tridT

    !***********************************************************************
    SUBROUTINE InitializeTurbulenceModel
        !***********************************************************************
        !
        !  Purpose: To define initial conditions for turbulent quantities.
        !
        !  Replace l (index) by m so that there ins not any overlap with variable
        !  L (length scale) used in gotm
        !
        !-----------------------------------------------------------------------

        ! ... Local variables
        INTEGER :: i,j,k,m, ll
        INTEGER :: kms, k1s
        REAL    :: zfromb0, zfromt0, htot
        REAL, DIMENSION (km1) :: mixl, zfromb, zfromt
        INTEGER, PARAMETER :: namlst = 10

        SELECT CASE (iturb)

            ! ... GOTM
            CASE(-1)

                CALL init_turbulence(namlst,'gotmturb.nml',km)
                CALL init_tridiagonal(km)
                ! Note - if we do not use KPP methods we do not have to initializeq
                ! equation of state or kpp -
                DO m = 1, lm
                    si3dtke(k1:km1,m) = tke
                    si3deps(k1:km1,m) = eps
                    si3dlen(k1:km1,m) = L
                    Av     (k1:km1,m) = num
                    Dv     (k1:km1,m) = nuh
                ENDDO

            ! ... Fixed eddy viscosity & diffusivity
            CASE (0)

                !.....Loop horizontally over pressure-points.....
                DO ll = 1, lm

                    ! .... Map 3D-(i,j) from 2L-l index
                    i = l2i(ll); j = l2j(ll);

                    ! ... Define top & bottom layers
                    kms = kmz(i,j); k1s = k1z(i,j)

                    ! Set eddy coefficients at the free surface to zero
                    Av(k1s,ll) = 0.0; Dv(k1s,ll) = 0.0   ! unnecessary
                    IF ( kms > k1s ) THEN   ! Two or more layers
                        DO k = k1s+1, kms
                            Av(k,ll) = Av0;
                            Dv(k,ll) = Dv0;
                        END DO
                    END IF
                    ! Set eddy coefficients at the bottom boundary to zero
                    Av(kms+1,ll) = 0.0; Dv(kms+1,ll) = 0.0
                END DO
     
            ! ... Original si3D implementation
            CASE (1)

                ! ... Fundamental parameters & combinations for the quasi-equilibrium
                !     stability functions of Kantha & Clayson (1994).
                t_1 = A_2 * (1.-6.*A_1/B_1)
                t_2 = 3.*A_2*(6.*A_1+B_2*(1.- C_3))
                t_3 = B_1**(-1./3.)
                t_4 = 9.*A_1*(2.*A_1+A_2*(1.-C_2))
                t_5 = 9.*A_1*A_2

                q2   = q2_min;
                q2p  = q2_min;
                q2pp = q2_min;

                DO m = 1, lm
                    ! ... 3D-(i,j) indexes for l
                    i = l2i(m);
                    j = l2j(m);
                    ! Calculate mixing length
                    kms = kmz(i,j);
                    k1s = k1z(i,j);
                    htot = SUM(hp(:,m))
                    IF ( kms > k1s ) THEN
                        ! Compute the array of vertical distances from the bottom
                        ! of the water column to the top of each layer
                        zfromb0 = 0.0
                        DO k = kms, k1s+1, -1
                            zfromb0   = zfromb0 + hp(k,m)
                            zfromb(k) = zfromb0
                        END DO
                        ! Compute the array of vertical distances from the top
                        ! of the water column to the top of each layer
                        zfromt0 = 0.0
                        DO k = k1s+1, kms
                            zfromt0   = zfromt0 + hp(k-1,m)
                            zfromt(k) = zfromt0
                        END DO
                        DO k = k1s+1, kms
                            mixl(k) = kappaS*zfromt(k)*SQRT(MAX(1.-(zfromt(k)/htot),0.0E0))
                            q2l(k,m) = q2(k,m)*mixl(k)
                        END DO
                        q2l(k1s,m) = 0.0; q2l(kms+1,m) = 0.0;
                    END IF
                END DO
                !q2l   = q2l_min;
                q2lp  = q2l;
                q2lpp = q2l;

        END SELECT

    END SUBROUTINE InitializeTurbulenceModel

    !***********************************************************************
    SUBROUTINE settrap_2EqTVars
        !***********************************************************************
        !
        !  Purpose: To setup the arrays at the  n  and  n+1/2  time levels
        !           for use in the first iteration of the  trapezoidal step.
        !           q2 and q2l at next time level are only solved once
        !           in a iteration cycle.
        !
        !  Replace l (index) by m so that there ins not any overlap with variable
        !  L (length scale) used in gotm
        !
        !-----------------------------------------------------------------------

        INTEGER :: i,j,k,m,k1s,kms

        DO m = 1,lm

            ! ... 3D-(i,j) indexes for l (= m)
            i = l2i(m); j = l2j(m);

            ! ... Top & bottom cells
            kms = kmz(i,j); k1s = k1z(i,j);

            DO k = k1, kms

                q2pp (k,m) = q2p (k,m)
                q2lpp(k,m) = q2lp(k,m)
                q2p  (k,m) = 0.5*(q2 (k,m)+q2pp (k,m))
                q2lp (k,m) = 0.5*(q2l(k,m)+q2lpp(k,m))

            ENDDO

            ! ... Define values at subsurface if wetting occurs
            k = k1s;
            IF ( hpp(k,m) <= ZERO ) THEN
                q2lp(k+1,m)=q2lp(k+2,m);
                q2p (k+1,m)=q2p (k+2,m);
            ENDIF

        ENDDO


    END SUBROUTINE settrap_2EqTVars


    !***********************************************************************
    SUBROUTINE settrap2_2EqTVars
        !***********************************************************************
        !
        !  Purpose: To setup the arrays at the  n  and  n+1/2  time levels
        !           for use in the first iteration of the  trapezoidal step.
        !           q2 and q2l at next time level are only solved once
        !           in a iteration cycle.
        !
        !  Replace l (index) by m so that there ins not any overlap with variable
        !  L (length scale) used in gotm
        !
        !-----------------------------------------------------------------------

        INTEGER :: i,j,k,m,k1s,kms

        DO m = 1,lm

            ! ... 3D-(i,j) indexes for l (= m)
            i = l2i(m); j = l2j(m);

            ! ... Top & bottom cells
            kms = kmz(i,j); k1s = k1z(i,j);

            DO k = k1, kms
                q2p (k,m) = 0.5*(q2 (k,m)+q2pp (k,m))
                q2lp(k,m) = 0.5*(q2l(k,m)+q2lpp(k,m))
            ENDDO

            ! ... Redefine q2&q2l at subsurface cell if wetting occurs
            !     from n+1/2(n) to n+1
            k = k1s;
            IF ( hpp(k,m) <= ZERO ) THEN
                q2lp(k+1,m)=q2lp(k+2,m);
                q2p (k+1,m)=q2p (k+2,m);
            ENDIF


        ENDDO

    END SUBROUTINE settrap2_2EqTVars

    !***********************************************************************
    SUBROUTINE save_2EqTvars
        !***********************************************************************
        !
        !  Purpose: To save solution of variables implicated in 2-Equation
        !           turbulence models
        !
        !  Replace l (index) by m so that there ins not any overlap with variable
        !  L (length scale) used in gotm
        !
        !-----------------------------------------------------------------------
   
        INTEGER :: i,j,k,m,k1s,kms

        SELECT CASE (istep)
            CASE (1)
    
                q2pp  = q2p ;
                q2p   = q2  ;
                q2lpp = q2lp;
                q2lp  = q2l ;

                DO m = 1,lm
                    ! ... 3D-(i,j) indexes for l (= m)
                    i = l2i(m); j = l2j(m);
                    ! ... Top & bottom cells
                    k1s = k1z(i,j);
                    ! ... Redo calculations if wetting occurs in
                    IF ( hpp(k1s,m) <= ZERO ) THEN
                        q2lp(k1s+1,m)=q2lp(k1s+2,m);
                        q2p (k1s+1,m)=q2p (k1s+2,m);
                    ENDIF
                ENDDO

            CASE (2)

                q2p   = q2  ;
                q2lp  = q2l ;
 
                DO m = 1,lm
                    ! ... 3D-(i,j) indexes for l (= m)
                    i = l2i(m); j = l2j(m);
                    ! ... Top & bottom cells
                    k1s = k1z(i,j);
                    IF ( km1 - k1s > 1 ) THEN
                        ! ... Redo calculations if wetting occurs in
                        IF ( hpp(k1s,m) <= ZERO ) THEN
                            q2lp(k1s+1,m)=q2lp(k1s+2,m);
                            q2p (k1s+1,m)=q2p (k1s+2,m);
                        ENDIF
                    ENDIF
                ENDDO


        END SELECT
 
    END SUBROUTINE save_2EqTvars



END MODULE si3d_mixing
