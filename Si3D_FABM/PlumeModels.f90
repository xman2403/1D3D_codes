Subroutine rk4(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, y, dydx, nn, x, h, yout)

Integer i, nn, nmax
Parameter (nmax=50)
Real *8 e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, h, x, dydx(nn), y(nn), yout(nn), h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)
External derivs
hh = h*0.5
h6 = h/6.
xh = x + hh
Do i = 1, nn
yt(i) = y(i) + hh*dydx(i)
End Do
Call derivs(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, xh, yt, dyt)
Do i = 1, nn
yt(i) = y(i) + hh*dyt(i)
End Do
Call derivs(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, xh, yt, dym)
Do i = 1, nn
yt(i) = y(i) + h*dym(i)
dym(i) = dyt(i) + dym(i)
End Do
Call derivs(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, x+h, yt, dyt)
Do i = 1, nn
yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.*dym(i))
End Do
Return
End Subroutine rk4
!
!----------------------------------------------------------------------
!
Subroutine derivs(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, x, y, dydx)
Real *8 e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, x, y(8), dydx(8)
!     Right-hand side of differential equations for Runge-Kutta solution
dydx(1) = e
dydx(2) = (densea-densew)/densep*g*l*2.*b + (densew-densep)/densep*g*lambda*2.*b*(l-2.*b*(1.-lambda))
dydx(3) = e*tamb
dydx(4) = e*(salamb*gamma/dense20)*densea
dydx(5) = (e*doamb/32.+4.0*pi*rb**2*n/(v+vb)*kolo*(ho2*po-comgp/32.))
dydx(6) = (e*dnamb/28.+4.0*pi*rb**2*n/(v+vb)*koln*(hn2*pn-cnmgp/28.))
dydx(7) = -4.0*pi*rb**2*n/(v+vb)*kolo*(ho2*po-comgp/32.)
dydx(8) = -4.0*pi*rb**2*n/(v+vb)*koln*(hn2*pn-cnmgp/28.)
!print *, 'dydx' , dydx(1:8)
Return
End Subroutine derivs
!
!----------------------------------------------------------------------
!
Subroutine linint(xtab, ytab, ntab, x, y)
Integer i, ntab
Real *8 x, y, xtab(ntab), ytab(ntab)
If (x<xtab(1) .Or. x>xtab(ntab)) Then
Write (*, *) 'X = ', x, ' IS OUT OF TABLE RANGE'
Pause
End If
Do i = 2, ntab
If (x<=xtab(i)) Goto 200
End Do
200 i1 = i - 1
wx = (x-xtab(i1))/(xtab(i1+1)-xtab(i1))
y = (1.-wx)*ytab(i1) + wx*ytab(i1+1)
Return
End Subroutine linint


!     *******************************************************************************
Subroutine doubplu_v1(year, julday, wsel, diffel, layers, ldiff, lambnot, salamb, patm, diamm, lake, qscfm, frconot, laydiff, hcell, locat, te, co2m, elevob, qwob, tplumb, comgpb, qwd, bwd, laytop, layob, oteff)
!     *******************************************************************************
!     THIS SUBROUTINE IS WRITTEN TO PREDICT THE PERFORMANCE OF A DOUBLE BUBBLE PLUME.
!
!     VERSION 1
!     THE MODEL IS BASED ON SOCOLOSKY ET AL. (2008) DOUBLE PLUME MODEL.
!     This version includes the following:
!
!     June , 2018
!     Revised By: Xiamei Man
!     VARIABLES
!
!     ALPHAI=ENTRAINMENT COEFFICIENT OF INNER PLUME
!     ALPHAO=ENTRAINMENT COEFFICIENT OF OUTER PLUME
!     ALPHAA=ENTRAINMENT COEFFICIENT OF AMBIENT WATER
!     B=1/2 DIFFUSER WIDTH (m)
!     BH(ID)=HALF INNER PLUME WIDTH AT DEPTH OF ID (m)
!     BO=HALF OUTER PLUME WIDTH (m)
!     BOH(ID)=HALP OUTER PLUME WIDTH AT DEPTH OF ID(m)
!     BAVG=AVERAGE 1/2 DIFFUSER WIDTH (m)
!     CO2=DISSOLVED OXYGEN (DO) CONCENTRATION (mol/m3)
!     CO2H(ID)=DISSOLVED OXYGEN (DO) CONCENTRATION WHERE DEPTH = ID(mol/m3)
!     CO2M=DO CONCENTRATION PROFILE FOR INPUT BOUNDARY CONDITION (g/m3)
!     CN2=DISSOLVED NITROGEN CONCENTRATION (mol/m3)
!     CN2H(ID)=DISSOLVED NITROGEN CONCENTRATION WHERE DEPTH = ID(mol/m3)
!     COMG=DISSOLVED OXYGEN CONCENTRATION (g/m3)
!     COMGP=DISSOLVED OXYGEN CONCENTRATION OF INNER PLUME (g/m3)
!     COMGPH(ID)=DISSOLVED OXYGEN CONCENTRATION OF INNER PLUME WHERE DEPTH = ID (g
!     /m3)
!     COMGPT=DO CONCENTRATION OF PLUME DETRAINMENT AT TOP OF PLUME (g/m3)
!     COMGPB=DO CONCENTRATION OF PLUME DETRAINMENT AT BOTTOM OF OUTER PLUME (g/m3)
!     CNMG=DISSOLVED NITROGEN CONCENTRATION (g/m3)
!     CNMGP=DISSOVED NITROGEN CONCENTRATION OF INNER PLUME (g/m3)
!     CNMGPH(ID)=DISSOLVED NITROGEN CONCENTRAION OF INNER PLUME WHERE DEPTH = ID (g/m3)
!     DENSE20=DENSITY OF WATER AT 20 C
!     DENSEA=AMBIENT WATER DENSITY (kg/m3)
!     DENSAH(ID)=AMBIENT WATER DENSITY WHERE DEPTH = ID (kg/m3)
!     DENSEP=DENSITY OF THE PLUME (kg/m3)
!     DENSEW=WATER DENSITY IN PLUME (kg/m3)
!     DENSEWH(ID)=WATER DENSITY IN PLUME WHERE DEPTH = ID (kg/m3)
!     DENSEWO=WATER DENSITY IN OUTER PLUME (kg/m3)
!     DENSEWOH(ID)=WATER DENSITY IN OUTER PLUME WHERE DEPTH = ID (kg/m3)
!     DIAMM=BUBBLE DIAMETER (mm)
!     DIFFEL=DIFFUSER ELEVATION (m)
!     DMPR=DEPTH OF MAXIMUM PLUME RISE (m)
!     DNAMB=AMBIENT DISSOLVED NITROGEN CONCENTRATION (g/m3)
!     DOAMB=AMBIENT DISSOLVED OXYGEN CONCENTRATION (g/m3)
!     E=ENTRAINMENT FACTOR OF INNER PLUME (m3/s)
!     EH(ID)=ENTRAINMENT FACTOR OF INNER PLUME WHERE DEPTH = ID (m3/s)
!     EO=ENTRAINMENT FACTOR OF OUTER PLUME(m3/s)
!     EOH(ID)=ENTRAINMENT FACTOR OF OUTER PLUME WHERE DEPTH = ID (m3/s)
!     EA=ENTRAINMENT FACTOR OF AMBIENT WATER(m3/s)
!     ELEV=ELEVATION (m)
!     ELEVT=TERMINAL ELEVATION OF PLUME IN SEGMENT (m)
!     ELEVOB=TERMINAL ELEVATION OF OUTER PLUME IN SEGMENT (m)
!     FDO=DISSOLVED OXYGEN FLUX (mol/s)
!     FDOO=DISSOLVED OXYGEN FLUX IN OUTER PLUME (mol/s)
!     FDN=DISSOLVED NITROGEN FLUX (mol/s)
!     FDNO=DISSOLVED NITROGEN FLUX IN OUTER PLUME (mol/s)
!     FRACO=MOLE FRACTION OF OXYGEN (-)
!     FRACN=MOLE FRACTION OF NITROGEN (-)
!     FRCONOT=INITIAL MOLE FRACTION OF OXYGEN IN DIFFUSER GAS SUPPLY, 0.21 OR 0.965 (-)
!     FSAL=SALINITY FLUX (kg/s)
!     FSALO=SALINITY FLUX OF OUTER PLUME (kg/s)
!     FTEMP=TEMPERATURE FLUX (m3/s)
!     FTEMPO=TEMPERATURE FLUX IN OUTER PLUME (m3/s)
!     FGO=GASEOUS OXYGEN FLUX (mol/s)
!     FGONOT=INITIAL GASEOUS OXYGEN FLUX (mol/s)
!     FGOO=GASEOUS OXYGEN FLUX IN OUTER PLUME (mol/s)
!     FGN=GASEOUS NITROGEN FLUX (mol/s)
!     FGNO=GASEOUS NITROGEN FLUX IN OUTER PLUME (mol/s)
!     FRCNOT=FRACTION OF NITROGEN IN ATMOSPHERE (-)
!     FRNOT=INITIAL FROUDE NUMBER (-)
!     FRNOTO=INITIAL FROUDE NUMBER FOR OUTER PLUME
!     GAMMA=SALINITY CONVERSION FACTOR [kg/m3/(uS/cm)]
!     GROSSMT=GROSS MASS TRANSFER OF OXYGEN FROM PLUME (kg/d)
!     HCELL=HEIGHT OF CELL IN GRID (m)
!     HO2=SOLUBILITY CONSTANT FOR OXYGEN (mol/m3/Pa)
!     HO2H(ID)=SOLUBILITY CONSTANT FOR OXYGEN WHERE DEPTH = ID (mol/m3/Pa)
!     HN2=SOLUBILITY CONSTANT FOR NITROGEN (mol/m3/Pa)
!     HN2H(ID)=SOLUBILITY CONSTANT FOR NITROGEN WHERE DEPTH = ID (mol/m3/Pa)
!     HWITH=HEIGHT OF WITHDRAWAL/ENTRAINMENT ZONE (m)
!     JULDAY=JULIAN DAY IN GIVEN YEAR
!     KOLO=MASS TRANSFER COEFFICIENT FOR OXYGEN (m/s)
!     KOLOH(ID)=MASS TRANSFER COEFFICIENT FOR OXYGEN WHERE DEPTH = ID (m/s)
!     KOLN=MASS TRANSFER COEFFICIENT FOR NITROGEN (m/s)
!     KOLNH(ID)=MASS TRANSFER COEFFICIENT FOR NITROGEN WHERE DEPTH = ID (m/s)
!     L=DIFFUSER LENGTH (m)
!     LAKE=LAKE AND DIFFUSER TYPE (SHR AND LINEAR=1 OR AMISK AND RECTANGULAR=2) FOR SELECTION OF LAMBDA
!     LAYERS=NUMBER OF LAYERS/DATA POINTS IN BOUNDARY CONDITION PROFILES (-)
!     LAYDIFF=GRID LAYER CORRESPONDING TO DIFFUSER DEPTH (-)
!     LAYTOP=GRID LAYER CORRESPONDING TO TOP OF PLUME (-)
!     LAYOB=GRID LAYER CORRESPONDING TO BOTTOM OF OUTER PLUME
!     LDIFF=LENGTH OF DIFFUSER (m)
!     LNOT=INITIAL DIFFUSER LENGTH (m)
!     LAMBDA=FRACTION OF PLUME OCCUPIED BY BUBBLES (-)
!     LAMBNOT=LAMBDA x INITIAL PLUME RADIUS;EQUAL TO DIFFUSER RADIUS (m)
!     LOCAT=DEPTHS FOR INPUT BOUNDARY CONDITION PROFILES (m)
!     MOMENT=MOMENTUM (m4/s)
!     MOMENTO=MOMENTUM OF OUTER PLUME (m4/s)
!     N=NUMBER OF BUBBLES PER SECOND (1/s)
!     OTEFF=OXYGEN TRANSFER EFFICEINCY (%)
!     PATM=ATMOSPHERIC PRESSURE AT AVERAGE WSEL (Pa)
!     PSTD=STANDARD PRESSURE (Pa)
!     PO=IN SITU PARTIAL PRESSURE OF GAS PHASE OF OXYGEN(Pa)
!     POH(ID)=IN SITU PARTIAL PRESSURE OF GAS PHASE OF OXYGEN WHERE DEPTH = ID (Pa)
!     PN=IN SITU PARTIAL PRESSURE OF GAS PHASE OF NITROGEN (Pa)
!     PNH(ID)=IN SITU PARTIAL PRESSURE OF GAS PHASE OF NITROGEN WHERE DEPTH = ID (Pa)
!     QSCFM=STANDARD GAS FLOW RATE (scfm), TOTAL GAS FLOW RATE TO DIFFUSER
!     QSCMS=STANDARD GAS FLOW RATE (scms)
!     QW=FLOWRATE OF WATER (m3/s)
!     QWO=FLOWRATE OF WATER IN OUTER PLUME (m3/s)
!     QWT=TOTAL DETRAINMENT FLOW RATE AT THE TOP OF THE PLUME (m3/s)
!     QWOB=TOTAL DETRAINMENT FLOW RATE AT THE BOTTOM OF THE OUTER PLUME (m3/s)
!     RB=BUBBLE RADIUS (m)
!     RBH(ID)=BUBBLE RADIUS WHERE DEPTH = ID (m)
!     RGAS=IDEAL GAS CONSTANT (J/mol/K)
!     SAL=SALINITY (uS/cm)
!     SALAMB=SALINITY OF AMBIENT WATER (uS/cm)
!     SALPL=SALINITY OF THE PLUME (uS/cm)
!     SALPLH(ID)=SALINITY OF THE PLUME WHERE DEPTH = ID (uS/cm)
!     SALPLO=SALINITY OF THE OUTER PLUME (uS/cm)
!     SAPLOH(ID)=SALINITY OF THE OUTER PLUME WHERE DEPTH = ID (uS/cm)
!     TAMB=AMBIENT WATER TEMPERATURE (C)
!     TAMBH(ID)=AMBIENT WATER TEMPERATURE WHERE DEPTH = ID (C)
!     TAVG=AVERAGE AMBIENT WATER TEMPERATURE (C)
!     TDS=TOTAL DISSOLVED SOLIDS (g/m3) [0.64 conversion factor from Chapra book]
!     TE=TEMPERATURE PROFILE FOR INPUT BOUNDARY CONDITION (C)
!     TPLUME=PLUME TEMPERATURE (C)
!     TPLUMEI=PLUME TEMPERATURE OF INNER PLUME (C)
!     TPLUMEH(ID)=PLUME TEMPERATURE WHERE DEPTH = ID (C)
!     TPLUMET=DETRAINMENT PLUME TEMPERATURE AT THE TOP OF THE PLUME (C)
!     TPLUMEOB=DETRAINMENT PLUME TEMPERATURE AT THE BOTTOM OF THE OUTER PLUME (C)
!     TPLUMEO=PLUME TEMPERATURE OF OUTER PLUME (C)
!     TPLUMEOH(ID)=PLUME TEMPERATURE OF OUTER PLUME WHERE DEPTH = ID (C)
!     TSTD=STANDARD TEMPERATURE (K)
!     V=WATER VELOCITY (m/s)
!     VH(ID)=WATER VELOCITY WHERE DEPTH = ID (m/s)
!     VI=WATER VELOCITY OF INNER PLUME (m/s)
!     VO=WATER VELOCITY OF OUTER PLUME (m/s)
!     VOH(ID)=WATER VELOCITY OF OUTER PLUME WHERE DEPTH = ID (m/s)
!     VAVG=AVERAGE WATER VELOCITY (m/s)
!     VB=BUBBLE RISE VELOCITY (m/s)
!     VBH(ID)=BUBBLE RISE VELOCITY AT DEPTH WHERE DEPTH = ID (m/s)
!     VBUB=BUBBLE VOLUME (m3)
!     VGUESS=GUESSED INITIAL WATER VELOCITY (m/s)
!     VGUESO=GUESSED INITIAL WATER VELOCITY OF OUTER PLUME (m/s)
!     WSEL=WATER SURFACE ELEVATION (m)
!     ID=NUMBER OF TIMES FOR RK4 SOLUTION PROCEDURE IN ONE ITERATION
!     IDM=MAMXIMUM NUMBER OF ITERATION FOR INNER PLUME(IN THE FIRST ITERATION)
!     YO2=GASEOUS OXYGEN CONCENTRATION (mol/m3)
!     YN2=GASEOUS NITROGEN CONCENTRATION (mol/m3)
!     Z=DEPTH TO DIFFUSER (m)
!     VAMB = Ambient Velocity (m) used to calculate VORTEX entrainment (FJR
!     PWD  = Perimeter (m) (FJR)
!     VUP  = Upward velocity (m/2)
!     AWD  = Area (m2) of plume
!     BWD  = Width
!
Real *8 alpha, alphaa, alphai, alphao, area, b, bo, co2, comg, cn2, cnmg, ds, dz, & 
densea, densep,denseo, densew, diamm, e, eo, eoh, fdo, fdoo, fdn, &
fdno,  fraco, fracn, fsal, fsalo, ftemp, ftempo, fgo, fgoo, fgn, fgno, g, gamma, ho2, hn2, &
 koln, kolo, l, lambda, moment, momeno, n, pi, po, pn, &
pstd, pz, qscfm, qscms, qw, qwo, qgas, rb, rgas, salamb, salpl, salplo, &
saploh, tamb, tplume, tplumo, tplumoh, tstd, v, vi, vo, voh, vb, &
vbub, vg, vguess, vgueso, yo2, yn2, z, aa, bb, cc, bnot, lnot, elev, dt
!Real , Allocatable , Dimension (:) :: bh, co2h , cn2h, densah , denswh , doambh , eh , ho2h , hn2h , kolnh, &
!koloh,poh, pnh, rbh , salplh , tambh , tplumeh , vh , vbh , comgph , cnmgph 
Real *8 bh(100),co2h(100),cn2h(100),densah(100),denswh(100),doambh(100),eh(100),ho2h(100),hn2h(100), &
kolnh(100),koloh(100),poh(100),pnh(100),rbh(100),salplh(100),tambh(100),tplumh(100),vh(100),vbh(100), &
comgph(100),cnmgph(100) 
Real *8 te(50), xloc, depth, co2m(50), locat(50), doamb, comgp, cnmgp, wsel, patm,&
 sal(50), diffel, grossmt, fgonot, dnamb, dmpr, tavg, sumtemp, sumsal, lambnot, frconot, rbnot, h, ho, &
 dydx(8), y(8), yout(8), dense20, oteff, frnot, frnoto, vdiff, fr, buoy, dco2, qgfrac, ldiff, tds(100), &
jday, el(70), deltac, comgnot, elevt, elevob, qwt, qwob, tplumt, tplumb, comgpt, comgpb, qwd(50), &
pwd(50), bwd(50), awd(50), hwith, hcell, julday, btop, bob, bequiv
Integer ii, ij, ik, in, jj, ll, neqn, nn, mi, jk, jl, lake, laytop, layob, layers, km, kn, ko, &
kp, rows, kq, kr, ku, kv, kw, kx, ky, kz, ks, laydiff, year
Integer nels, id, idm
!
!allocate(bh(100))
!allocate(co2h(100))
!allocate(cn2h(100))
!allocate(densah(100))
!allocate(denswh(100))
!allocate(doambh(100))
!allocate(eh(100))
!allocate(ho2h(100))
!allocate(hn2h(100))
!allocate(kolnh(100))
!allocate(koloh(100))
!allocate(poh(100))
!allocate(pnh(100))
!allocate(rbh(100))
!allocate(salplh(100))
!allocate(tambh(100))
!allocate(tplumeh(100))
!allocate(vh(100))
!allocate(vbh(100))
qgfrac = 1.0
sumtemp = 0.0
Do ks = 1, layers
sumtemp = sumtemp + te(ks)
End Do
tavg = sumtemp/layers
!
!     Assume that gas bubbles are composed of oxygen and nitrogen only.
fraco = frconot
fracn = 1.0 - fraco
depth = wsel - diffel
z = & !from free surface to diffusers, commented by Chris
depth
elev = diffel
x = 0.
!
!     Interpolate input profiles to obtain line plume initial conditions
xloc = depth - x
Call linint(locat, te, layers, xloc, tamb)
Call linint(locat, co2m, layers, xloc, comg)
comgp = comg
doamb = comg
co2 = comg/32.
comgnot = comg
salpl = salamb
tplume = tamb
vguess = 0.07
vgueso = -0.07
v = vguess
!
!     CONSTANTS
alpha = & ! Schladow 1993 original 0.11
0.083
g = 9.80665
gamma = 6.9E-4
If (lake==1) Then
lambda = 0.93
frnot = 2.0
Else If (lake==2) Then
lambda = 0.8
frnot = 1.6
frnoto = 0.1
End If
pi = acos(-1.0)
pstd = 101325.
rgas = 8.314
tstd = 293.15
dense20 = 998.2
frcnatm = 0.79
bnot = lambnot/lambda
b = bnot
!     Revised LNOT to account for additional length due to spreading of velocity/water plume beyond bubble plume (4-16-09)
lnot = ldiff + 2.0*bnot*(1.0-lambda)
l = lnot
bequiv = 0.5*(4.*ldiff*2.*lambnot/pi)**0.5
!
!     AMBIENT AND AVERAGE WATER DENSITIES
densea = (0.059385*tamb**3-8.56272*tamb**2+65.4891*tamb)*0.001 + 999.84298 + (gamma)*salamb
densew = densea
!
!     BUBBLE PROPERTIES
!     Gas flow rate per segment asssumed to be proportional to fraction of total diffuser length.
qscms = qgfrac*qscfm/3.281**3/60.0
qgas = pstd*qscms*(tamb+273.15)/((patm+densea*g*z)*tstd)
!     For diffuser in SHR, use correlation by McGinnis and Little (2000) for initial bubble size.
rb = diamm/2000.
rbnot = rb
If (rb<=(7.5E-4)) Then
vb = 1189.0*rb**1.1945
Else If (rb>(7.5E-4) .And. rb<(4.8E-3)) Then
vb = 0.22
Else
vb = 2.995*rb**0.489
End If
!
kolo = 0.6*rb
If (kolo>(4.0E-4)) Then
kolo = 4.0E-4
End If
koln = kolo
!
ho2 = (2.125-0.05023*tplume+5.7714E-4*tplume**2)/100000.
hn2 = (1.042-0.02457*tplume+3.1714E-4*tplume**2)/100000.
!
!     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
 cn2 = (patm*frcnatm)*hn2
 cnmg = cn2*28.0
 cnmgp = cnmg
dnamb = cnmg

!     CALCULATION OF INITIAL WATER VELOCITY USING FROUDE NUMBER
vbub = 4./3.*pi*rb**3
n = qgas/vbub
9 vg = qgas/((vguess+vb)*(2.*lambda*b)*(l-2.0*b*(1.0-lambda)))
densep = (1.0-vg)*densew
If (lake==1) Then
v = frnot*(2.0*lambda*b*g*(densea-densep)/densep)**0.5
Else If (lake==2) Then
v = frnot*(2.0*bequiv*g*(densea-densep)/densep)**0.5
!         VLS: For testing purposes, assume that characteristic length is equal to rectangle width (6-4-09)
!          V=FRNOT*(2.0*LAMBDA*B*G*(DENSEA-DENSEP)/DENSEP)**0.5
End If
vdiff = abs(v-vguess)
If (vdiff>1.0E-6) Then
vguess = v
Goto 9
End If
!
!     VARIABLE TRANSFORMATION
e = 2.*(l+2.*b)*alpha*v
qw = 2.*l*b*v
moment = 2.*l*b*v**2
ftemp = qw*tplume
fsal = qw*(salpl*gamma/dense20)*densew
!     Previous equation corrected to account for salinity units conversion.
fdo = qw*co2
fdn = qw*cn2
fgo = pstd*qscms/(rgas*tstd)*fraco
fgonot = fgo
fgn = pstd*qscms/(rgas*tstd)*fracn
!     Revised gaseous flux equations.
yo2 = fgo/(lambda*2.*b*(l-2.*b*(1.-lambda))*(v+vb))
yn2 = fgn/(lambda*2.*b*(l-2.*b*(1.-lambda))*(v+vb))
pz = patm + (densea*g*z)
po = pz*fraco
pn = pz*fracn
buoy = (g*(densea-densep)/densep*qw)/lnot
tds = salpl*0.64
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
jj = 0
qwd(laydiff) = qw
pwd(laydiff) = 2.*(l+2.*b)
awd(laydiff) = l*b
bwd(laydiff) = sqrt(awd(laydiff)/3.141592)
!
!    SOLUTION PROCEDURE
dz = 0.001
h = 0.001
ho = -0.001
nels = 0
id = 1
idm = 0
hwith = 0.0
10 z = z - dz
x = x + dz
elev = elev + dz
nels = nels + 1
idm = idm + 1
 if ((MOD(idm,500)) == 0) Then
 id = id + 1
 end if
!Print *, id
!
!     Interpolate input profiles to obtain line plume boundary conditions
xloc = depth - x
Call linint(locat, co2m, layers, xloc, comg)
doamb = comg
doambh(id) = doamb
!print *,doambh(id)
Call linint(locat, te, layers, xloc, tamb)
tambh(id) = tamb
!print *, tambh(id)
!     Use subroutines for Runge Kutta method solution
neqn = 8
y(1) = qw
y(2) = moment
y(3) = ftemp
y(4) = fsal
y(5) = fdo
y(6) = fdn
y(7) = fgo
y(8) = fgn
Call derivs(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, z, y, dydx)
Call rk4(e, densea, densew, densep, g, l, b, lambda, tamb, salamb, gamma, dense20, doamb, pi, rb, n, v, vb, kolo, ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, y, dydx, neqn, z, h, yout)
!
qw = yout(1)
moment = yout(2)
ftemp = yout(3)
fsal = yout(4)
fdo = yout(5)
fdn = yout(6)
fgo = yout(7)
fgn = yout(8)
If (moment<0.0) Then
tplume = ftemp/qw
tplumh(id) = tplume
!print *, tplumh(id)
salpl = fsal/(qw*densew)/(gamma/dense20)
salplh(id) = salpl
!print *, salplh(id)
!    Previous equation corrected to consistently express salinity in uS/cm
 co2 = fdo/qw
 co2h(id) = co2
!print *, co2h(id)
 cn2 = fdn/qw
 cn2h(id) = cn2
!print *, cn2h(id)
Goto 20
End If
v = moment/qw
vh(id) = v
!print *,vh(id)
area = qw/v
!     SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
 aa = 1.0
 bb = 2.*bnot - lnot
 cc = -1.0*area
l = (-1.0*bb+(bb**2-4.0*aa*cc)**(0.5))/(2.0*aa)
If (l<0.0) Then
l = (-1.0*bb-(bb**2-4.0*aa*cc)**(0.5))/(2.0*aa)
End If
b = area/(2.0*l)
bh(id) = b
!print *,bh(id)
e = 2.*(l+2.*b)*alpha*v
eh(id) = e
!print *,eh(id)
!     Add incremental entrainment to total cell entrainment/withdrawal
qwd(laydiff-jj) = qwd(laydiff-jj) + e*dz
pwd(laydiff-jj) = pwd(laydiff-jj) + 2.*(l+2.*b)
awd(laydiff-jj) = awd(laydiff-jj) + l*b*2.
bwd(laydiff-jj) = bwd(laydiff-jj) + sqrt(l*b*2./3.141592)
hwith = hwith + dz
If (hwith>hcell) Then
pwd(laydiff-jj) = pwd(laydiff-jj)/nels
awd(laydiff-jj) = awd(laydiff-jj)/nels
bwd(laydiff-jj) = bwd(laydiff-jj)/nels
Print *, qwd(laydiff-jj), bwd(laydiff-jj), v
jj = jj + 1
hwith = 0.0
qwd(laydiff-jj) = 0.0
pwd(laydiff-jj) = 0.0
bwd(laydiff-jj) = 0.0
awd(laydiff-jj) = 0.0
nels = 0
End If
tplume = ftemp/qw
tplumh(id) = tplume
!print *,tplumeh(id)
salpl = fsal/(qw*densew)/(gamma/dense20)
salplh(id) = salpl
!print *,salplh(id)
!     Previous equation corrected to consistently express salinity in uS/cm
 co2 = fdo/qw
 co2h(id) = co2
!print *,co2h(id)
 cn2 = fdn/qw
 cn2h(id) = cn2
!print *,cn2h(id)
 comgp = co2*32.
 comgph(id) = comgp
!print *,comgph(id)
 cnmgp = cn2*28.
 cnmgph(id) = cnmgp
!print*, cnmgph(id)
!     Revised gaseous flux equations.
yo2 = fgo/(lambda*2.*b*(l-2.*b*(1.-lambda))*(v+vb))
yn2 = fgn/(lambda*2.*b*(l-2.*b*(1.-lambda))*(v+vb))
!
pz = patm + (densea*g*z)
qgas = (fgo+fgn)*rgas*(tplume+273.15)/pz
vbub = qgas/n
!
vg = vbub*n/((v+vb)*(2.*lambda*b)*(l-2.*b*(1.-lambda)))
!     Previous equation revised to account for correct plume cross-sectional area occupied by bubbles.
rb = (3.*qgas/(4.*pi*n))**(1./3.)
If (rb<0.0) Then
rb = 1.0E-8
End If
rbh(id) = rb
!print *, rbh(id)
fraco = fgo/(fgo+fgn)
fracn = 1.0 - fraco
po = pz*fraco
poh(id) = po
!print *, poh(id)
pn = pz*fracn
pnh(id) = pn
!print *,pnh(id)
densea = (0.059385*tamb**3-8.56272*tamb**2+65.4891*tamb)*0.001 + 999.84298 + (gamma)*salamb
densah(id) = densea
!print *,densah(id)
densew = (0.059385*tplume**3-8.56272*tplume**2+65.4891*tplume)*0.001 + 999.84298 + (gamma)*salpl
denswh(id) = densew
!print *,denswh(id)
!     Previous equation re-revised to account for correct salinity units (uS/cm) in density calculations.
densep = (1.0-vg)*densew
!
!    BUBBLE PROPERTIES
If (rb<=(7.5E-4)) Then
vb = 1189.0*rb**1.1945
Else If (rb>(7.5E-4) .And. rb<(4.8E-3)) Then
vb = 0.22
Else
vb = 2.995*rb**0.489
End If
vbh(id) = vb
!print *,vbh(id)
kolo = 0.6*rb
If (kolo>(4.0E-4)) Then
kolo = 4.0E-4
End If
koloh(id) = kolo
!print *,koloh(id)
koln = kolo
kolnh(id) = koln
!print *, kolnh(id)
ho2 = (2.125-0.05023*tplume+5.7714E-4*tplume**2)/100000.
ho2h(id) = ho2
!print *, ho2h(id)
hn2 = (1.042-0.02457*tplume+3.1714E-4*tplume**2)/100000.
hn2h(id) = hn2
!print *, hn2h(id)
fr = v/(2.*lambda*b*g*(densea-densep)/densep)**0.5
dco2 = ho2*po - co2
If (v>1.E-6) Then
If (z>0.0) Then
Goto 10
End If
End If
!
!     CALCULATION OF AVERAGE NET OXYGEN MASS TRANSFER FOR DAY
20 grossmt = (fgonot-fgo)*32./1000.*86400.
oteff = (fgonot-fgo)/fgonot*100.
deltac = comgp - comgnot
!
Print *, 'O2 Trans Eff:', oteff
elevt = elev
qwt = qw
!print *, qwt
tplumt = tplume
!print *, tplumt
comgpt = & !*0.9 ! Added '*0.9' on 26 Nov to try to decrease O2
comgp
!print *,comgp
laytop = laydiff - jj
btop = b
!    Start outer plume calculation
salplo = salamb
tplumo = tamb
vi = v
vgueso = -0.05
vo = vgueso
!
!     CONSTANTS
alphai = 0.055
alphao = 0.11
alphaa = 0.11
!     Assume initial ambient dissolved nitrogen conc. equals saturated conc. at surface.
 cnmg = cn2*28.0
dnamb = cnmg
!     calculation of initial water velocity of outer plume using froude
vbub = 4./3.*pi*rb**3
n = qgas/vbub
bo = 0.4521
!11 bo = (qwt**2./(pi*qwt*vgueso+b**2.))**0.5
!print *,bo
!vo = frnoto*((bo-b)*g*abs(densea-densep)/densep)**0.5
!print *,vo
!vdiff = abs(vo-vgueso)
!If (vdiff>1.0E-6) Then
vo = vgueso
!Goto 11
!ENDIF
!
!     VARIABLE TRANSFORMATION
ei = 2.*pi*b*alphai*vi
eo = -1.*pi*b*alphao*vo
ea = -1.*pi*bo*alphaa*vo
qwo = qwt
!print *, qwo
salplo =salpl
!print *, 'salplo,salpl', salplo
momeno = -2.*l*bo*vo**2
denseo = densep
!print *, 'denseo', denseo
ftempo = qwo*tplumo
!print *,ftempo
fsalo = qwo*(salplo*gamma/dense20)*densew
!     Previous equation corrected to account for salinity units conversion.
fdoo = qwo*co2
fdno = qwo*cn2
fgoo = pstd*qscms/(rgas*tstd)*fraco
fgno = pstd*qscms/(rgas*tstd)*fracn
!     Revised gaseous flux equations.
pz = patm + (densea*g*z)
po = pz*fraco
pn = pz*fracn
buoy = (g*(densea-densep)/densep*qwo)/lnot
tds = salplo*0.64
!     Initialize lateral withdrawal flowrate for first/lowest cell in column/segment
jj = 0
qwd(laytop) = qwo
pwd(laytop) = 2.*(l+2.*bo)
awd(laytop) = l*bo
bwd(laytop) = sqrt(awd(laytop)/3.141592)
!
!    SOLUTION PROCEDURE
dz = 0.001
h = 0.001
nels = 0
hwith = 0.0
idm = idm +1
id = id - 1
30 z = z + dz
elev = elev - dz
nels = nels + 1
idm =idm - 1
if ((MOD(idm,500))==0) Then
id = id - 1
end if
!print *, id
!
!     Retrieve inner plume boundary conditions
 comgp = comgph(id)
!print *,comgph(1:id)
 cnmgp = cnmgph(id)
!print *,cnmgph(1:id)
tplumi = tplumh(id)
!print *,tplumh(1:id)
tamb = tambh(id)
!print *,tambh(1:id)
doamb = doambh(id)
!print *,doambh(1:id)
salpl = salplh(id)
!print *,salplh(1:id)
bi = bh(id)
!print *,bh(1:id)
vi = vh(id)
!print *,vh(1:id)
densea = densah(id)
!print *,densah(1:id)
densew = denswh(id)
!print *,denswh(1:id)
kolo = koloh(id)
!print *,koloh(1:id)
koln = kolnh(id)
!print *,kolnh(1:id)
ho2 = ho2h(id)
!print *,ho2h(1:id)
hn2 = hn2h(id)
!print *,hn2h(1:id)
po = poh(id)
!print *,poh(1:id)
pn = pnh(id)
!print *,pnh(1:id)
ei = eh(id)
!print *,eh(1:id)
rb = rbh(id)
!print *,rbh(1:id)
vb = vbh(id)
!print *,vbh(1:id)
co2 = co2h(id)
!print *,co2h(1:id)
 cn2 = cn2h(id)
!print *,cn2h(1:id)
!     Use subroutines for Runge Kutta method solution
y(1) = qwo 
y(2) = momeno
y(3) = ftempo
y(4) = fsalo
y(5) = fdoo
y(6) = fdno
y(7) = fgoo
y(8) = fgno
print *, ea,eo,ei
!print *, bo,bi,densew,densea,vo,vi
!print *, tplumo, tplumi ,tamb, salplo,denseo
print *,'do', doamb,co2 
!rb,vb,koln,ho2,po,comgp
call derivs_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, tplumo, &
 tplumi, salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, &
 ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, z, y, dydx)
call rk4_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, tplumo, &
 tplumi, salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, &
 ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, y, dydx, neqn, z, h, yout)
qwo = yout(1)
momeno = yout(2)
ftempo = yout(3)
fsalo = yout(4)
fdoo = yout(5)
fdno = yout(6)
fgoo = yout(7)
fgno = yout(8)
If (momeno>-0.0001) Then
tplumo = ftempo/qwo
salplo = fsalo/(qwo*densew)/(gamma/dense20)
!      Previous equation corrected to consistently express salinity in uS/cm
co2 = fdoo/qwo
 cn2 = fdno/qwo
Goto 40
End If
vo = momeno/qwo
print *, 'momeno,vo',momeno,vo
area = qwo/vo
!      SOLVE FOR DIMENSIONS USING L^2+(2Bo-Lo)L-AREA=0 USING QUADRATIC EQN.
aa = 1.0
bb = 2.*bnot - lnot
 cc = -1.0*area
l = (-1.0*bb+(bb**2-4.0*aa*cc)**(0.5))/(2.0*aa)
If (l<0.0) Then
l = (-1.0*bb-(bb**2-4.0*aa*cc)**(0.5))/(2.0*aa)
End If
bo = area/(2.0*l)
ea = -2.*bo*alphaa*vo
eo = -2.*bh(id)*alphao*vo
! add incremental entrainment to total cell entrainment/withdrawal
qwd(laytop+jj) = qwd(laytop+jj) + ea*dz + eo*dz - ei*dz
pwd(laytop+jj) = pwd(laytop+jj) + 2.*(l+2.*bo)
awd(laytop+jj) = awd(laytop+jj) + l*bo*2.
bwd(laytop+jj) = bwd(laytop+jj) + sqrt(l*bo*2./3.141592)
hwith = hwith + dz
If (hwith>hcell) Then
pwd(laytop+jj) = pwd(laytop+jj)/(2.*nels)
awd(laytop+jj) = awd(laytop+jj)/(2.*nels)
bwd(laytop+jj) = bwd(laytop+jj)/(2.*nels)
Print *, qwd(laytop+jj)!, bwd(laytop+jj), vo
jj = jj + 1
hwith = 0.0
nels = 0
End If
tplumo = ftempo/qwo
print *,'ftempo,qwo', ftempo, qwo
salplo = fsalo/(qwo*densew)/(gamma/dense20)
denseo = (0.059385*tplumo**3-8.56272*tplumo**2+65.4891*tplumo)*+999.84298 + (gamma)*salplo

If (vo<-1.E-5) Then
If (z<depth) Then
Goto 30
End If
End If

!calculation of average net oxygen mass transfer for day
40 grossmt = (fgonot-fgoo)*32./1000.*86400.
oteff = (fgonot-fgoo)/fgonot*100.

Print *, 'O3 Trans Eff:', oteff
elevob = elev
qwob = qwo
tplumb = tplumo
comgpb = & !*0.9 ! Added '*0.9' on 26 Nov to try to decrease
comgp
layob = laytop + jj
bob = bo
Return
End Subroutine doubplu_v1
!------------------------------------------------------------------------------
!
Subroutine rk4_2(e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, y, dydx, nn, x, h, yout)
Integer i, nn, nmax
Parameter (nmax=50)
Real *8 e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, h, x, dydx(nn), y(nn), yout(nn), h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)
External derivs_2
hh = h*0.5
h6 = h/6.
xh = x + hh
Do i = 1, nn
yt(i) = y(i) + hh*dydx(i)
End Do
Call derivs_2(e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, xh, yt, dyt)
Do i = 1, nn
yt(i) = y(i) + hh*dyt(i)
End Do
Call derivs_2(e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, xh, yt, dym)
Do i = 1, nn
yt(i) = y(i) + h*dym(i)
dym(i) = dyt(i) + dym(i)
End Do
Call derivs_2(e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, x+h, yt, dyt)
Do i = 1, nn
yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.*dym(i))
End Do
Return
End Subroutine rk4_2
!
!----------------------------------------------------------------------
!
Subroutine derivs_2(e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, x, y, dydx)
Real *8 e, densea, densew, g, b, tamb, salamb, gamma, dense20, doamb, pi, v, comgp, dnamb, cnmgp, x, y(8), dydx(8)
!     Right-hand side of differential equations for Runge-Kutta solution
dydx(1) = e
dydx(2) = 2.*pi*b**2*g*(densea-densew)/densew
dydx(3) = e*tamb
dydx(4) = e*(salamb*gamma/dense20)*densea
dydx(5) = e*doamb/32.
dydx(6) = e*dnamb/28.
Return
End Subroutine derivs_2
!----------------------------------------------------------------------
Subroutine rk4_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, tplumo, tplumi, &
salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, ho2, po, &
comgp, dnamb, koln, hn2, pn, cnmgp, y, dydx, nn, x, ho, yout)

Integer i, nn, nmax
Parameter (nmax=50)
Real *8 ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, tplumo, tplumi, &
salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, ho2, po, &
comgp, dnamb, koln, hn2, pn, cnmgp, ho, x, dydx(nn), y(nn), yout(nn), &
h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)
External derivs_3
hh = ho*0.5
h6 = ho/6.
xh = x + hh
Do i = 1, nn
yt(i) = y(i) + hh*dydx(i)
End Do
Call derivs_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, &
tplumo,tplumi, salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, &
ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, xh, yt, dyt)
Do i = 1, nn
yt(i) = y(i) + hh*dyt(i)
End Do
Call derivs_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, & 
tplumo, tplumi, salamb, salplo, salpl, gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, &
ho2, po, comgp, dnamb, koln, hn2, pn, cnmgp, xh, yt, dym)
Do i = 1, nn
yt(i) = y(i) + ho*dym(i)
dym(i) = dyt(i) + dym(i)
End Do
Call derivs_3(ea, eo, ei, densea, densew, denseo, densep, g, l, bi, bo, tamb, tplumo, &
tplumi,salamb, salplo, salpl,gamma, dense20, doamb, pi, rb, n, vi, vo, vb, kolo, &
 ho2, po,comgp, dnamb, koln,  hn2, pn, cnmgp, x+ho, yt, dyt)
Do i = 1, nn
yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.*dym(i))
End Do
Return
End Subroutine rk4_3
!----------------------------
SUBROUTINE derivs_3(ea,eo,ei,densea,densew,denseo,densep,g,l,bi,bo,tamb,tplumo, &
tplumi,salamb,salplo,salpl,gamma,dense20,doamb,pi,rb,n,vi,vo,vb,kolo, &
ho2,po,comgp,dnamb,koln,hn2,pn,cnmgp,x,y,dydx)
REAL*8 ea,eo,ei,densea,densew,denseo,densep,g,l,bi,bo,tamb,tplumo, &
tplumi,salamb,salplo,salpl,gamma,dense20,doamb,pi,rb,n,vi,vo,vb, &
kolo,ho2,po,comgp,dnamb,koln,hn2,pn,cnmgp,x,y(8),dydx(8)
!     Right-hand side of differential equations for Runge-Kutta solution
dydx(1)=(ea+eo-ei)/100000
dydx(2)=pi*g*(bo**2.-bi**2.)*(densew-densea)/1.1+(ei*vo-eo*vi)/100000
!dydx(3)=ei*tplumo-eo*tplumi-ea*tamb
dydx(3)=(eo*tplumi+ea*tamb-ei*tplumo)/100000
!dydx(4)=ei*(salplo*gamma/dense20)*denseo- &
!eo*(salpl*gamma/dense20)*densew &
!-ea*(salamb*gamma/dense20)*densea
dydx(4)=-1*ei/100000*(salplo*gamma/dense20)*denseo+ &
eo*(salpl*gamma/dense20)*densew/100000 &
+ea*(salamb*gamma/dense20)*densea/100000
dydx(5)=ei*doamb/32. -eo*co2 &
+4.0*pi*rb**2*n/(vi+vb)*koln*(ho2*po-comgp/32.) 
dydx(6)=ei*dnamb/28.-eo*cn2 &
+4.0*pi*rb**2*n/(vi+vb)*koln*(hn2*pn-cnmgp/28.) 
dydx(7)=0
dydx(8)=0
print *, 'dydx', dydx(1:8)
Return
End Subroutine derivs_3


