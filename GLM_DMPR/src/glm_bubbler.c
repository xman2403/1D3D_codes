/******************************************************************************
 *                                                                            *
 * glm_bubbler.c                                                              *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Earth & Environment                                          *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aed.see.uwa.edu.au/                                             *
 *                                                                            *
 * Copyright 2015, 2016 -  The University of Western Australia                *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  GLM is distributed in the hope that it will be useful,                    *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "glm.h"
#include "glm_globals.h"
#include "glm_util.h"
#include "glm_layers.h"
#include "glm_mixu.h"
#include "glm_deep.h"
#include "glm_input.h"
#include "glm_const.h"
#include "glm_bubbler.h"

#include "aed_time.h"
#include "aed_csv.h"

#define SQRT_Pi sqrt(Pi_)
#define SQRT_Two sqrt(2.0)

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#ifdef _VISUAL_C_
#define ALLOCATE1D_I(dst, cnt) dst=malloc(sizeof(int) * (cnt))
#define ALLOCATE1D_R(dst, cnt) dst=malloc(sizeof(AED_REAL) * (cnt))
#define ALLOCATE2D(dst, cnt1, cnt2) dst=malloc(sizeof(AED_REAL) * (cnt1) * (cnt2))
#endif


/******************************************************************************/

LOGICAL     bubbler_on        = FALSE;
char       *bubbler_data_file = NULL;
AED_REAL    bubbler_aflow     = 0.0;
int         bubbler_nports    = 0;
AED_REAL    bubbler_bublen    = 0.0;
AED_REAL    bubbler_bdepth    = 0.0;
LOGICAL     bubbler_opopt     = FALSE;
AED_REAL    bubbler_ton       = 0.0;
AED_REAL    bubbler_toff      = 0.0;
LOGICAL     bubbler_intopt    = FALSE;
char       *bubbler_start     = NULL;
char       *bubbler_stop      = NULL;
LOGICAL     bubbler_eff       = FALSE;
char       *bubbler_eff_file  = NULL;



//# Bubbler efficiency variables
static AED_REAL tMean, tMeanS;
static AED_REAL tEffic1, tEffic2, tEffic3, tPn, tTDiff, tHTRise, tFlowT;
static AED_REAL aEffic1, aEffic2,          aPn, aTDiff, aHTRise, aFlowT;

static AED_REAL timeOn;
static int  bSecs, nTimeSteps;
static FILE *ueff = NULL, *uair = NULL;

static AED_REAL PeI, PeF, xMeI, xMeF, zCI, zCF;
static AED_REAL cee, WISO, effic1, effic2, effic3;
static AED_REAL xPerI, xPerF, xVolR;
static AED_REAL htRise;

static AED_REAL lm2, lambda, pWidth, alp2, ub2d, ub2;

static const AED_REAL one = 1.0, ten = 10.0;

/*############################################################################*/

static AED_REAL aFlow;
static AED_REAL bDepth;
static AED_REAL bubLen;
static AED_REAL tOn;
static AED_REAL tOff;
static AED_REAL cdifBeg;

static int nPorts;
static int iDayStart;
static int iDayFinish;
static int bLevel;

static int ibub;
static LOGICAL opopt, eff, intopt, bubOn, bubprint;

static char *fair = NULL;
static char *feff = NULL;

#define nStepY   10

/*============================================================================*/
static void pe3(AED_REAL *_xp, AED_REAL *_zc, AED_REAL *_xmas, AED_REAL BDE,
                                       AED_REAL *_xper, AED_REAL *_xvolr);

static void inBubble(LOGICAL intopt, int MaxNSy);

static void InActRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM);

static void InDepRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, AED_REAL E, AED_REAL DENSEA,
              AED_REAL DENSEW, AED_REAL DENSEP, AED_REAL LAMBDA, AED_REAL L,
              AED_REAL TAMB, AED_REAL DOAMB, AED_REAL DNAMB, AED_REAL RB,
              AED_REAL NN, AED_REAL V, AED_REAL VB, AED_REAL KOLO, 
              AED_REAL KOLN, AED_REAL HO2, AED_REAL HN2, AED_REAL PO,
              AED_REAL PN, AED_REAL COMGP, AED_REAL CNMGP);

static void InActDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM);

static void InDepDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, AED_REAL E, AED_REAL DENSEA,
              AED_REAL DENSEW, AED_REAL DENSEP, AED_REAL LAMBDA, AED_REAL L,
              AED_REAL TAMB, AED_REAL DOAMB, AED_REAL DNAMB, AED_REAL RB,
              AED_REAL NN, AED_REAL V, AED_REAL VB, AED_REAL KOLO, 
              AED_REAL KOLN, AED_REAL HO2, AED_REAL HN2, AED_REAL PO,
              AED_REAL PN, AED_REAL COMGP, AED_REAL CNMGP);

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_bubbler()
{
    int secs1, secs2;
    bSecs = 0;

    //# Variables for changes to the efficiency files
    nTimeSteps = 0;
    bubprint = FALSE;

    ibub = -1;

    //# Initialise bubbler details
    intopt = FALSE;
    timeOn = 0.0;
    aFlow = 0.0001;
    bDepth = 0.0;
    bubLen = 1.0;
    nPorts = 1;

    bubOn = TRUE;   // 13/12/19

    if (bubOn) {
        printf("Bubbler is ON\n");
        bubprint = TRUE;

        bubOn  = bubbler_on;
        fair   = bubbler_data_file;
        aFlow  = bubbler_aflow;
        nPorts = bubbler_nports;
        bubLen = bubbler_bublen;
        bDepth = bubbler_bdepth;
        opopt  = bubbler_opopt;
        tOn    = bubbler_ton;
        tOff   = bubbler_toff;
        intopt = bubbler_intopt;
        read_time_string(bubbler_start, &iDayStart, &secs1);
        read_time_string(bubbler_stop, &iDayFinish, &secs2);
        eff    = bubbler_eff;
        feff   = bubbler_eff_file;

        fair = "bubbler.csv";   // 14/12/19

        iDayStart = 2457174;
        iDayFinish = 2457176;

        printf(" Air Flow = %10.4f No. Ports = %5d Time On = %8d days\n",    // 14/12/19
                                  aFlow, nPorts, iDayStart);
        if (opopt) printf ("opopt check") ; // 14/12/19

        ibub = 0;

        //# Bubbler efficiency file setup
        if (eff) {
            ueff = fopen(feff, "w");

            fprintf(ueff, "JDAY, EFFECT, T DIFF, C, HT RISE, aFlow, EFFIC\n");
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void do_bubbler(int jday, int iclock)
{
    AED_REAL tDiff;
    int i;
    
    bubOn = TRUE;      // 17/12/19
  //  ibub = 0;

    if ( ! bubOn ) return;

    //# Bubbler calculations
    if (ibub >= 0) {
        /**********************************************************************
         * Install bubble, get pe etc                                         *
         * Test for operation of bubbler                                      *
         * Test that JDAY falls within period of bubbler operation            *
         * Calculate the time step                                            *
         **********************************************************************/
        for (i = surfLayer; i > botmLayer; i--) {
            if (Lake[i].Height >= bDepth) bLevel = i;
        }

        if (bLevel == surfLayer)
            tDiff = 0.0;
        else
            tDiff = Lake[surfLayer].Temp - Lake[bLevel+1].Temp;

        
        printf ("before optimisation");

        opopt = FALSE; // 14/12/19

        /**********************************************************************
         * Bubbler optimisation                                               *
         * if tDiff GT Delta T (BUBBLER ON) { TURN BUBBLER ON                 *
         * if tDiff LT Delta T (BUBBLER OFF) { TURN BUBBLER OFF               *
         * } else { LEAVE AS FOR PREVIOUS DAY                                 *
         * Note: pe's from s/r pe3 are based on area calculations             *
         *       from TABLE.  Units of area (in TABLE) are 10**6 sq.m         *
         *       Therefore this multiplier needs to be used on                *
         *       everything that comes out of pe3                             *
         **********************************************************************/
          if (jday >= iDayStart && jday <= iDayFinish) {
              if (opopt) {
                  if (! bubOn && tDiff > tOn) bubOn = TRUE;
                  if (tDiff < tOff) bubOn = FALSE;
              } else {
                  bubOn = TRUE;
              }
              if (aFlow <= 0.0) bubOn = FALSE;
          } else {
              bubOn = FALSE;
          }

       if (bubOn) printf ("after optimisation");

        if (iclock == 0) {
            bSecs = 0;
            tEffic1 = 0.0;
            tEffic2 = 0.0;
            tEffic3 = 0.0;
            tPn = 0.0;
            tTDiff = 0.0;
            tHTRise = 0.0;
            tFlowT = 0.0;
        }

        if (bubOn) {
            effic1 = 0.;
            effic2 = 0.;
            effic3 = 0.;
            ibub = ibub+1;

            pe3(&PeI, &zCI, &xMeI, bDepth, &xPerI, &xVolR);

            if (ibub == 1)
                cdifBeg = PeI;

            inBubble(intopt, nStepY*MaxLayers);

            check_layer_stability();
            check_layer_thickness();
            pe3(&PeF, &zCF, &xMeF, bDepth, &xPerF, &xVolR);

            //# Change in PE as we proceed
            //# Mechanical efficiency
            //# Isothermal compression work of compressor
            WISO = 2.303 * (1.0) * aFlow * noSecs *
                              log10((Lake[surfLayer].Height+10.2)/10.2) * 1.E5;
            effic2 = (xPerF-xPerI) * 9.81 * 1.E6 / WISO;

            //# change in pe based on initial stratification
            //# effic3 = (xPerF*1.e6-peibeg)/(pembeg-peibeg)
            //# try using difference between c.o.m and c.o.v
            effic3 = PeF/cdifBeg;

            bSecs = bSecs + noSecs;
            tEffic1 = tEffic1 + effic1 * (noSecs);
            tEffic2 = tEffic2 + effic2 * (noSecs);
            tEffic3 = tEffic3 + effic3 * (noSecs);
            tPn = tPn + cee * (noSecs);
            tTDiff = tTDiff + tDiff * (noSecs);
            tHTRise = tHTRise + htRise * (noSecs);
            tFlowT = tFlowT + aFlow * (noSecs);

            //# Calculate mean and variance of tDiff
            timeOn = timeOn + (noSecs) / 86400.0;
            nTimeSteps = nTimeSteps + 1;
            tMean = tMean + tDiff;
            tMeanS = tMeanS + tDiff * tDiff;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_bubbler(int jday)
{
    //# Input from bubbler file if required
    if ( bubOn && jday >= iDayStart && jday <= iDayFinish) {
        read_bubble_data(jday, &aFlow, &nPorts, &bDepth, &bubLen);
        if (aFlow <= 0.0) bubOn = FALSE;
    }

   printf(" Air Flow = %10.4f No. Ports = %5d Time On = %8d days\n",
                                  aFlow, nPorts, jday);   // 21/01/20 jday

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void write_bubbler(int jday)
{
    if ( ! bubOn ) return;

    //# Bubbler efficiency file writing
    //# Write daily bubbler stuff if required
    if (eff && ueff!=NULL) {
        aEffic1 = tEffic1/(bSecs);
        aEffic2 = tEffic2/(bSecs);
        aPn = tPn/(bSecs);
        aTDiff = tTDiff/(bSecs);
        aHTRise = tHTRise/(bSecs);
        aFlowT = tFlowT;

        fprintf(ueff, "%d, %12.3f, %12.3f, %12.3f, %12.3f, %12.3f, %12.3f\n",
                          jday, aEffic1, aTDiff, aPn, aHTRise, aFlowT, aEffic2);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void print_bubbler()
{
    if ( ! bubOn ) return;

    if (bubprint && eff) {
        //# Check whether an efficiency file exists, and is opened
        //# If so, then need to write out all the data read in
        if ( ueff != NULL ) {
            //# Check whether writing bubbler data from keyboard or file
            if (uair != NULL )
                fprintf(ueff,
                       "Bubbler Information\nInput File is \"%s\" Time On = %8.2f days\n",
                                  fair, timeOn);
            else
                fprintf(ueff, " Air Flow = %10.4f No. Ports = %5d Time On = %8.2f days\n",
                                  aFlow, nPorts, timeOn);
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void done_bubbler()
{
    if ( ! bubOn ) return;

    fclose(uair);
    fclose(ueff);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *  This gets the PE stored above the level of BDEP (Includes the layer       *
 *  containing BDEP)                                                          *
 *                                                                            *
 * NOTE: All totals are relative to BDEP                                      *
 *                                                                            *
 * Variable Dictionary                                                        *
 *                                                                            *
 *     XP=(Centre of mass-Centre of vol)*g                                    *
 *     ZC=Centre of volume                                                    *
 *     XMAS=Total mass above bubbler                                          *
 *     XVOL=Total vol above bubbler                                           *
 *     XPE=Potential energy/g                                                 *
 *                                                                            *
 ******************************************************************************/
void pe3(AED_REAL *_xp, AED_REAL *_zc, AED_REAL *_xmas, AED_REAL BDE,
                                       AED_REAL *_xper, AED_REAL *_xvolr)
{
    AED_REAL  XP = *_xp, ZC = *_zc, XMAS = *_xmas;
    AED_REAL  XPER = *_xper, xVolR = *_xvolr;

    DOUBLETYPE X, Y, ZCVP;
    DOUBLETYPE AB, HB, DA, HT, DENS;
    DOUBLETYPE YB, XVOL, XVOLP, XMASSP;
    DOUBLETYPE XPE, ZCV, XMASS, BDEP;

    int i, ij, il, ix, L, IBDEP, ik;

/*----------------------------------------------------------------------------*/

    IBDEP = 0;
    DENS = 0.;

    BDEP = BDE;

    //# Get the layer containing BDEP
    for (i = surfLayer; i > botmLayer; i--)
        if (Lake[i].Height > BDEP) IBDEP = i;

    /**************************************************************************
     * Get centre of volume and potential energy, total mass above bdep etc   *
     * locate surface in relation to storage table                            *
     * ij is the last storage entry before the surface, y is the distance     *
     * to the surface from there                                              *
     **************************************************************************/
    X = Lake[surfLayer].Height * ten;
    Y = AMOD(X, one);
    ij = (X - Y);
    if (ij >= X) ij--;
    if (ij > Nmorph) {
        Y += (ij - Nmorph);
        ij = Nmorph;
    }
    ij--; // offset for 0 based index
    Y /= ten;

    //# Locate BDEP in relation to storage table
    //# ik is the last storage entry before BDEP, YB is the distance to
    //# BDEP from there; NOTE ik could = 0
    X = BDEP * 10.0;
    YB = AMOD(X, one);
    ik = (X - YB);
    if (ik > Nmorph) {
        YB += (ik - Nmorph);
        ik = Nmorph;
    }
    ik--; // offset for 0 based index
    YB /= 10.0;

    //# Initialise the loop, il is the current layer number,
    //# DA the area gradient
    //# starting at IBDEP
    //#
    ZCV = 0.0;
    XPE = 0.0;
    XMASS = 0.0;
    XVOL = 0.0;
    il = IBDEP;
    ix = 0;
    if (ik < 0) {
        DA = (MphLevelArea[0]) * 10.0;
        AB = DA * YB;
    } else {
        DA = (dMphLevelArea[ik]) * 10.0;
        AB = (MphLevelArea[ik]) + DA * YB;
    }
    HB = 0.0;

    //# Loop
    if (ik != Nmorph) {
        ik++;
        for (i = ik; i < ij; i++) {
            ix = 1;
            while (ix == 1) {
                DENS = (Lake[il].Density);
                HT = (i)/10.0 - BDEP;
                if (i != 0) DA = (dMphLevelArea[i-1])*10.0;
                ix = 0;
                if (HT>=((Lake[il].Height) - BDEP)) {
                    HT = (Lake[il].Height) - BDEP;
                    il = il+1;
                    ix = 1;
                }
                ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
                XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
                XMASSP = XVOLP*(1000.0+DENS);
                ZCV = ZCV+ZCVP;
                XPE = XPE+(1000.0+DENS)*ZCVP;
                XMASS = XMASS+XMASSP;
                XVOL = XVOL+XVOLP;
                AB = DA*(HT-HB)+AB;
                HB = HT;
            }
        }
    }

    //# This takes us to the last storage table entry
    //# Now go to the surface
    if (il == surfLayer) {
        HT = (Lake[surfLayer].Height) - BDEP;
        ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
        XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
        XMASSP = XVOLP*(1000.0+DENS);
        ZCV = ZCV+ZCVP;
        XPE = XPE+(1000.0+DENS)*ZCVP;
        XMASS = XMASS+XMASSP;
        XVOL = XVOL+XVOLP;
    } else {
        for (L = il; L <= surfLayer; L++) {
           DENS = (Lake[L].Density);
           HT = (Lake[L].Height) - BDEP;
           ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
           XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
           XMASSP = XVOLP*(1000.0+DENS);
           ZCV = ZCV+ZCVP;
           XPE = XPE+(1000.0+DENS)*ZCVP;
           XMASS = XMASS+XMASSP;
           XVOL = XVOL+XVOLP;
           AB = (Lake[L].LayerArea);
           HB = HT;
       }
    }

    //# Note that the units of XP are millions of kg M**2/sec**2,
    //# units of XMAS are millions of kg
    ZC = (ZCV/XVOL);
    XP = ((XPE/XMASS-ZCV/XVOL)*9.810);
    XMAS = (XMASS);
    XPER = (XPE);
    xVolR = (XVOL);

    // set returned values
    *_xp = XP; *_zc = ZC; *_xmas = XMAS;
    *_xper = XPER; *_xvolr = xVolR;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int inner_loop(int iLevel, int level, int *pCountI, int *retI, int NSY, int *iTop,
                                AED_REAL *yDepth, AED_REAL *yDepthM, AED_REAL *TH,
                                AED_REAL *W, AED_REAL *B, AED_REAL *ENT, AED_REAL *yDen)
{
    const AED_REAL SQLAM=0.09, US=0.3E0, HPRES=10.2E0;
    const AED_REAL R2=1000., ALP=0.083;
   //# const int threeI=3;
    const int threeI = 8; // 28/04/20

     const AED_REAL DIAMM = 1.14, QGFRAC = 1.0, ALPHA = 0.083, G = 9.80665;
     const AED_REAL GAMMA = 6.9E-4, LAMBDA=0.93, FRNOT=2.0 , PSTD=101325., PATM = 96261.0;
     const AED_REAL   RGAS=8.314, TSTD=293.15 ,DENSE20=998.2;

    int DIM0, DIM;

    AED_REAL GRADMO;
 //#   AED_REAL DYDX[3], Y[3], YOUT[3];

    AED_REAL DYDX[8], Y[8], YOUT[8];

    AED_REAL BIGH, BIGM, Z, DZ, SMALLZ;
    AED_REAL B0, W0, TH0;
    AED_REAL LPP;
    AED_REAL Q0, BU0, XM0, DQ, DRH, DW;
    AED_REAL B2, W2, TH2, XM2, XM1, BU01, BU02;
    AED_REAL DQDZ, DXMDZ, DBUDZ, ZZ, DW0DZ;
    AED_REAL PA, UB, HLAY;
    AED_REAL FLOW;
    AED_REAL RHO0;

     AED_REAL DENSEA, DENSEW, TAMB, QSCFM, QSCMS, QGAS, RB, RBNOT, VB, VG, VGUESS, BW, L, VDIFF;
     AED_REAL DENSEP, E, V, QW, MOMENT, AA, BB, CC, AREA, VBUB, NN, DOAMB, DNAMB, KOLO,KOLN,HO2,HN2;  //# B is changed to BW to avoid redifination
     AED_REAL PO, PN, COMGP, CNMGP, CO2, COMGNOT, TPLUME, CN2, CNMG, FTEMP, FSAL, FDO, FDN, FGO;
     AED_REAL FGONOT, FGN, YO2, YN2, PZ, FRACO, FRACN, DCO2, FR; 

    int i, pCount = *pCountI;    

    pCount = 1;

///#--------     02/04/20  update density 
 
     FRACO = 0.8;
     FRACN = 0.2;   //#  28/04/20
     COMGP = 300;  //# _WQ_Vars(jj, ii)     
     DOAMB = COMGP;
     CO2 = COMGP/32;
     COMGNOT = COMGP;
//#------
      BW = 0.012/0.93;
      QSCFM = 25;
      L = 150 + 2.0*BW*0.07;
         HLAY = yDepth[NSY] - bDepth;       //# 11/03/20
      DENSEA = Lake[level].Density;  //#RHO0?

      if (DENSEA <= 992)
          DENSEA = 999.43;

      DENSEW = 999.43;
      TAMB = Lake[level].Temp;
  
      if (TAMB <= 8 ) 
          TAMB = 15;
       
      TPLUME = TAMB;
//#-------     BUBBLE PROPERTIES
//# C     Gas flow rate per segment asssumed to be proportional to fraction of total diffuser length.
      QSCMS=pow(QGFRAC*QSCFM/3.281,3)/60.0;
      QGAS=PSTD*QSCMS*(TAMB+273.15)/((PATM+DENSEA*G*5.5)*TSTD);  //# 5.5 is diffuser depth
//# C     For diffuser in SHR, use correlation by McGinnis and Little (2000) for initial bubble size.     
      RB=DIAMM/2000.;
      RBNOT=RB;
        if (RB <= (7.5E-4))
              VB=pow(1189.0*RB,1.1945);
        else { if (RB >= (7.5E-4) && RB <= (4.8E-3))
              VB=0.22;
        else
              VB=2.995*pow(RB,0.489);}
//#--------
    KOLO = 0.6*RB;
    if (KOLO > (4.0E-4))
         KOLO = 4.0E-4;
    KOLN = KOLO;

    HO2 = (2.125 - 0.05023*TPLUME+5.7714E-4*pow(TPLUME,2))/100000.;
    HN2 = (1.042 - 0.02457*TPLUME+3.1714E-4*pow(TPLUME,2))/100000.;

    CN2 = (PATM*0.79)*HN2;
    CNMG = CN2*28.0;
    CNMGP = CNMG;
    DNAMB = CNMG;

      VGUESS = 0.07;
//#--------
      VBUB = 4./3.*3.1415*pow(RB,3);
      NN = QGAS/VBUB;
//#--------
  do {   VG=QGAS/((VGUESS+VB)*(2.*LAMBDA*BW)*(L-2.0*BW*(1.0-LAMBDA)));
      DENSEP=(1.0-VG)*DENSEW;
        V=FRNOT*pow((2.0*LAMBDA*BW*G*(DENSEA-DENSEP)/DENSEP),0.5);     
     VDIFF=abs(V-VGUESS);
     if (VDIFF > 1.0E-6)
         VGUESS=V;
          else break;
     }   while (1);
         
        V = 0.083;  //# 12/03/20
        DENSEP = (1.0-VG)*DENSEW;
        DENSEP = 992.6;  //# 12/03/20
//#---------

//  HLAY = Lake[surfLayer].Height - BD;

    *retI = level;

    LPP = bubLen/(2.0*(nPorts));
    RHO0 = Lake[iLevel].Density;
    PA = RHO0 * g * HPRES;
    FLOW = aFlow/(nPorts);
    UB = US * (1.E0 + SQLAM);

    /********************%%%%%%%%%%%%%%%%%%%%%%%%%%%%**********************/
    //# Start of outer loop
    //# Initilize parameters for plume
    do  {         //  18/12/19
        bDepth = yDepth[level];

        DIM = 1;
    
        //# Origin of bubbles is level; locate depth and set initial conditions
        //# Based on series solutions McDougall
        HLAY = yDepth[NSY] - bDepth;

 //       HLAY = 4.5;    // 16/12/19

        BIGH = HLAY + HPRES;

   //     BIGH = 14.70;  // 16/12/19        

        BIGM = FLOW * PA * (1.E0 + SQLAM) /
                       (4.E0 * Pi * SQR(ALP) * RHO0 * SQR(BIGH) * SQR(UB));

       // BIGM = 545.53;  //  16/12/19  for flow 10
        

        //#  Change initial height to 1 layer depth
        if (level != NSY)
            Z = yDepth[level+1]-yDepth[level];
        else
            Z = yDepth[level]-yDepth[level-1];

      //  Z =  0.1 ;   // 16.12.19
        
        ZZ = Z/BIGH;
        B0 = 2.E0 * ALP * Z *
               (0.6E0 + 0.01719E0 * pow(ZZ/BIGM, 0.333E0) -
                  0.002527E0 * pow(ZZ/BIGM, 0.667E0) +
                  ZZ * (-0.04609E0+0.000031E0/BIGM));
        W0 = UB * pow(BIGM/ZZ, 0.333E0) * (1.609E0 - 0.3195E0 * pow(ZZ/BIGM, 0.333E0) +
              0.06693E0 * pow(ZZ/BIGM, 0.667E0) + ZZ * (0.4536E0 - 0.0105E0 / BIGM));

        if (intopt)  //# NOTE TH0 no longer has g in it
            TH0 = FLOW * PA * (1.E0+SQLAM)/(SQLAM*(1.-ZZ)*BIGH*Pi*RHO0*(W0+UB)*B0*B0);
        else
            TH0 = FLOW * PA * (1.E0+SQLAM)/(SQLAM*(1.-ZZ)*BIGH*Pi*RHO0*(W0+UB)*B0*B0*g);

        //# Use 0 for initial gradient
        DW0DZ = 0.0;

        if (intopt) {
            //# Initial check for plume interaction
            if (B0>=LPP)
                DIM = 2;
            else
                DIM = 1;

            //#  **  Initial conditions for line source conditions
            //#       however presently it is assumed to be similar to point source
            //#       but if aerator is set up with single ports this should be ok
            //#       because it will not go through these calculations anyhow.
            //#       should be modified to be consistent.
            if (DIM==2) {
                BIGM = FLOW * PA * (1.0+lm2) / (4.0 * Pi * SQR(alp2) * RHO0 * SQR(BIGH) * CUBE(UB));
                B0 = 2.0 * alp2 * Z * (0.6E0 + 0.01719E0 * pow(ZZ/BIGM, 0.333E0) -
                       0.002527E0 * pow(ZZ / BIGM, 0.667E0) + ZZ * (-0.04609E0 + 0.000031E0 / BIGM));
                W0 = ub2 * pow(BIGM/ZZ, 0.333E0) * (1.609E0 - 0.3195E0 *
                       pow(ZZ/BIGM, 0.333E0) + 0.06693E0 * pow(ZZ/BIGM, 0.667E0) +
                                  ZZ * (0.4536E0 - 0.0105E0 / BIGM));

                //# NOTE TH0 no longer has g in it
                TH0 = FLOW*PA*(1.0+lm2)/(lm2*(1.-ZZ)*BIGH*Pi* RHO0*(W0+UB)*B0*B0);

                Q0 = B0*W0*pWidth* sqrt(Pi);
                BU0 = RHO0*lambda*B0*TH0* sqrt(Pi)*W0*(pWidth)/sqrt(1.E0+lm2);
                XM0 = RHO0*W0*W0*sqrt(Pi)*B0*(pWidth)/sqrt(2.E0);
            } else {
                //# Initial fluxes

                Q0 = B0 * B0 * W0 * Pi;
                BU0 = RHO0 * SQLAM * B0 * B0 * TH0 * Pi * W0 / (1.E0 + SQLAM);
                XM0 = RHO0 * W0 * W0 * Pi * B0 * B0 / 2.E0;
            }
        } else {
            //# Initial fluxes
            Q0 = B0*B0*W0;
            BU0 = Q0*TH0;
            XM0 = Q0*W0;
        }

         E=2.*(L+2.*BW)*ALPHA*V;    
         QW=2.*L*BW*V;
         MOMENT=2.*L*BW*pow(V,2);
         FTEMP = QW*TPLUME;
         FSAL = 0;
         FDO = QW*CO2;
         FDN = QW*CN2;
         FGO = PSTD*QSCMS/(RGAS*TSTD)*0.21;
         FGONOT = FGO;
         FGN = PSTD*QSCMS/(RGAS*TSTD)*0.79;

         YO2 = FGO/(LAMBDA*2.*BW*(L-2.*BW*(1.-LAMBDA))*(V+VB));
         YN2 = FGN/(LAMBDA*2.*BW*(L-2.*BW*(1.-LAMBDA))*(V+VB));
         PZ = PATM+(DENSEA*9.8* bDepth);
         PO = PZ*FRACO;
         PN = PZ*FRACN;
     //# BUOY = (G*(DENSEA-DENSEP)/DENSEP*QW)/150; //#LNOT = 150 02/04/20

        // pCount = 1; // 14/12/19

        //# Start inner loop

      //  iTop[pCount] = level; // 14/12/19
 
        iTop[*pCountI] = level;

        for (i = level; i <= NSY; i++) {
            if (W0 <= 0.0) {   //# W0 is v  30/01/20
                if (i == level) return FALSE;
                break;
            }

            if (i == 1) {
                DZ = yDepth[0];
                DRH = (yDen[1] - yDen[0]) / (yDepthM[1] - yDepthM[0]);
            } else {
                DZ = yDepth[i] - yDepth[i-1];
           
            if (DZ <= 0)
                DZ = 0.1;   //# 12/03/20
            //#    DZ = 0.1;   //#  19/12/19

                if (i == NSY) {
                    DRH = (yDen[NSY]-yDen[NSY-1]) / (yDepthM[NSY]-yDepthM[NSY-1]);
                } else {
                    DRH = (yDen[i+1]-yDen[i-1]) / (yDepthM[i+1]-yDepthM[i-1]);
            //       DRH = 0.18;  // 19/12/19
                }
            }

            SMALLZ = BIGH - (yDepth[i] - bDepth);
            DRH = DRH / RHO0;

            //#  Calculate initial derivatives
            if (intopt) {
                if (DIM == 2) {
                    DQDZ = 2.0 * alp2 * pWidth * W0;
                    DXMDZ = sqrt(Pi) * RHO0 * TH0 * lambda * B0 * pWidth;
                    DBUDZ = RHO0 * sqrt(Pi) * B0 * pWidth * W0 * g * DRH + FLOW*PA *
                              (DW0DZ*(1.E0-(W0 / (W0+ub2d)))+(W0/SMALLZ))/
                              (SMALLZ*(W0+ub2d));
                } else {
                    DQDZ = 2.0 * ALP * B0 * W0 * Pi;
                    DXMDZ = Pi * RHO0 * TH0 * SQLAM * B0 * B0;
                    DBUDZ = RHO0 * Pi * B0 * B0 * W0 * g * DRH + FLOW * PA *
                              (DW0DZ*(1.-W0/(W0+UB))+ W0/SMALLZ)/(SMALLZ*(W0+UB));
                }
            } else {
                DQDZ = 2.0 * ALP * B0 * W0;
                DXMDZ = 2.0 * g * SQLAM * B0 * B0 * TH0;
                DBUDZ = (1.0 + SQLAM) * DRH * Q0 / SQLAM +
                          ((1.+SQLAM)*FLOW*PA/(g*Pi*SQLAM*RHO0*SMALLZ*(W0+UB))) *
                                             (W0/SMALLZ + (DW0DZ)*(1.-W0/(W0+UB)));
            }


//#
//#

            DYDX[0] = DQDZ;
            DYDX[1] = DXMDZ;
            DYDX[2] = DBUDZ;
//#
//#
           DYDX[0]=E;         
           DYDX[1]=(DENSEA-DENSEW)/DENSEP*G*L*2.*BW+(DENSEW-DENSEP)/
        +DENSEP*G*LAMBDA*2.*BW*(L-2.*BW*(1.-LAMBDA));
           DYDX[2]=E*TAMB;
           DYDX[3]=0;
           DYDX[4]=(E*DOAMB/32.+4*3.14*RB*pow(RB,2)*NN/(V+VB)*KOLO*(HO2*PO-COMGP/32.));
           DYDX[5]=(E*DNAMB/28.+4*3.14*RB*pow(RB,2)*NN/(V+VB)*KOLN*(HN2*PN-CNMGP/32.)); 
           DYDX[6]= -4*3.14*pow(RB,2)*NN/(V+VB)*KOLO*(HO2*PO-COMGP/32.);
           DYDX[7]= -4*3.14*pow(RB,2)*NN/(V+VB)*KOLN*(HN2*PN-CNMGP/28.);

//#
          Y[0]=QW;
          Y[1]=MOMENT;
          Y[2]=FTEMP;
          Y[3]=0;
          Y[4]=FDO;
          Y[5]=FDN;
          Y[6]=FGO;
          Y[7]=FGN;
   
//#
         //#   Y[0] = Q0;
         //#   Y[1] = XM0;
         //#   Y[2] = BU0;

            if (intopt) {
                if (i != level){
                    if (DIM==2) {
                        W0 = sqrt(2.0)*Y[1]/(RHO0*Y[0]);
                        TH0 = Y[2]*(sqrt(1.0+lm2))/(Y[0]*lambda*R2);
                    } else {
                        W0 = 2.0*Y[1]/(RHO0*Y[0]);
                        TH0 = Y[2]*(1.0+SQLAM)/(SQLAM*Y[0]*R2);
                    }
                }
            } else {
                if (i != level){
                    W0 = Y[1]/Y[0];
                    TH0 = Y[2]/Y[0];
                }
            }

            B0 = BW;    // 11/03/20 

            if (intopt)
                InActRK4(Y, DYDX, threeI, DZ, YOUT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);
            else
                InDepRK4(Y, DYDX, threeI, DZ, YOUT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, E, DENSEA, DENSEW, DENSEP,LAMBDA,L,TAMB, DOAMB, DNAMB, RB, NN, V, VB, KOLO,  KOLN, HO2, HN2, PO, PN, COMGP, CNMGP);

            if (YOUT[1] <= 0.000001) break;

            Q0 = YOUT[0];
            XM0 = YOUT[1];
            BU0 = YOUT[2];

//#
           QW=YOUT[0];
           MOMENT=YOUT[1];
           FTEMP = YOUT[2];
           FSAL = YOUT[3];
           FDO = YOUT[4];
           FDN = YOUT[5];
           FGO = YOUT[6];
           FGN = YOUT[7];

            DQ = YOUT[0]-Y[0];

            //# Use test for change in sign of momentum flux and also
            //# rate of change of momentume flux.  the latter is necessary as
            //# erroneous results will be produced if the momentum just falls
            //# short of going negative
            GRADMO = Y[1]/YOUT[1];
            GRADMO = fabs(GRADMO);

            if (XM0 > 0.000001 && Q0 > 0.000001) {
                if (intopt) {
                    if (DIM==2)
                        W[i] = sqrt(2.0)*XM0/(RHO0*Q0);
                    else
                        W[i] = 2.0*XM0/((RHO0)*Q0);
                } else
                    W[i] = XM0 / Q0;
            } else
                break;

            if (GRADMO > 1.5)
                break;

            if (W[i] > 0.000001 && Q0 > 0.000001) {
                if (intopt) {
                    if (DIM==2)
                        B[i] = Q0/(W[i]*pWidth*sqrt(Pi));
                    else
                        B[i] = sqrt(Q0/(W[i]*Pi));
                } else
                    B[i] = sqrt(Q0 / W[i]);
            } else //#  Change to avoid divison by zero and sqrt(-ve)
                B[i] = B0;  //# break here? 30/01/20

            DW = W[i] - W0;


//#-------------
         V = MOMENT/QW;

         if (V <= 0.000001)
          break;

         AREA = QW/V;
         AA=1.0;
         BB=2.*(5/0.93)-(10+2*(5/0.93)*0.07);
         CC=-1.0*AREA;
         L=(-1.0*BB+pow((pow(BB,2)-4.0*AA*CC),0.5))/(2.0*AA);
         if(L < 0.0)
             L=(-1.0*BB-pow((pow(BB,2)-4.0*AA*CC),0.5))/(2.0*AA);
         BW =AREA/(2.0*L);
         E = 2*ALPHA*(BW+2*L)*V; 

         TPLUME = FTEMP/QW;
         CO2=FDO/QW;
         CN2=FDN/QW;
         COMGP = CO2*32;
         CNMGP = CN2*28;
         YO2 = FGO/(LAMBDA*2.*BW*(L-2.*BW*(1-LAMBDA))*(V+VB));
         YN2 = FGN/(LAMBDA*2.*BW*(L-2.*BW*(1-LAMBDA))*(V+VB)); 

          if (yDepth[i] == 0)
              yDepth[i] = 5;
         PZ = PATM*(DENSEA*G*yDepth[i]);
         QGAS = (FGO+FGN)*RGAS*(TPLUME+273.15)/PZ;
         VBUB = QGAS/NN;

         VG=VBUB*NN/((V+VB)*(2.*LAMBDA*BW)*(L-2.*BW*(1.-LAMBDA)));      
         RB=pow((3.*QGAS/(4.*3.1415*NN)),(1./3.));
          if(RB < 0.0)
             RB=1.0E-8;
         FRACO=FGO/(FGO+FGN);
         FRACN=1.0-FRACO;
	
         PO=PZ*FRACO;
         PN=PZ*FRACN;
         DENSEA=Lake[i].Density;
        if (DENSEA < 995)
              DENSEA = 999.5;
         DENSEP=(1.0-VG)*DENSEW;


         if(RB <= (7.5E-4))
              VB=1189.0*pow(RB,1.1945);
         else { if((RB > (7.5E-4)) & (RB <(4.8E-3)))
                  VB=0.22;
                else
                  VB=2.995*pow(RB,0.489);
               }

         KOLO=0.6*RB;
         if(KOLO > (4.0E-4))
               KOLO=4.0E-4;
          KOLN=KOLO;

         HO2=(2.125-0.05023*TPLUME+5.7714E-4*pow(TPLUME,2))/100000.;
         HN2=(1.042-0.02457*TPLUME+3.1714E-4*pow(TPLUME,2))/100000.;

        FR=V/pow((2.*LAMBDA*BW*G*(DENSEA-DENSEP)/DENSEP),0.5);
        DCO2=HO2*PO-CO2;

       //#  DENSEP = 992.6 +  ((yDepth[NSY] - yDepth[i])

//#     DENSEA DENSEP  //# 24/02/20
            if (intopt) {
                if (Q0 != 0.0) {
                    if (DIM==2)
                        TH[i] = BU0*(sqrt(1.0+lm2))/(Q0*lambda*R2);
                    else
                        TH[i] = BU0*(1.0+SQLAM)/(SQLAM*Q0*R2);
                } else {
                    if (DIM==2)
                        TH[i] = BU0*(sqrt(1.0+lm2))/(Y[0]*lambda*R2);
                    else
                        TH[i] = BU0*(1.0+SQLAM)/(SQLAM*Y[0]*R2);
                }

                if (i==level)
                    ENT[i] = E;
                else
                    ENT[i] = E;

                DIM0 = DIM;

                if (B[i]>=LPP) DIM = 2;

                if (DIM0==1 && DIM==2) {
                    B2 = (sqrt(2.0*Pi))*(SQR(B[i]))/pWidth;
                    W2 = W[i]/sqrt(2.0);
                    TH2 = TH[i] * SQLAM * (sqrt(1.0 + lambda * lambda))/(lambda * (SQLAM + 1.0));
                    XM2 = W2 * W2 * sqrt(Pi) * B2 * pWidth / sqrt(2.0);
                    XM1 = W[i] * W[i] * Pi * B[i] * B[i] / 2.0;
                    XM0 = XM0 * XM2 / XM1;
                    BU01 = RHO0 * SQLAM * B[i] * B[i] * TH[i] * Pi * W[i] / (1.0 + SQLAM);
                    BU02 = RHO0 * lambda * B2 * TH2 * sqrt(Pi) * W2 * pWidth / (sqrt(1.0+lm2));
                    BU0 = BU02 / BU01 * BU0;
                    B[i] = B2;
                    W[i] = W2;
                    TH[i] = TH2;
                }
            } else {
                if (Q0 != 0.0)
                    TH[i] = BU0 / Q0;
                else
                    TH[i] = BU0 / Y[i];

                if (i==level)
                //#    ENT[i] = Pi * Q0;
                    ENT[i] = E;  //# 24/02/20
                else
                //#    ENT[i] = Pi * DQ;
                    ENT[i] = E*DZ; // # 24/02/20
            }

            B0 = B[i];
            W0 = W[i];
            TH0 = TH[i];
            DW0DZ = DW/DZ;

        }

        if ( i < NSY ) {
            //# Here when plume has reached top
            //# set iTop(pCount) to level before last one processed
            //# this is because solution is inaccurate as W --> 0
            //# height of rise of initial plume
            iTop[pCount] = i-1;

            if (pCount == 1)
                htRise = yDepth[i-1];

        } else {
            //# Top of water column met
            iTop[pCount] = NSY;
            if (pCount == 1) htRise = yDepth[NSY];
        }

        //# plume has stopped at layer iTop(pCount)
        //# total discharge is Pi*Q0, with ENT(j) from
        //# each layer between level and iTop(pCount).
        //# Q0 includes entrainment from iTop
        //# NOTE that ENT(j) comes from each port; need nPorts*ENT(j) from layer j

        //# Check to see if iTop(pCount) equals surfLayer (IE. top of water column)
        //# If not set level to iTop(pCount) + 1
        if (iTop[pCount] < NSY) {
          //  level = iTop[pCount] + 1;  // 16/12/19
          //  pCount++;                  // 16/12/19 disabled 
            *pCountI = pCount;
            break;         //    16/12/19 test the cycle
        } else
            break;
    } while (1);
    /********************%%%%%%%%%%%%%%%%%%%%%%%%%%%%**********************/
    *retI = i;
    return TRUE;
}


/******************************************************************************
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void inBubble(LOGICAL intopt, int MaxNSy)
{
//CONSTANTS
    const AED_REAL US=0.3E0, HPRES=10.2E0;

    int level, iLevel, i, ii, j, jj, iRiver, ll, IUP;
    int NINT, NINT0, pCount, iSkip;
    int NSY, indexY;

    AED_REAL BSL;
    AED_REAL FLOW;
    AED_REAL RHO0;
    AED_REAL QVOLT, QVOL;

#ifndef _VISUAL_C_
    int iTop[MaxLayers];
    AED_REAL T[MaxLayers], S[MaxLayers], QF[MaxLayers];
    AED_REAL PLMDEN[MaxLayers], BFSQ[MaxLayers];
    AED_REAL yDepth[MaxNSy], yDepthM[MaxNSy], YVOL[MaxNSy];
    AED_REAL ENT[MaxNSy], W[MaxNSy], B[MaxNSy], TH[MaxNSy];
    AED_REAL yTemp[MaxNSy], ySal[MaxNSy], yDen[MaxNSy];
    AED_REAL ywqual[Num_WQ_Vars * MaxNSy];
    AED_REAL xw[Num_WQ_Vars * MaxLayers];
    AED_REAL wqx[Num_WQ_Vars];
#else
    int *iTop;
    AED_REAL *T, *S, *QF;
    AED_REAL *PLMDEN, *BFSQ;
    AED_REAL *yDepth, *yDepthM, *YVOL;
    AED_REAL *ENT, *W, *B, *TH;
    AED_REAL *yTemp, *ySal, *yDen;
    AED_REAL *xw, *wqx, *ywqual;
#endif

#define _ywqual(i, j) ywqual[_IDX_2d(Num_WQ_Vars, MaxNSy, i, j)]
#define _xw(i, j)         xw[_IDX_2d(Num_WQ_Vars, MaxLayers, i, j)]

    AED_REAL BD, XMom0, HL, HU, DENTOP, DENBOT;

    AED_REAL tWidth;
    LOGICAL  done;

/*----------------------------------------------------------------------------*/

    //# Calculated parameters and initialization
    IUP = 0;
    HU = 0.;
    HL = 0.;

#ifdef _VISUAL_C_
    ALLOCATE1D_I(iTop, MaxLayers);

    ALLOCATE1D_R(T, MaxLayers);      ALLOCATE1D_R(S, MaxLayers);    ALLOCATE1D_R(QF, MaxLayers);
    ALLOCATE1D_R(PLMDEN, MaxLayers); ALLOCATE1D_R(BFSQ, MaxLayers);
    ALLOCATE1D_R(yDepth, MaxNSy);
    ALLOCATE1D_R(yDepthM, MaxNSy); ALLOCATE1D_R(YVOL, MaxNSy); ALLOCATE1D_R(ENT, MaxNSy);
    ALLOCATE1D_R(W, MaxNSy);     ALLOCATE1D_R(B, MaxNSy);    ALLOCATE1D_R(TH, MaxNSy);
    ALLOCATE1D_R(yTemp, MaxNSy); ALLOCATE1D_R(ySal, MaxNSy); ALLOCATE1D_R(yDen, MaxNSy);
    ALLOCATE1D_R(wqx,  Num_WQ_Vars);

    ALLOCATE2D(ywqual, Num_WQ_Vars, MaxNSy);
    ALLOCATE2D(xw,     Num_WQ_Vars, MaxLayers);
#endif

    htRise = 0.;
    cee = 0.;
    DENTOP = 0.;
    DENBOT = 0.;
    lambda = 0.3;
    alp2 = 0.102;
    iSkip = 0;
    if (Lake[surfLayer].Temp == Lake[botmLayer].Temp) iSkip = 1;
    lm2 = lambda*lambda;
    pWidth = bubLen/(nPorts);

    BD = bDepth;
    FLOW = aFlow/(nPorts);

    //# Basic bubbler time step = 15 Min
    NINT0 = noSecs;
    NINT = 900;

    //# Initialize sublayer variables
    for (i = 0; i < MaxNSy; i++) {
        yDepth[i] = 0.0;
        yDepthM[i] = 0.0;
        YVOL[i] = 0.0;
        yTemp[i] = 0.0;
        ySal[i] = 0.0;
        yDen[i] = 0.0;
    }

    ub2 = US * (1.0+lm2);
    ub2d = US * sqrt(1.0+lm2);

    do  {
        //# Find center of buoyancy; XMoment1 is 1st moment, XMom0 is 0th moment
        //# use this to estimate two-layer equivalent structure for purposes
        //# of calculating cee
        if (iSkip != 1) {
            for (i = 0; i < MaxLayers; i++) BFSQ[i] = 0.;

            XMoment1 = 0.;
            XMom0 = 0.;

            for (i = 1; i <= surfLayer; i++) {
                BFSQ[i] = gprime(Lake[i].Density, Lake[i-1].Density)/(Lake[i].MeanHeight-Lake[i-1].MeanHeight);
                XMoment1 += Lake[i-1].Height * BFSQ[i] * (Lake[i].LayerVol + Lake[i-1].LayerVol) / 2.;
                XMom0 +=                       BFSQ[i] * (Lake[i].LayerVol + Lake[i-1].LayerVol) / 2.;
            }

            //# Define length scales, XMoment1 is the center of buoyancy ( M above bot
            //#                       HU is the thickness of the upper mixed layer
            if (XMom0 != 0.) {
                HL = XMoment1 / XMom0;
                HU = Lake[surfLayer].Height - HL;
            } else {
                HL = Lake[surfLayer].Height;
                HU = 0.;
            }

            //# Find layer containing the centre of buoyancy (layer i)
            for (i = botmLayer; i <= surfLayer; i++)
                if (Lake[i].Height > HL) break;

            DENBOT = 0.;
            DENTOP = 0.;

            IUP = i;
            if (IUP != surfLayer) {
                for (i = botmLayer; i < IUP; i++)
                    DENBOT += Lake[i].Density*Lake[i].LayerVol;
                DENBOT = DENBOT/Lake[IUP].Vol1;

                for (i = IUP; i <= surfLayer; i++)
                    DENTOP += Lake[i].Density*Lake[i].LayerVol;
                DENTOP = DENTOP/(Lake[surfLayer].Vol1-Lake[IUP].Vol1);
            } else
                cee = 0.0;
        }

        //# Set first value of level
        iLevel = surfLayer;

        for (i = surfLayer; i >= botmLayer; i--)
            if (Lake[i].Height >= bDepth) iLevel = i;

        //# Set conditions at starting level
        RHO0 = Lake[iLevel].Density + 1000.0;

        //# pCount is number of plumes
        pCount = 1;
        level = iLevel;

        //# plume number for 2-layer approximation - assumes bubbler on bottom
        if (IUP == surfLayer)
            cee = 0.;
        else {
            if (DENBOT > DENTOP) {
                cee = pow(((DENBOT-DENTOP)/RHO0)*HL, 1.5)*Lake[surfLayer].Height *
                                      pow(HU/Lake[surfLayer].Height, 3.)/FLOW;
                cee = cee * ((Lake[surfLayer].Height + HPRES) / HPRES) * RHO0 / 1000.;
            } else
                cee = 0.0;
        }

        //# Now get into bubbler internal substeps
        NSY = surfLayer * nStepY;

        for (i = surfLayer; i > iLevel+1; i--) {
            for (j = 1; j <= nStepY; j++) {
                indexY = i*nStepY-j+1;
                yDepth[indexY] = Lake[i].Height - (Lake[i].Height - Lake[i-1].Height) *
                                    (j-1)/(nStepY);
                yDepthM[indexY] = Lake[i].Height - (Lake[i].Height - Lake[i-1].Height) *
                                    ((j-1) + 0.5)/(nStepY);
                YVOL[indexY] = Lake[i].LayerVol/(nStepY);
                yTemp[indexY] = Lake[i].Temp;
                ySal[indexY] = Lake[i].Salinity;
                yDen[indexY] = Lake[i].Density;

                for (ii = 0; ii < Num_WQ_Vars; ii++)
                    _ywqual(ii, indexY) = _WQ_Vars(ii, i);
            }
        }

        for (i = iLevel; i > 1; i--) {
            for (j = 1; j < nStepY; j++) {
                indexY = i*nStepY-j+1;
                if (i != botmLayer) {
                    yDepth[indexY] = Lake[i].Height - (Lake[i].Height - Lake[i-1].Height) *
                                                (j-1)/(nStepY);
                    yDepthM[indexY] = Lake[i].Height - (Lake[i].Height - Lake[i-1].Height) *
                                                ((j-1) + 0.5)/(nStepY);
                } else {
                    yDepth[indexY] = Lake[i].Height-(Lake[i].Height)*
                                             (j-1)/(nStepY);
                    yDepthM[indexY] = Lake[i].Height-(Lake[i].Height)*
                                             ((j-1) + 0.5)/(nStepY);
                }
                YVOL[indexY] = Lake[i].LayerVol/(nStepY);
                yTemp[indexY] = Lake[i].Temp;
                ySal[indexY] = Lake[i].Salinity;
                yDen[indexY] = Lake[i].Density;

                for (ii = 0; ii < Num_WQ_Vars; ii++)
                    _ywqual(ii, indexY) = _WQ_Vars(ii, i);
            }
        }

        level = level * nStepY;

        for (i = 0; i < MaxNSy; i++) {
            ENT[i] = 0.0;
            W[i] = 0.0;
            B[i] = 0.0;
            TH[i] = 0.0;
        }

        done = FALSE;
        do  {
            if ( !done && inner_loop(iLevel, level, &pCount, &i, NSY, iTop, yDepth, yDepthM, TH, W, B, ENT, yDen) ) {
                //# Multiply ENT(j) by the number of ports
                for (j = 0; j < NSY; j++)
                    ENT[j] = (nPorts) * ENT[j];

                //# Check if new time step required
                if (NINT0 < NINT) NINT = NINT0;

                //# Fix entrainments
                for (i = 0; i < NSY; i++) {
                    ENT[i] *= (NINT)/1000.0;    //  18/12/19
                //   ENT[i] *= (NINT)/900.0;    //  18/12/19
                    YVOL[i] -= ENT[i];

                    if (YVOL[i] <= 0.0) {
                        YVOL[i] += ENT[i];
                        if (i != NSY-1)
                            ENT[i+1] += (ENT[i] - 0.9*YVOL[i])*1000./(NINT);  // 18/12/19

                      //      ENT[i+1] += (ENT[i] - 0.9*YVOL[i])*900./(NINT);    // 18/12/19 

                        ENT[i] = 0.9 * YVOL[i];
                        YVOL[i] -= ENT[i];
                    }
                }

                //# Fix properties of plume
                //# -incoming
                i = 1;
                level = 1;
                done = TRUE;
            }

            if ( done ) {
                T[i] = 0.0;
                S[i] = 0.0;
                PLMDEN[i] = 0.0;

                for (ii = 0; ii < Num_WQ_Vars; ii++) _xw(ii, i) = 0.0;

                //# -entrainment
                QF[i] = 0.0;

                for (j = level; j < iTop[i]+1; j++) {   //   j<iTop[i] 18/12/19
                    if (ENT[j] > 0.0) {
                        T[i] = combine(T[i], QF[i], PLMDEN[i], yTemp[j], ENT[j], yDen[j]);
                        S[i] = combine(S[i], QF[i], PLMDEN[i],  ySal[j], ENT[j], yDen[j]);
                        PLMDEN[i] = calculate_density(T[i], S[i]);

                        for (ii = 0; ii < Num_WQ_Vars; ii++)
                            _xw(ii, i) = combine_vol(_xw(ii, i), QF[i], _ywqual(ii, j), ENT[j]);

                        QF[i] += ENT[j];
                    }
                }
                done = FALSE;
            }

            //# Check for next plume
            //# Restart at layer iTop+1 and go to top of next plume
            if (i < pCount) {
                level = iTop[i] + 1;
                i++;
                done = TRUE;
            }
        } while (done);

        //# Now swap back to the regular layer size
        //# Volume is only physical property we need change
          for (ii = botmLayer; ii <= surfLayer; ii++) {
          //   Lake[ii].LayerVol = 0.;
          //   for (jj = 0; jj < nStepY; jj++)                     // 17/12/19
          //      Lake[ii].LayerVol += YVOL[ii*nStepY-jj+1];

            Lake[ii].Temp = yTemp[ii*nStepY];
            Lake[ii].Salinity = ySal[ii*nStepY];
            Lake[ii].Density =yDen[ii*nStepY];

            for (jj = 0; jj < Num_WQ_Vars; jj++)
                _WQ_Vars(jj, ii) = _ywqual(jj, ii*nStepY);
        }

        //# Find the level of intrusion of the water when it falls back
        //# and insert it using the routine from the DoInFlow routine.
        //# assume the basin length at the bubbler is determined by the slope of
        //# the major DoInFlow valley, I.E river 1.
        iRiver = 0;
        BSL = Inflows[iRiver].Phi * Pi / 180.0;

        //# Make insertion WIDTH 0, then insert will use reservoir width at discha
        QVOLT = 0.0;

   //     for (i = 0; i < pCount; i++) {   //  19/12/19 should be pcount +1 
        for (i = 1; i < pCount+1; i++) {   //  20/01/29
            for (jj = 0; jj < Num_WQ_Vars; jj++) wqx[jj] = _xw(jj, i);

            QVOLT = QVOLT+QF[i];
            QVOL = QF[i];

            tWidth = 1.E-8;
            insert(QVOL, PLMDEN[i], BSL, T[i], S[i], wqx, NINT, &tWidth, &ll);
        }

        //# Fix depths with resint
        resize_internals(2, 1);

        //# Fix layers
        check_layer_thickness();

        //# Reinitialize bubbler depth to original (needed by pe3)
        bDepth = BD;

        //# Check if number of seconds left
        NINT0 = NINT0 - NINT;
        if (NINT0 > 0)
            iSkip = 1;
        else
            break;
    } while (1);//# Finished!

#ifdef _VISUAL_C_
    free(iTop);
    free(T);       free(S);      free(QF);
    free(PLMDEN);  free(BFSQ);   free(yDepth);
    free(yDepthM); free(YVOL);   free(ENT);
    free(W);       free(B);      free(TH);
    free(yTemp);   free(ySal);   free(yDen);
    free(wqx);     free(ywqual); free(xw);
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InActRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM)
{
    AED_REAL YT[3], DYT[3], DYM[3], HH, H6;
    int i;

/*----------------------------------------------------------------------------*/

    HH = H*0.5;
    H6 = H/6;

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYDX[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYT[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYM, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++) {
        YT[i] = Y[i] + H * DYM[i];
        DYM[i] = DYT[i] + DYM[i];
    }
    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++)
        YOUT[i] = Y[i] + H6*(DYDX[i] + DYT[i] + 2.*DYM[i]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InActDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM)
{
//CONSTANTS
    const AED_REAL R2 = 1000.;

//LOCALS
    AED_REAL W, B, TH;

/*----------------------------------------------------------------------------*/

    if (fabs(Y[0]) <= 1.E-7) {
        W = W0;
        TH = TH0;
    } else {
        if (DIM==2) {
            W = (sqrt(2.0)) * Y[1] / (RHO0 * Y[0]);
            TH = Y[2] * (sqrt(1.0+lm2)) / (Y[0] * lambda * R2);
        } else {
            W = 2.0 * Y[1] / (RHO0*Y[0]);
            TH = Y[2] * (1.0 + SQLAM) / (SQLAM * Y[0] * R2);
        }
    }

    if (W > 0. && Y[0] > 0.) {
        if (DIM==2)
            B = Y[0] / (W * pWidth * sqrt(Pi));
        else
            B = sqrt(Y[0] / (W * Pi));
    } else
        B = B0;

    if (DIM==2) {
        DYDX[0] = 2.0 * alp2 * pWidth * W;
        DYDX[1] = sqrt(Pi) * RHO0 * TH * lambda * B * pWidth;
        DYDX[2] = sqrt(Pi) * B * pWidth * W * g * DRH * RHO0 +
                  FLOW * PA * (((W-W0)/X) * (1. - W/(W + ub2d)) + W/SMALLZ) /
                  (SMALLZ * (W+ub2d));
    } else {
        DYDX[0] = 2.0 * ALP * B * W * Pi;
        DYDX[1] = Pi * RHO0 * TH * SQLAM * B * B;
        DYDX[2] = Pi * B * B * W * g * DRH * RHO0 +
                  FLOW * PA * (((W-W0)/X) * (1. - W/(W + UB)) + W/SMALLZ) /
                  (SMALLZ * (W+UB));
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InDepRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0,
              AED_REAL E, AED_REAL DENSEA,
              AED_REAL DENSEW, AED_REAL DENSEP, AED_REAL LAMBDA, AED_REAL L,
              AED_REAL TAMB, AED_REAL DOAMB, AED_REAL DNAMB, AED_REAL RB,
              AED_REAL NN, AED_REAL V, AED_REAL VB, AED_REAL KOLO, 
              AED_REAL KOLN, AED_REAL HO2, AED_REAL HN2, AED_REAL PO,
              AED_REAL PN, AED_REAL COMGP, AED_REAL CNMGP)
{
    AED_REAL YT[8], DYT[8], DYM[8], HH, H6;   //# 3 before 28/04/20

    int i;

/*----------------------------------------------------------------------------*/

    HH = H/2;
    H6 = H/6;

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYDX[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, E, DENSEA, DENSEW, DENSEP,LAMBDA,L, TAMB, DOAMB, DNAMB, RB, NN, V, VB, KOLO,  KOLN, HO2, HN2, PO, PN, COMGP, CNMGP);

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH*DYT[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(HH, YT, DYM, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, E, DENSEA, DENSEW, DENSEP, LAMBDA,L, TAMB, DOAMB, DNAMB, RB, NN, V, VB, KOLO,  KOLN, HO2, HN2, PO, PN, COMGP, CNMGP);

    for (i = 0; i < N; i++) {
        YT[i] = Y[i] + H*DYM[i];
        DYM[i] = DYT[i] + DYM[i];
    }
    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(H, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, E, DENSEA, DENSEW, DENSEP, LAMBDA, L, TAMB, DOAMB, DNAMB, RB, NN, V, VB, KOLO,  KOLN, HO2, HN2, PO, PN, COMGP, CNMGP);

    for (i = 0; i < N; i++)
        YOUT[i] = Y[i] + H6*(DYDX[i] + DYT[i] + 2.*DYM[i]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//#void InDepDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
//#              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
//#              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
//#              AED_REAL RHO0, AED_REAL TH0)
void InDepDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, AED_REAL E, AED_REAL DENSEA,
              AED_REAL DENSEW, AED_REAL DENSEP, AED_REAL LAMBDA, AED_REAL L,
              AED_REAL TAMB, AED_REAL DOAMB, AED_REAL DNAMB, AED_REAL RB,
              AED_REAL NN, AED_REAL V, AED_REAL VB, AED_REAL KOLO, 
              AED_REAL KOLN, AED_REAL HO2, AED_REAL HN2, AED_REAL PO,
              AED_REAL PN, AED_REAL COMGP, AED_REAL CNMGP)

//# REAL N is changed to NN // 28/04/20

{
    AED_REAL W, B, TH;
 //# AED_REAL 
 //# 24/02/20
/*----------------------------------------------------------------------------*/

    if (fabs(Y[0]) <= 1.E-7) {
        W = W0;
        TH = TH0;
    } else {
        W = Y[1]/Y[0];
        TH = Y[2]/Y[0];
    }

    if (W > 0. && Y[0] > 0.)
        B = sqrt(Y[0] / W);
    else
        B = B0;
       
        B = B0;    //# 11/03/20  

    //#    B = BW; //# 11/03/20

    DYDX[0] = 2.*ALP*B*W;
    DYDX[1] = 2.*9.81*SQLAM*SQR(B)*TH;
 //#   DYDX[2] = (1.+SQLAM)*DRH*Y[0]/SQLAM
 //#            +((1.+SQLAM)*FLOW*PA/(g*Pi*SQLAM*RHO0*SMALLZ*(W+UB)))*(W/SMALLZ +
   //#                                                   ((W-W0)/X)*(1.-W/(W+UB)));
//#--------------------------------
      DYDX[0]= E;     //# * E,B,L,V is updated; E = 2*ALPHA*(B+2L)*V 
      DYDX[1]= (DENSEA-DENSEW)/DENSEP*9.8*L*2.*B+(DENSEW-DENSEP)/
     +DENSEP*9.8*LAMBDA*2.*B*(L-2.*B*(1.-LAMBDA));
      DYDX[2]=E*TAMB;
      DYDX[3]=0;
      DYDX[4]=(E*DOAMB/32.+4*3.14*RB*pow(RB,2)*NN/(V+VB)*KOLO*(HO2*PO-COMGP/32.));
      DYDX[5]=(E*DNAMB/28.+4*3.14*RB*pow(RB,2)*NN/(V+VB)*KOLN*(HN2*PN-CNMGP/32.)); 
      DYDX[6]= -4*3.14*pow(RB,2)*NN/(V+VB)*KOLO*(HO2*PO-COMGP/32.);
      DYDX[7]= -4*3.14*pow(RB,2)*NN/(V+VB)*KOLN*(HN2*PN-CNMGP/28.);
///#---------------------------------
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
