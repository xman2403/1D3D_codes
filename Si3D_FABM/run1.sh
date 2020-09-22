########################################################################
########################################################################
###	   This file compile SI3D with some subroutines of GOTM

###	Check permissions of run.sh
###     if permission denied --> type 'sudo chmod 755 run.sh'
###     To compile --> './run.sh'

###     MAC

########################################################################
########################################################################

SHELL=$!/bin/sh

export FORTRAN_COMPILER=IFORT
export SI3DDIR=/media/sf_Shared/Si3D_FABM/v12_tracerinit2_WQ_lightmodel/
export GOTMDIR=/home/chris/gotm-4.0.0
export FABMDIR=/home/chris/local/fabm/Si3D
export FABMSRC=/media/sf_Shared/fabm/fabm-code/src
export MODDIR=$GOTMDIR/modules/IFORT
export FABMINC=$FABMDIR/include
export FABMLIB=$FABMDIR/lib
export INCDIR=$GOTMDIR/include
export BINDIR=$GOTMDIR/bin
export LIBDIR=$GOTMDIR/lib/IFORT

#IF IFGOTM=false --> The modules 'util' and 'turbulence' are not compiled
export IFGOTM=false
#IF IFSI3D=false --> Only SI3D is compiled
export IFSI3D=true

if $IFGOTM
then
	cd $GOTMDIR/src/util
	make clean
	make
	cd $GOTMDIR/src/turbulence
	make clean
	make
fi

if $IFSI3D
then
	cd $SI3DDIR
	make all
	rm *.o
	cd $MODDIR
	rm si3d*
else
	cd $SI3DDIR
	make si3d
	rm *.o
fi
