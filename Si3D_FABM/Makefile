all:
	ifort  -module $(MODDIR) -I$(FABMINC) -I$(INCDIR) nspcg.f si3d_types.f90 si3d_boundaryconditions.f90 si3d_mixing_sod.f90 PlumeModels.f si3d_ecomod.f90 si3d_tracer.f90 si3d_fabm_v3_nonuniforminit.f90 si3d_procedures_waterjet_v7_nonuniforminit.f90 si3d.f90 -c -O2
	ifort  -o si3d *.o -module $(MODDIR) -I$(INCDIR) -L$(LIBDIR) -lturbulence_prod -lutil_prod -L$(FABMLIB) -lfabm -O2
si3d:	ifort  nspcg.f PlumeModels.f si3d_types.f90 si3d_ecomod.f90 si3d_boundaryconditions.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -o si3d -O2
