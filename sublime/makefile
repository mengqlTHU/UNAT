#fc  	=  tau_f90.sh
fc  	=  swgfortran
exe		=  mbt

#flag	=  -O0 -Wall -Wextra -Wline-truncation -Wcharacter-truncation \
		   -Wsurprising -Waliasing -Wunused-parameter -Wconversion\
		   -fwhole-file -fcheck=all -fbacktrace -fPIC -J mod -std=legacy
#flag	=  -O2 -march=native -fwhole-file -fPIC -J mod -std=legacy

#flag	=  -O0 -module mod -check all -warn all -warn unused -fpe3 -traceback \
		   -gen-interface -check noarg_temp_created -fpic
flag	=  -g -fpic -cpp -DSW_slave -DSW_time
flag_mod = 
#libdir	=  -L/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/cgnslib_2.5/INSTALL -lcgns -L/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/ParMetis-3.2.0 -lparmetis -lmetis
libdir	=  -L/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/cgnslib_2.5/LINUX64 -lcgns \
		   -L/home/export/online1/systest/swrh/installer/parmetis-4.0.3/SWInstall_32_64/lib -lparmetis \
		   -L/home/export/online1/systest/swrh/installer/parmetis-4.0.3/SWInstall_32_64/lib -lmetis \
		   -L/home/export/online1/swmore/opensource/swlu/lib -lswlu_mpi \
		   -L./SWComponent/lib -linterface -lUNAT \
		   -lstdc++ \
		   -L/usr/sw-mpp/mpi2/mpiswgcc/lib -lmpi 
incdir	=  -I/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/cgnslib_2.5/INSTALL \
		   -I/home/export/online1/systest/swrh/installer/parmetis-4.0.3/include \
		   -I/home/export/online1/swmore/opensource/swlu/include \
		   -I/use/sw-mpp/mpi2/mpiswgcc/include 
#incdir	=  -I/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/cgnslib_2.5/INSTALL -I/home/export/online1/systest/swrh/lhb/Tsinghua_UNS/ParMetis-3.2.0
inc		=  include/var_def.f90

other	=  setup.o parallel.o sfci.o sfcs.o sfcd.o linear_solver.o prec.o \
		   fluid.o sem.o
main	=  main.o $(other)
mesh	=  mesh.o mesh_before_balance.o mesh_after_balance.o mesh_io.o \
		   mesh_hanging_node.o mesh_structured.o mesh_gmg.o bnd.o mesh_geo.o \
		   mesh_line.o mesh_wall.o mesh_patch.o mesh_ls.o
fr		=
fv		=  fv.o fv_bnd.o fv_wall.o fv_io.o fv_rhs.o fv_slv_exp.o fv_slv_imp.o\
		   fv_entropy.o fv_gmg.o fv_rhs_upw.o fv_rhs_cen.o\
		   fv_rans.o fv_transition.o fv_les.o fv_rans_sa.o fv_rans_wa.o fv_rans_sst.o \
		   fv_slv_monitor.o fv_reot.o fv_uns.o fv_uns_io.o fv_imex.o \
		   fv_slv_gmres.o
fv_adj	=
heat	=
blas	=  blas.o mc60.o
cht		=

object	=  $(main) $(mesh) $(blas) $(fr) $(fv) $(fv_adj) $(heat) $(cht)

edit:		dir_mod $(object)
			$(fc) $(flag) -mhybrid -o $(exe) $(object) $(libdir) $(incdir)

dir_mod:
			mkdir -p mod
$(main):	%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(mesh):	%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(blas):	%.o: %.f
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(fr):		%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(fv):		%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(fv_adj):	%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(heat):	%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
$(cht):		%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)

clean:
			rm -f $(exe) *.o *.a *.so mod/*.* *.mod

f2pyopt		=  --fcompiler=intelem --f90exec=mpif90 --noopt --noarch
py_obj		=  $(var) $(other) $(mesh) $(blas) $(fv) $(heat) $(cht) $(py_face)
py_face		=  py_face.o
libf		=  -L/home/suxr/Desktop/fv -L/usr/local/lib/ -lpumbt -lcgns \
			   -lparmetis -lmetis
$(py_face):	%.o: %.f90
			$(fc) $(flag) $(flag_mod) -c $< -o $@ $(libdir) $(incdir)
py:			dir_mod $(py_obj)
			ar rcs libpumbt.a $(py_obj)
			f2py $(f2pyopt) -m solid -c py_main.f90 $(libf)
			f2py $(f2pyopt) -m fluid -c py_main.f90 $(libf)
