.PHONY: all

FC := gfortran -I/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/include -L/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib -O3 -march=native
#FC := ifort -I/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/include -L/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib -O3 -xHost -ipo -no-prec-div -unroll -inline
#FC := ifort -I/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/include -L/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib -lnetcdff -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib -O3 -xHost -ipo -no-prec-div -unroll -inline

all: compute_optdepth compute_vsedi

compute_optdepth: compute_optdepth.f90
	$(FC) -o $@ $^ -lnetcdff

compute_vsedi: compute_vsedi.f90
	$(FC) -o $@ $^ -lnetcdff
