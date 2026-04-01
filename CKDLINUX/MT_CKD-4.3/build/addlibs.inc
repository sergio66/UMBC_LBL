# here the following variables must be defined
# NCL netCDF lib location
# NCI netCDF include location
# ADDLIB any additional libs required by netCDF lib to compile
ADDLIB = 
NETCDF = no

ifeq ($(FC),gfortran)
	NCL = /usr/lib
	NCI = /usr/include
	NCI = /usr/ebuild/installs/software/netCDF-Fortran/4.6.1-iimpi-2023a/include
	INCLUDES:= -I. -I$(NCI) 
	NETCDF = yes
endif

