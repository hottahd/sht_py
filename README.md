# sht_py
Spherical Harmonic Expansion routine in python+fortran


## modules
- sht_py: pure python functions
- sht_f2py: functions using f2py (faster)
- sht_ctype: functions using ctype (fastest)

## Install
### sht_py
### sht_f2py
please use a command at sht/ directory
```
f2py --fcompiler=gfortran -m legendre -c --f90flags='-O3 -fopenmp' -lgomp legendre.f90
```
### sht_ctype
please just make at sht/ directory
```
make
```
You might need to set enviroment variable and increase stack size
```
export OMP_STACKSIZE=512000
export OMP_NUM_THREADS=4
ulimit -s unlimited
```
For Mac OS X, it does not work and we cannot use large size of quantities for the transformation.

## functions
please see ```help(sht_py)```