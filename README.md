# sht_py
Spherical Harmonic Expansion routine in python+fortran


## modules
- sht_py: pure python functions
- sht_f2py: funcstions using f2py (faster)

If you want to use sht_f2py, please use command.
```
f2py --fcompiler=gfortran -m legendre -c --f90flags='-O3 -fopenmp' -lgomp legendre.f90
```

## functions
please see ```help(sht_py)```