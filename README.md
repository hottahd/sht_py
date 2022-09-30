# sht_py
Spherical Harmonic Expansion routine in python+fortran


## modules
- sht_py: pure python functions
- sht_ctype: functions using ctype (faster and parallerized)

## Install
### sht_py
You just need ```numpy``` and ```scipy```.

### sht_ctype
please just make at sht/ directory
```
make
```
You might need to set enviroment variable and increase stack size
```
export OMP_STACKSIZE=512000
ulimit -s unlimited
```

The number of threads is defined as
```
export OMP_NUM_THREADS=4

```

For MacOS X, it does not work and we cannot use large size of quantities for the transformation.

## functions
please see ```help(sht_py)``` and ```help(sht_ctype)```