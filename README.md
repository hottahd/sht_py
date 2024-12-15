# sht_py
Spherical Harmonic Expansion routine in python+fortran

## Files
- sht/sht_py: pure python functions
- sht/sht_ctype: functions using ctype (faster and parallerized)

## Install
### Requirement
You need `make` and`gfortran`.

Debian-based (Debian, Ubuntu, Mint, etc…)

```shell
sudo apt update
sudo apt-get install -y build-essential gfortran
```

You just need ```numpy``` and ```scipy```.

### sht_ctype
please just make at sht/ directory
```
make
```

## OpenMP Usage

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