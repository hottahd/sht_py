# sht_py
Spherical Harmonic Expansion routine in python+fortran

## Files
- sht.sht: functions using ctype (faster and parallerized)
- sht.sht_py: pure python functions

## Install
### Requirement
You need `make` and`gfortran`.

Debian-based (Debian, Ubuntu, Mint, etc…)

```shell
sudo apt update
sudo apt-get install -y build-essential gfortran
```

```shell
pip install .
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
please see ```help(sht.sht_py)``` and ```help(sht.sht)```