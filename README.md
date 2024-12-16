# sht_py
Spherical Harmonic Expansion routine in python+fortran

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

## Functions
please see ```help(sht.sht_py)``` and ```help(sht.sht)```

## Example

```python
import numpy as np
import sht

th = np.linspace(0,  np.pi, 128)
ph = np.linspace(0,2*np.pi, 256)

TH, PH = np.meshgrid(th, ph, indexing='ij')

qq = np.sin(TH)*np.cos(PH)

# forward expansion
fqq = sht.sht(qq, th, ph)
# backward expansion
qqb = sht.sht(fqq, th, ph, direction = -1)


```