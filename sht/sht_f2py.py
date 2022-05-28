#
# Spherical harmonics transformation
#
import numpy as np
from numpy.fft import fftn, ifftn
from . import legendre

def sht(qq,y,z,N=1,direction=1):
    jx, kx = y.shape[0], z.shape[0]

    if direction == 1:
        fqq = fftn(qq,axes=[1],norm='forward')
        ffqq = legendre.forward(N,fqq,y,z,jx,kx)
    else:
        fqq = legendre.backward(N,qq,y,z,jx,kx)
        ffqq = ifftn(fqq,axes=[1],norm='forward')
    return ffqq

    
