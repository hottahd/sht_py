#
# Spherical harmonics transformation
#
import numpy as np
#from numpy.fft import fftn, ifftn
from scipy.fft import rfftn, irfftn
from . import legendre
import time

def sht(qq,y,z,N=1,direction=1):
    jx, kx = y.shape[0], z.shape[0]
    
    if direction == 1:
        fft_start = time.time()
        fqq = rfftn(qq,axes=[1],norm='forward')
        fft_end = time.time()
        ffqq = legendre.forward(N,fqq,y,jx,kx)
        led_end = time.time()
        print(fft_end - fft_start,led_end - fft_end)
    else:
        fft_start = time.time()
        fqq = legendre.backward(N,qq,y,kx)
        fft_end = time.time()
        ffqq = irfftn(fqq,axes=[1],norm='forward')
        led_end = time.time()
        print(fft_end - fft_start,led_end - fft_end)
    return ffqq
