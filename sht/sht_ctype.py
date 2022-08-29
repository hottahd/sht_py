from scipy.fft import rfftn, irfftn
import ctypes
import numpy as np
import os

def sht(qq,y,z,N=1,direction=1):
    """
    Spherical Harmonic transform

    Parameters
    ----------
    qq: float [jx,kx] for forward [N*jx,kx] for backward
        a variable for transform
    y: float [jx]
        colatitude 0<y<pi [radian]
    z: float [kx]
        longitude 0<z<2pi/N [radian]
    N: int
       ratio of the longitude extent to 2pi

    Return
    ----------
    ffqq: complex [N*jx,kx] for forward and float [jx,kx] for backward
        a transformed variable
    """
    
    jx, kx = y.shape[0], z.shape[0]
    libdir = os.path.dirname(os.path.realpath(__file__))
    legendre = np.ctypeslib.load_library('legendre.so',libdir)

    N_C = ctypes.byref(ctypes.c_int64(N))
    jx_C = ctypes.byref(ctypes.c_int64(jx))
    kx_C = ctypes.byref(ctypes.c_int64(kx))
        
    if direction == 1:
        legendre.forward.argtypes = [
            ctypes.POINTER(ctypes.c_int64), # N
            np.ctypeslib.ndpointer(dtype=np.complex128),  # qqg
            np.ctypeslib.ndpointer(dtype=np.float64),  # yg
            ctypes.POINTER(ctypes.c_int64),  # jxg
            ctypes.POINTER(ctypes.c_int64),  # kx
            np.ctypeslib.ndpointer(dtype=np.complex128)  # fqqg
        ]
        legendre.forward.restype = ctypes.c_void_p
           
        fqq = rfftn(qq,axes=[1],norm='forward').reshape(kx//2+1,jx,order='F')
        ffqq = np.zeros((kx//2+1,N*kx//2+1),dtype=np.complex128)
            
        legendre.forward(N_C,fqq,y,jx_C,kx_C,ffqq)
    else:
        legendre.backward.argtypes = [
            ctypes.POINTER(ctypes.c_int64), # N
            np.ctypeslib.ndpointer(dtype=np.complex128),  # qq
            np.ctypeslib.ndpointer(dtype=np.float64),  # yg
            ctypes.POINTER(ctypes.c_int64),  # jxg
            ctypes.POINTER(ctypes.c_int64),  # kx
            np.ctypeslib.ndpointer(dtype=np.complex128)  # fqqg
        ]
        legendre.backward.restype = ctypes.c_void_p        
            
        fqq = np.zeros((kx//2+1,jx),dtype=np.complex128)
        legendre.backward(N_C,qq,y,jx_C,kx_C,fqq)

        ffqq = irfftn(fqq,axes=[0],norm='forward')
        
    return ffqq.T
