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
    direction: int
       1: Forward transform
      -1: Backward transform
    Return
    ----------
    ffqq: complex [N*jx,kx] for forward and float [jx,kx] for backward
        a transformed variable
    """

    os.environ['OMP_STACKSIZE'] = str(512000)
    
    jx, kx = y.shape[0], z.shape[0]
    libdir = os.path.dirname(os.path.realpath(__file__))
    lib = np.ctypeslib.load_library('legendre_associated_function.so',libdir)

    N_C = ctypes.byref(ctypes.c_int64(N))
    jx_C = ctypes.byref(ctypes.c_int64(jx))
    kx_C = ctypes.byref(ctypes.c_int64(kx))
        
    if direction == 1:
        lib.forward.argtypes = [
            ctypes.POINTER(ctypes.c_int64), # N
            np.ctypeslib.ndpointer(dtype=np.complex128),  # qqg
            np.ctypeslib.ndpointer(dtype=np.float64),  # yg
            ctypes.POINTER(ctypes.c_int64),  # jxg
            ctypes.POINTER(ctypes.c_int64),  # kx
            np.ctypeslib.ndpointer(dtype=np.complex128)  # fqqg
        ]
        lib.forward.restype = ctypes.c_void_p
           
        fqq = rfftn(qq,axes=[1],norm='forward').astype(np.complex128).reshape(kx//2+1,jx,order='F')
        ffqq = np.zeros((kx//2+1,N*kx//2+1),dtype=np.complex128)
            
        lib.forward(N_C,fqq,y,jx_C,kx_C,ffqq)
    else:
        lib.backward.argtypes = [
            ctypes.POINTER(ctypes.c_int64), # N
            np.ctypeslib.ndpointer(dtype=np.complex128),  # qq
            np.ctypeslib.ndpointer(dtype=np.float64),  # yg
            ctypes.POINTER(ctypes.c_int64),  # jxg
            ctypes.POINTER(ctypes.c_int64),  # kx
            np.ctypeslib.ndpointer(dtype=np.complex128)  # fqqg
        ]
        lib.backward.restype = ctypes.c_void_p        
            
        fqq = np.zeros((kx//2+1,jx),dtype=np.complex128)
        lib.backward(N_C,qq,y,jx_C,kx_C,fqq)

        ffqq = irfftn(fqq,axes=[0],norm='forward')
                
    return ffqq.T

def legendre_polynomial(y,l):
    """
    calculate Legendre polynomial

    Parameters
    ----------
    y: float [jx]
        colatitude 0<y<pi [radian]
    l: int
        degree of Legendre polynomial
        
    Return
    ----------
    pm: float [jx]
        Legendre polynomial of degree l
    """
    
    jx = y.shape[0]
    libdir = os.path.dirname(os.path.realpath(__file__))
    lib = np.ctypeslib.load_library('legendre_polynomial.so',libdir)

    l_C = ctypes.byref(ctypes.c_int64(l))
    jx_C = ctypes.byref(ctypes.c_int64(jx))

    lib.polynomial.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),  # y
        ctypes.POINTER(ctypes.c_int64), # jx
        ctypes.POINTER(ctypes.c_int64), # l
        np.ctypeslib.ndpointer(dtype=np.float64),  # pm
    ]
    
    lib.polynomial.restype = ctypes.c_void_p
    
    pm = np.zeros((jx),dtype=np.float64)
    lib.polynomial(y,jx_C,l_C,pm)
    
    return pm

def legendre_transform(qq,y,direction=1):
    """
    Legendre polynoial transform

    Parameters
    ----------
    qq: float [jx]
        a variable for transform
    y: float [jx]
        colatitude 0<y<pi [radian]
    direction: int
         1: Forward transform
        -1: Backward transform
    
    Return
    ----------
    """
    
    jx = y.shape[0]
    libdir = os.path.dirname(os.path.realpath(__file__))
    lib = np.ctypeslib.load_library('legendre_polynomial.so',libdir)

    jx_C = ctypes.byref(ctypes.c_int64(jx))

    if direction == 1:
        lib.forward.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),  # qq
        np.ctypeslib.ndpointer(dtype=np.float64),  # y
        ctypes.POINTER(ctypes.c_int64), # jx
        np.ctypeslib.ndpointer(dtype=np.float64),  # fqq
        ]
        lib.forward.restype = ctypes.c_void_p    

        fqq = np.zeros((jx),dtype=np.float64)
        lib.forward(qq,y,jx_C,fqq)
    else:
        lib.backward.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),  # fqq
            np.ctypeslib.ndpointer(dtype=np.float64),  # y
            ctypes.POINTER(ctypes.c_int64), # jx
            np.ctypeslib.ndpointer(dtype=np.float64),  # qq            
        ]
        lib.backward.restype = ctypes.c_void_p
        
        fqq = np.zeros((jx),dtype=np.float64)
        lib.backward(qq,y,jx_C,fqq)
    
    return fqq
