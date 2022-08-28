from scipy.fft import rfftn, irfftn
import ctypes
import numpy as np

def sht(qq,y,z,N=1,direction=1):
    legendre = np.ctypeslib.load_library('legendre.so','./sht')
    legendre.forward.argtypes = [
        ctypes.POINTER(ctypes.c_int64), # N
        np.ctypeslib.ndpointer(dtype=np.complex128),  # qqg
        np.ctypeslib.ndpointer(dtype=np.float64),  # yg
        ctypes.POINTER(ctypes.c_int64),  # jxg
        ctypes.POINTER(ctypes.c_int64),  # kx
        np.ctypeslib.ndpointer(dtype=np.complex128)  # fqqg
    ]
    legendre.forward.restype = ctypes.c_void_p

    jx, kx = y.shape[0], z.shape[0]
    fqq = rfftn(qq,axes=[1],norm='forward').reshape(kx//2+1,jx,order='F')
    ffqq = np.zeros((kx//2+1,N*kx//2+1),dtype=np.complex128)

    N_C = ctypes.byref(ctypes.c_int64(N))
    jx_C = ctypes.byref(ctypes.c_int64(jx))
    kx_C = ctypes.byref(ctypes.c_int64(kx))

    legendre.forward(N_C,fqq,y,jx_C,kx_C,ffqq)

    return ffqq.T
