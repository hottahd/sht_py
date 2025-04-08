import numpy as np
import scipy.special as sp

import sht

def define_geometry(ith, jph):
    """
    Define the geometry of the grid in spherical coordinates.
    Parameters
    ----------
    ith: int
        Number of grid points in the colatitude direction.
    jph: int
        Number of grid points in the longitude direction.
        
    Returns
    -------
    th: ndarray
        Colatitude grid points.
    ph: ndarray
        Longitude grid points.
    """
    dth = np.pi/ith
    dph = 2*np.pi/jph

    th = np.zeros(ith)
    th[0] = 0.5*dth
    for i in range(1, ith):
        th[i] = th[i-1] + dth

    ph = np.zeros(jph)
    ph[0] = 0.5*dph
    for j in range(1, jph):
        ph[j] = ph[j-1] + dph

    return th, ph

def test_simple_forward():
    
    th, ph = define_geometry(ith=128, jph=512)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    
    l = 20 # degree of harmonic
    m = 0 # order of harmnic
    qq = sp.sph_harm_y(l, m, TH, PH).real
    fqq = sht.sht(qq, th, ph)
    
    if m == 0:
        sc = 1.0
    else:
        sc = 0.5
    print(abs(abs(fqq)[l,m]*np.sqrt(2*np.pi) - sc))
    assert abs(abs(fqq)[l,m]*np.sqrt(2*np.pi) - sc) < 1.e-2

def test_forward_backward():
    
    th, ph = define_geometry(ith=128, jph=256)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    
    l = 10 # degree of harmonic
    m = 7 # order of harmnic
    
    qq = sp.sph_harm_y(l, m, TH, PH).real
    fqq = sht.sht(qq, th, ph, direction = 1)
    qqb = sht.sht(fqq, th, ph, direction = -1)
    
    assert np.max(np.abs(qq-qqb)) < 1.e-10
    
    
    