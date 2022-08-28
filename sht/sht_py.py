#
# Spherical harmonics transformation
#
import numpy as np
from scipy.fft import fftn, ifftn

##################################################################
##################################################################

def legendre_init(y,kx,N):
    """
    Initialize associate Legendre function

    Parameters
    ----------
    y: float [jx]
        colatitude [radian]
    kx: int
        number of grid point in longitude
    N: int
        ratio of the longitude extent to 2pi
    
    Return
    ----------
    siny: float [jx]
        sin(y)
    cosy: float [jx]
        cos(y)
    sindy: float [jx]
        sin(y)*dy, where dy is the grid spacing in latitude
    epm: float [N*kx//2]
        factor for recurrent formula for P_m^m
    faca, facb: float [N*kx//2]
        factor for recurrent formula for P_l^m
    """

    jx = len(y)

    epm = np.zeros(N*kx//2)
    faca = np.zeros((N*kx//2,N*kx//2))
    facb = np.zeros((N*kx//2,N*kx//2))

    dy = y[1] - y[0]

    epm[0] = 0

    for m in range(1,N*kx//2):
        fm = float(m)
        epm[m] = -np.sqrt((2.*fm+1.)/(2.*fm))

    for m in range(0,N*kx//2):
        fm = float(m)
        for l in range(m+1,N*kx//2-1):
        
            fl = float(l)
            faca[l,m] = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
            facb[l,m] = np.sqrt( \
                    (2.*fl+1.)/(2.*fl-3.)* \
                    (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
                )

    cosy = np.cos(y)
    siny = np.sin(y)
    sindy = np.sin(y)*dy

    return(siny,cosy,sindy,epm,faca,facb)

##################################################################
##################################################################

def legendre_m_up(m,siny,epm_m,pm):
    """
    Calculate associate Legendre function P_m^m from P_{m-1}^{m-1}

    Parameters
    ----------
    m: int
        order for associate Legendre function for return
    y: float [jx]
        colatitude [radian]
    pm: float [jx]
        associate Legendre function P_{m-1}^{m-1}
    
    Return
    ----------
    pm: float [jx]
        associate Legendre function P_m^m for (m>0)
    pn: float [jx]
        associate Legendre function P_m^m for (m<0)
    """
    
    # l = m case
    pm = epm_m*siny*pm
    pm[np.where(np.abs(pm) < 1.e-300)] = 0.e0
    # negative m
    pn = (-1)**m*pm

    return(pm,pn)

##################################################################
##################################################################

def legendre_l_up(l,m,cosy,faca_lm,facb_lm,pm1,pm2):
    """
    Calculate associate Legendre function P_l^m from P_{l-1}^{m} and P_{l-2}^m

    Parameters
    ----------
    l: int
        order for associate Legendre function for return
    m: int
        order for associate Legendre function for return
    y: float [jx]
        colatitude [radian]
    pm1: float [jx]
        associate Legendre function P_{l-1}^{m}
    pm2: float [jx]
        associate Legendre function P_{l-2}^{m}

    
    Return
    ----------
    pm0: float [jx]
        associate Legendre function P_l^m for (m>0)
    pm1: float [jx]
        associate Legendre function P_l^m for (m<0)
    """

    pm0 = faca_lm*cosy*pm1 - facb_lm*pm2
    pn0 = (-1)**m*pm0

    return(pm0,pn0)
    
##################################################################
##################################################################

def sht(qq,y,z,N=1,direction=1):
    """
    Spherical Harmonic transform

    Parameters
    ----------
    qq: float or complex [jx,kx] for forward [N*jx,kx] for backward
        a variable for transform
    y: float [jx]
        colatitude 0<y<pi [radian]
    z: float [kx]
        longitude 0<z<2pi/N [radian]
    N: int
       ratio of the longitude extent to 2pi

    Return
    ----------
    ffqq: complex [N*jx,kx] for forward and [jx,kx] for backward
        a transformed variable
    """

    jx, kx = y.shape[0], z.shape[0]
    dy, dz = y[1] - y[0], z[1] - z[0]

    if direction !=1 and direction != -1:
        "Please use direction = 1 or -1"
        return

    # Original transformation
    if direction == 1:
        # Fourier transformation
        fqq = fftn(qq,axes=[1],norm='forward')
        ffqq = np.zeros((N*kx//2,kx),dtype=np.complex64)

        siny,cosy,sindy,epm,faca,facb = legendre_init(y,kx,N)

        m = 0
        pm = 1.0

        ffqq[m,m//N] = np.sum(fqq[:,m//N]*pm*sindy)

        pm1 = pm
        pm2 = 0.0

        for l in range(m+1,N*kx//2):
            # P_l^m from P_{l-1}^{m} and P_{l-2}^m
            pm0,pn0 = legendre_l_up(l,m,cosy,faca[l,m],facb[l,m],pm1,pm2)
                
            # Integration
            ffqq[l,m//N] = np.sum(fqq[:,m//N]*pm0*sindy)
                
            pm2 = pm1
            pm1 = pm0

        # Legendre transformation
        for m in range(1,N*kx//2):
            # l = m case
            # from P_{m-1}^{m-1} to P_m^m
            pm, pn = legendre_m_up(m,siny,epm[m],pm)
            
            # Integration
            if m % N == 0:
                ffqq[m,   m//N] = np.sum(fqq[:,   m//N]*pm*sindy)
                ffqq[m,kx-m//N] = np.sum(fqq[:,kx-m//N]*pn*sindy)

                pm1 = pm
                pm2 = 0.e0
                for l in range(m+1,N*kx//2):
                    # P_l^m from P_{l-1}^{m} and P_{l-2}^m
                    pm0,pn0 = legendre_l_up(l,m,cosy,faca[l,m],facb[l,m],pm1,pm2)
                
                    # Integration
                    ffqq[l,   m//N] = np.sum(fqq[:,   m//N]*pm0*sindy)
                    ffqq[l,kx-m//N] = np.sum(fqq[:,kx-m//N]*pn0*sindy)
                
                    pm2 = pm1
                    pm1 = pm0

        # additional factor for the spherical harmonic expansion
        ffqq = 0.5*ffqq
    # Inverse transformation
    if direction == -1:
        fqq = np.zeros((jx,kx),dtype=np.complex64)
        
        siny,cosy,sindy,epm,faca,facb = legendre_init(y,kx,N)

        m = 0
        pm = 1.0
        fqq[:,m//N] = fqq[:,m//N] + qq[m,m//N]*pm

        pm1 = pm
        pm2 = 0.e0
        for l in range(m+1,N*kx//2):
            # P_l^m from P_{l-1}^{m} and P_{l-2}^m
            pm0,pn0 = legendre_l_up(l,m,cosy,faca[l,m],facb[l,m],pm1,pm2)
                
            # Integration
            fqq[:,m//N] = fqq[:,m//N] + qq[l,m//N]*pm0
                
            pm2 = pm1
            pm1 = pm0


        # Inverse Legendre transformation
        for m in range(1,N*kx//2):
            # l = m case
            # from P_{m-1}^{m-1} to P_m^m            
            pm, pn = legendre_m_up(m,siny,epm[m],pm)
            
            if m % N == 0:
                # Integration
                fqq[:,m//N] = fqq[:,m//N] + qq[m,m//N]*pm

                if m != 0:
                    fqq[:,kx-m//N] = fqq[:,kx-m//N] + qq[m,kx-m//N]*pn

                pm1 = pm
                pm2 = 0.e0
                for l in range(m+1,N*kx//2):
                    # P_l^m from P_{l-1}^{m} and P_{l-2}^m
                    pm0,pn0 = legendre_l_up(l,m,cosy,faca[l,m],facb[l,m],pm1,pm2)
                
                    # Integration
                    fqq[:,m//N] = fqq[:,m//N] + qq[l,m//N]*pm0
                    if m != 0:
                        fqq[:,kx-m//N] = fqq[:,kx-m//N] + qq[l,kx-m//N]*pn0
                
                    pm2 = pm1
                    pm1 = pm0
        ffqq = ifftn(fqq,axes=[1],norm='forward')#*2*np.pi
                
    return(ffqq)
