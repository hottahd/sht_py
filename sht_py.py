#
# Spherical harmonics transformation
#
# INPUT:
# qq: variable in real space [jx,kx]
# y: colatitude from 0. to pi [jx]
# z: longitude from -pi to pi [kx]
#
# OUTPUT:
# ffqq: complex [kx/4+1,kx/4+1]
import numpy as np
from numpy.fft import fftn, ifftn


def sht(qq,y,z):
    jx, kx = y.shape[0], z.shape[0]
    dy, dz = y[1] - y[0], z[1] - z[0]
    lmax = kx//4
  
    fqq0 = fftn(qq,axes=[1],norm='forward')
    fqq = np.zeros((jx,kx//2+1),dtype=np.complex64)
    fqq[:,0] = fqq0[:,0]*2.*np.pi/float(kx)
    fqq[:,1:kx//2] = (fqq0[:,1:kx//2] + np.conj(fqq0[:,kx:kx//2:-1]))*2.*np.pi/float(kx)
    fqq[:,kx//2] = fqq0[:,kx//2]*2.*np.pi/float(kx)
    
    ffqq = np.zeros((lmax+1,lmax+1),dtype=np.complex64)

    for m in range(0,lmax+1):
        fm = float(m)
        if m == 0:
            pm = 1.e0/np.sqrt(4.e0*np.pi)
        else:
            pm = -np.sqrt((2.*fm+1.)/(2.*fm))*np.sin(y)*pm
            pm[np.where(np.abs(pm) < 1.e-100)] = 0.e0
            
        # Integration
        ffqq[m,m] = np.sum(fqq[:,m]*pm*np.sin(y)*dy)

        pm1 = pm
        pm2 = 0.e0
        for l in range(m+1,lmax+1):
            fl = float(l)

            faca = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
            facb = np.sqrt( \
                    (2.*fl+1.)/(2.*fl-3.)* \
                    (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
            )
            pm0 = faca*np.cos(y)*pm1 - facb*pm2

            # Integration
            ffqq[l,m] = np.sum(fqq[:,m]*pm0*np.sin(y)*dy)
            
            pm2 = pm1
            pm1 = pm0
            
    return(ffqq)

def shtm_restrict(qq,y,z,N,direction=1):
    jx, kx = y.shape[0], z.shape[0]
    dy, dz = y[1] - y[0], z[1] - z[0]

    if direction !=1 and direction != -1:
        "Please use direction = 1 or -1"
        return

    # Original transformation
    if direction == 1:
        # Fourier transformation
        fqq = fftn(qq,axes=[1],norm='forward')
        ffqq = np.zeros((N*jx,kx),dtype=np.complex64)

        # Legendre transformation
        for m in range(0,N*kx//2):
            fm = float(m)
            # l = m case
            if m == 0:
                pm = 1.e0/np.sqrt(4.e0*np.pi)
            else:
                pm = -np.sqrt((2.*fm+1.)/(2.*fm))*np.sin(y)*pm
                pm[np.where(np.abs(pm) < 1.e-100)] = 0.e0
            # negative m
            pn = (-1)**fm*pm
            
            # Integration
            if m % N == 0:
                ffqq[m,m//N] = np.sum(fqq[:,m//N]*pm*np.sin(y)*dy)
            
                if m != 0:
                    ffqq[m,kx-m//N] = np.sum(fqq[:,kx-m//N]*pn*np.sin(y)*dy)

                pm1 = pm
                pm2 = 0.e0
                for l in range(m+1,N*kx//2):
                    fl = float(l)

                    faca = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
                    facb = np.sqrt( \
                                (2.*fl+1.)/(2.*fl-3.)* \
                                (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
                    )
                    pm0 = faca*np.cos(y)*pm1 - facb*pm2
                    # negative m
                    pn0 = (-1)**fm*pm0
                
                    # Integration
                    ffqq[l,m//N] = np.sum(fqq[:,m//N]*pm0*np.sin(y)*dy)
                    if m != 0:
                        ffqq[l,kx-m//N] = np.sum(fqq[:,kx-m//N]*pn0*np.sin(y)*dy)
                
                    pm2 = pm1
                    pm1 = pm0

    # Inverse transformation
    if direction == -1:
        fqq = np.zeros((jx,kx),dtype=np.complex64)

        # Inverse Legendre transformation
        for m in range(0,N*kx//2):
            fm = float(m)
            # l = m case
            if m == 0:
                pm = 1.e0/np.sqrt(4.e0*np.pi)
            else:
                pm = -np.sqrt((2.*fm+1.)/(2.*fm))*np.sin(y)*pm
                pm[np.where(np.abs(pm) < 1.e-35)] = 0.e0
            # negative m
            pn = (-1)**fm*pm
            
            if m % N == 0:
                # Integration
                fqq[:,m//N] = fqq[:,m//N] + qq[m,m//N]*pm

                if m != 0:
                    fqq[:,kx-m//N] = fqq[:,kx-m//N] + qq[m,kx-m//N]*pn

                pm1 = pm
                pm2 = 0.e0
                for l in range(m+1,N*kx//2):
                    fl = float(l)

                    faca = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
                    facb = np.sqrt( \
                                (2.*fl+1.)/(2.*fl-3.)* \
                                (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
                                )
                    pm0 = faca*np.cos(y)*pm1 - facb*pm2
                    # negative m
                    pn0 = (-1)**fm*pm0
                
                    # Integration
                    fqq[:,m//N] = fqq[:,m//N] + qq[l,m//N]*pm0
                    if m != 0:
                        fqq[:,kx-m//N] = fqq[:,kx-m//N] + qq[l,kx-m//N]*pn0
                
                    pm2 = pm1
                    pm1 = pm0
        ffqq = ifftn(fqq,axes=[1],norm='forward')*2*np.pi
                
    return(ffqq)

def shtm(qq,y,z,direction=1):
    jx, kx = y.shape[0], z.shape[0]
    dy, dz = y[1] - y[0], z[1] - z[0]

    if direction !=1 and direction != -1:
        "Please use direction = 1 or -1"
        return

    # Original transformation
    if direction == 1:
        # Fourier transformation
        fqq = fftn(qq,axes=[1],norm='forward')
        ffqq = np.zeros((jx,kx),dtype=np.complex64)

        # Legendre transformation
        for m in range(0,kx//2):
            fm = float(m)
            # l = m case
            if m == 0:
                pm = 1.e0/np.sqrt(4.e0*np.pi)
            else:
                pm = -np.sqrt((2.*fm+1.)/(2.*fm))*np.sin(y)*pm
                pm[np.where(np.abs(pm) < 1.e-100)] = 0.e0
            # negative m
            pn = (-1)**fm*pm
            
            # Integration
            ffqq[m,m] = np.sum(fqq[:,m]*pm*np.sin(y)*dy)
            
            if m != 0:
                ffqq[m,kx-m] = np.sum(fqq[:,kx-m]*pn*np.sin(y)*dy)

            pm1 = pm
            pm2 = 0.e0
            for l in range(m+1,kx//2):
                fl = float(l)

                faca = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
                facb = np.sqrt( \
                            (2.*fl+1.)/(2.*fl-3.)* \
                            (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
                )
                pm0 = faca*np.cos(y)*pm1 - facb*pm2
                # negative m
                pn0 = (-1)**fm*pm0
                
                # Integration
                ffqq[l,m] = np.sum(fqq[:,m]*pm0*np.sin(y)*dy)
                if m != 0:
                    ffqq[l,kx-m] = np.sum(fqq[:,kx-m]*pn0*np.sin(y)*dy)
                
                pm2 = pm1
                pm1 = pm0

    # Inverse transformation
    if direction == -1:
        fqq = np.zeros((jx,kx),dtype=np.complex64)

        # Inverse Legendre transformation
        for m in range(0,kx//2):
            fm = float(m)
            # l = m case
            if m == 0:
                pm = 1.e0/np.sqrt(4.e0*np.pi)
            else:
                pm = -np.sqrt((2.*fm+1.)/(2.*fm))*np.sin(y)*pm
                pm[np.where(np.abs(pm) < 1.e-35)] = 0.e0
            # negative m
            pn = (-1)**fm*pm
            
            # Integration
            fqq[:,m] = fqq[:,m] + qq[m,m]*pm

            if m != 0:
                fqq[:,kx-m] = fqq[:,kx-m] + qq[m,kx-m]*pn

            pm1 = pm
            pm2 = 0.e0
            for l in range(m+1,kx//2):
                fl = float(l)

                faca = np.sqrt((2.*fl-1.)*(2.*fl+1.)/(fl+fm)/(fl-fm))
                facb = np.sqrt( \
                            (2.*fl+1.)/(2.*fl-3.)* \
                            (fl + fm - 1.)*(fl - fm - 1.)/(fl + fm)/(fl - fm) \
                )
                pm0 = faca*np.cos(y)*pm1 - facb*pm2
                # negative m
                pn0 = (-1)**fm*pm0
                
                # Integration
                fqq[:,m] = fqq[:,m] + qq[l,m]*pm0
                if m != 0:
                    fqq[:,kx-m] = fqq[:,kx-m] + qq[l,kx-m]*pn0
                
                pm2 = pm1
                pm1 = pm0
        ffqq = ifftn(fqq,axes=[1],norm='forward')*2*np.pi
                
    return(ffqq)

import legendre
def shtf(qq,y,z,direction=1):
    jx, kx = y.shape[0], z.shape[0]

    if direction == 1:
        fqq = fftn(qq,axes=[1],norm='forward')
        ffqq = legendre.forward(fqq,y,z,jx,kx)
    else:
        fqq = legendre.backward(qq,y,z,jx,kx)
        ffqq = ifftn(fqq,axes=[1],norm='forward')*2*np.pi
    return ffqq
