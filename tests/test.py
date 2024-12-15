import numpy as np
import sht


ith = 128
jph = 256

th = np.linspace(0,  np.pi, ith)
ph = np.linspace(0,2*np.pi, jph)

TH, PH = np.meshgrid(th, ph, indexing='ij')

qq = np.zeros((ith, jph))
for i in range(ith//2):
    for j in range(jph//2):
        amp = np.random.uniform(-1,1)/ith/jph
        qq += amp*np.sin(i*TH)*np.sin(j*PH)

fqq = sht.sht_py(qq, th, ph)
qqb = sht.sht_py(fqq, th, ph,direction=-1)

print((qq - qqb).mean())