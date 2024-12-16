import numpy as np
import scipy.special as sp
import sht


ith = 512
jph = 1024

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

TH, PH = np.meshgrid(th, ph, indexing='ij')

qq = sp.sph_harm(40, 200, PH, TH).real

# qq = np.zeros((ith, jph))
# for i in range(ith//2):
#     print(i)
#     for j in range(jph//2):
#         amp1 = np.random.uniform(-1,1)/ith/jph
#         amp2 = np.random.uniform(-1,1)/ith/jph
#         amp3 = np.random.uniform(-1,1)/ith/jph
#         amp4 = np.random.uniform(-1,1)/ith/jph
#         qq += (amp1*np.sin(i*TH) + amp2*np.cos(i*TH)) \
#             *(amp3*np.sin(j*PH) + amp4*np.cos(j*PH))
fqq = sht.sht(qq, th, ph)
qqb = sht.sht(fqq, th, ph, direction=-1)

print(np.sqrt(((qq - qqb)**2).mean()))