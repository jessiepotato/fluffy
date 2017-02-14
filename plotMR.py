import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Box1DKernel

data=np.genfromtxt('output.txt', names=['P_z', 'T_z', 'Q_c', 'Q_v', 'Qt', 'Qs'])
P = data['P_z']
qc = data['Q_c']
qv = data['Q_v']
qs = data['Qs']
T = data['T_z']
qt = data['Qt']
#K = data['K']
#L = data['L']
#H = data['H']
#tau = data['Tau']
#tau_a = data['Tau_a']
#dq = data['Dqt']

fig = plt.figure()

ax = fig.add_subplot(111)

ax.invert_yaxis()

#ax.plot(Qc,P,'.')

#ax.plot(T,P,'r.')

#qc_smooth = convolve(qc, Box1DKernel(10))
#qc_0 = qc_smooth==0
#qc_smooth[qc_0]=10E-10

#qc_0 = qc==0
#qc[qc_0]=10E-10

#qv_0 = qv==0
#qv[qv_0]=10E-10

fontsize=25
ax.semilogy(np.log10(qs),P,'c-', lw=3, label='Saturation')
ax.semilogy(np.log10(qv),P,'k--', lw=3, label='Vapour')
ax.semilogy(np.log10(qc),P,'r.', lw=3, label='Condensate')
#ax.semilogy((qc),P,'r-', lw=3, label='Condensate')
ax.semilogy(np.log10(qt),P,'m--', lw=3, label='Total')
#ax.semilogy((dq),P,'g.-', label='change in total')
#ax.semilogy(np.log10(qc+qv),P,'g.-', label='Qt calc')
#ax.semilogy((qs),P,'r.-', label='Qs')
#ax.semilogy((qv),P,'k.-', label='Qv')
#ax.semilogy((qc_smooth),P,'b.-', label='Qc')
#ax.semilogy((qt),P,'m.-', label='Qt')
#ax.semilogy((qc2),P,'g.-', label='Qc calc')
#ax.semilogy(np.log10(qc_smooth),P,'b.-')
#ax.semilogy(np.log10(qt),P,'r.')
#ax.semilogy(T,P,'m.-')
#ax.semilogy(H,P,'k-.')
#ax.semilogy(L,P,'b-.')
#ax.semilogy(tau,P,'k-.')
#ax.semilogy(tau_a,P,'b-.')
#ax.semilogy(qt,P,'r.-')
#ax.semilogy(np.log10(qv),P,'k.-')
#ax.semilogy(eta,P,'k.')
#ax.semilogy(30*(abs(dqt)),P,'k--')
ax.legend()
ax.set_xlim([-8, -4])
ax.set_ylim([1, .1])
#ax.set_xlabel('log10(Volume mixing ratio)', fontsize=fontsize)
#ax.set_ylabel('Pressure [bar]', fontsize=fontsize)

#ax2 = ax.twiny()

#ax.plot(np.log10(qc_smooth),P,'k--')

#ax2.set_xlim([110, 160])
fig.savefig('mixingR.png',bbox_inches='tight')

plt.show()
