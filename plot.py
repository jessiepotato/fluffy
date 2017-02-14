import matplotlib.pyplot as plt
import numpy as np

data=np.genfromtxt('output.txt', names=['P_z', 'x'])
P = data['P_z']
x = data['x']

P1 = [0.432299642, 0.430553093, 0.428729946, 0.42547262, 0.423732984, 0.420672143, 0.420582619, 0.414704886, 0.397618146, 0.386609633, 0.374601402, 0.334969449, 0.295441387, 0.249040687, 0.218897183, 0.213622456, 0.189758771, 0.169740681, 0.157768951, 0.148698874, 0.145143375, 0.141189382, 0.136867091, 0.134997436, 0.132216107, 0.132712423]
x1 = [-7.946700508, -7.470812183, -6.827411168, -6.279187817, -5.76142132, -5.53680203, -5.354060914, -5.270304569, -5.144670051, -5.038071066, -4.946700508, -4.935279188, -5.121827411, -5.426395939, -5.654822335, -5.711928934, -6.005076142, -6.28680203, -6.48857868, -6.652284264, -6.873096447, -7.158629442, -7.46319797, -7.653553299, -7.779187817, -7.996192893 ]
fig = plt.figure()

ax = fig.add_subplot(111)

ax.invert_yaxis()

fontsize=25
ax.semilogy(np.log10(x),P,'k.-', lw=3, label='Mine')
ax.semilogy(x1,P1,'m.-', lw=3, label='A&M')
#ax.semilogy(x,P,'k.-', lw=3, label='')

ax.legend()
ax.set_xlim([-8, -4])
ax.set_ylim([1, .1])
#ax.set_xlabel('Geometric mean radius [micron]', fontsize=fontsize)
#ax.set_xlabel('alpha', fontsize=fontsize)
ax.set_ylabel('Pressure [bar]', fontsize=fontsize)

#ax2 = ax.twiny()

#ax.plot(np.log10(qc_smooth),P,'k--')

#ax2.set_xlim([110, 160])
fig.savefig('mixingR_Scloud.png',bbox_inches='tight')

plt.show()
