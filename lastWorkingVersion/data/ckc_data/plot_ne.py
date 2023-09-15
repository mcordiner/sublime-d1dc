import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})

def Telec(r, Q, Tkin):
    Tmax = 1e4
    rcs = 1.125e6 * pow(Q/1e29,0.75)
    if (r < rcs):
        Te = Tkin
    elif (r > 2.*rcs):
        Te = Tmax
    else: 
        Te = Tkin + (Tmax - Tkin)*((r/rcs)-1.)
        
    return Te

def nelec(r, Q, vexp, Te, rH, xne):
    kion = 4.1e-7
    krec = 3e-13 * np.sqrt(300./Te) #/* Recombination rate From RATE12, accountiong for cm3 to m3 conversion */
    Rrec = 3.2e6 * np.sqrt(Q/1e29)
    #/* Equation 5 of Zakharov 2007 */
    ne = xne * np.sqrt(Q*kion/vexp/krec/(rH*rH)) * np.power((Te/300.),0.15) * (Rrec/(r*r)) * (1-np.exp(-r/Rrec)) + (5e6/(rH*rH))
    
    return ne

RH = [0.50,1.00,1.50,2.00,2.50]
Q0 = [1e27,5e27,1e28,5e28,1e29,5e29,1e30]
    
Tkin = 50.
vexp = 0.5

for Q in Q0:
	cid=0
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for rH in RH:
		TeCKC = np.genfromtxt("Te_Q%s_H%s.dat" %(Q,rH),skip_header=4)
		neCKC = np.genfromtxt("ne_Q%s_H%s.dat" %(Q,rH),skip_header=4)
		CKC = np.column_stack((TeCKC,neCKC[:,1]))
		i = 0
		while i < len(CKC[:,0]):
			CKC[i,0] = CKC[i,0]*1000
			CKC[i,2] = CKC[i,2]*1e6
			i += 1
		
		# Calculate values analytically for Xne = 0.2
		ANA02 = np.empty(shape = CKC.shape, dtype=float)
		i=0
		xne = 0.2
		while i < len(ANA02[:,0]):
			ANA02[i,0] = CKC[i,0]
			ANA02[i,1] = Telec(ANA02[i,0], Q, Tkin)
			ANA02[i,2] = nelec(ANA02[i,0], Q, vexp, ANA02[i,1], rH, xne)
			i += 1
		
		ax.plot(ANA02[:,0]/1000.,ANA02[:,2],'--', color = 'C%d'%(cid), label=r'Analytic: H=%s'%(rH))
		ax.plot(CKC[:,0]/1000.,CKC[:,2],'-', color = 'C%d'%(cid), label=r'CKC: H=%s'%(rH))
		cid += 1

	plt.legend(loc="best",ncol=2,numpoints=1)
	plt.xlim(3e0,1.1e5)

	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel('Radius (km)')
	ax.set_ylabel(r'$n_{elec}$ $(m^{-3})$')
	
	plt.title(r'$n_{elec}$ Comparison $Q = $%.1e'%(Q))

	filename = "Q%s" %(Q)
	plt.savefig('plots/ne_' + filename + '.png')
	#plt.show()
	
# for rH in RH:
# 	cid=0
# 	fig = plt.figure()
# 	ax = fig.add_subplot(111)
# 	for Q in Q0:
# 		TeCKC = np.genfromtxt("Te_Q%s_H%s.dat" %(Q,rH),skip_header=4)
# 		neCKC = np.genfromtxt("ne_Q%s_H%s.dat" %(Q,rH),skip_header=4)
# 		CKC = np.column_stack((TeCKC,neCKC[:,1]))
# 		i = 0
# 		while i < len(CKC[:,0]):
# 			CKC[i,0] = CKC[i,0]*1000
# 			CKC[i,2] = CKC[i,2]*1e6
# 			i += 1
# 		
# 		# Calculate values analytically for Xne = 0.2
# 		ANA02 = np.empty(shape = CKC.shape, dtype=float)
# 		i=0
# 		xne = 0.2
# 		while i < len(ANA02[:,0]):
# 			ANA02[i,0] = CKC[i,0]
# 			ANA02[i,1] = Telec(ANA02[i,0], Q, Tkin)
# 			ANA02[i,2] = nelec(ANA02[i,0], Q, vexp, ANA02[i,1], rH, xne)
# 			i += 1
# 		
# 		ax.plot(ANA02[:,0]/1000.,ANA02[:,2],'--', color = 'C%d'%(cid),label=r'Analytic: Q=%s'%(Q))
# 		ax.plot(CKC[:,0]/1000.,CKC[:,2],'-', color = 'C%d'%(cid), label=r'CKC: Q=%s'%(Q))
# 		cid += 1
# 
# 	plt.legend(loc="best",ncol=3,numpoints=1)
# 	plt.xlim(3e0,1.1e5)
# 
# 	ax.set_xscale("log")
# 	ax.set_yscale("log")
# 
# 	ax.set_xlabel('Radius (km)')
# 	ax.set_ylabel(r'$n_{elec}$ $(m^{-3})$')
# 	
# 	plt.title(r'$n_{elec}$ Comparison $H = $%.1e'%(rH))
# 
# 	filename = "H%s" %(rH)
# 	plt.savefig('plots/ne_' + filename + '.png')
# 	#plt.show()
