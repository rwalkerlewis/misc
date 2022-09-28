#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize as opt
plt.rcParams.update({'font.size': 22})

# ==============================================================================
# [Displacement, Trace Strain, Pressure]

N = np.array( [ [    129,      34,      32],
				[    515,      99,      96],
				[   2055,     325,     320],
				[   8207,    1161,    1152],
				[  32799,    4369,    4352] ] )

#       > L_2 convergence rate: [1.7, 1.2, 1.1]
L_2 = np.array( [ [0.013658, 0.0719259, 0.959572],
				  [0.00481491, 0.0457775, 0.609139],
				  [0.00159355, 0.0248645, 0.330599],
				  [0.00047527, 0.0129663, 0.171391],
				  [0.000123604, 0.00385394, 0.0534565] ] )
				  
slope = np.zeros(3)

slope_u = np.polyfit(np.log10(N[:,0]),np.log10(L_2[:,0]),1)
slope_e = np.polyfit(np.log10(N[:,1]),np.log10(L_2[:,1]),1)
slope_p = np.polyfit(np.log10(N[:,2]),np.log10(L_2[:,2]),1)

poly1d_fn_u = np.poly1d(slope_u) 
poly1d_fn_e = np.poly1d(slope_e) 
poly1d_fn_p = np.poly1d(slope_p) 

u_fit = lambda x: np.power(10,poly1d_fn_u(np.log10(x)))
e_fit = lambda x: np.power(10,poly1d_fn_e(np.log10(x)))
p_fit = lambda x: np.power(10,poly1d_fn_p(np.log10(x)))

mesh_title = "Terzaghi Problem Mesh Convergence"

					   
# ==============================================================================
# Mesh Plot
fig, axes = plt.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)

axes.loglog( N[:,0], L_2[:,0], '-b', label=r'$||u-u^{*}||_{2}$')
axes.loglog(N[:,0], u_fit(N[:,0]),'--b', label='slope: ' + np.str(np.round(slope_u[0],3)) )
axes.loglog( N[:,1], L_2[:,1], '-m', label=r'$||\epsilon-\epsilon^{*}||_{2}$')
axes.loglog(N[:,1], e_fit(N[:,1]),'--m', label='slope: ' + np.str(np.round(slope_e[0],3)) )
axes.loglog( N[:,2], L_2[:,2], '-g', label=r'$||p-p^{*}||_{2}$')
axes.loglog(N[:,2], p_fit(N[:,2]),'--g', label='slope: ' + np.str(np.round(slope_p[0],3)) )
axes.grid(True, which="both", ls="-")
axes.set_title(mesh_title)

axes.legend()
axes.set_xlabel(r'Field Size $N$')
axes.set_ylabel(r'Error $\epsilon$')
fig.tight_layout()

fig.savefig('Terzaghi_Problem_Mesh_Conv.png',dpi = 300)
plt.show()


