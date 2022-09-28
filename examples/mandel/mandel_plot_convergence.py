#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize as opt
plt.rcParams.update({'font.size': 22})

# ==============================================================================
# [Displacement, Trace Strain, Pressure]

N = np.array( [ [   2046,     289,     255],  
			 	[   8190,    1089,    1023],
				[  32766,    4225,    4095],
				[ 131070,   16641,   16383]] )

L_2 = np.array( [ [0.414314, 8.36831, 34.8587],
				  [0.13899, 4.68244, 18.0136],
				  [0.0390221, 1.89597, 7.40731],
				  [0.0168899, 1.39359, 3.13863] ] )

#       > L_2 convergence rate: [1.6, 0.93, 1.2]

slope_u = np.polyfit(np.log10(N[:,0]),np.log10(L_2[:,0]),1)
slope_e = np.polyfit(np.log10(N[:,1]),np.log10(L_2[:,1]),1)
slope_p = np.polyfit(np.log10(N[:,2]),np.log10(L_2[:,2]),1)

poly1d_fn_u = np.poly1d(slope_u) 
poly1d_fn_e = np.poly1d(slope_e) 
poly1d_fn_p = np.poly1d(slope_p) 

u_fit = lambda x: np.power(10,poly1d_fn_u(np.log10(x)))
e_fit = lambda x: np.power(10,poly1d_fn_e(np.log10(x)))
p_fit = lambda x: np.power(10,poly1d_fn_p(np.log10(x)))

mesh_title = "Mandel Problem Mesh Convergence"

					   
# ==============================================================================
# Mesh Plot
fig, axes = plt.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)

axes.loglog( N[:,0], L_2[:,0], '-b', label=r'$||u-u^{*}||_{2}$')
axes.loglog(N[:,0], u_fit(N[:,0]),'--b', label='slope: ' + str(np.round(slope_u[0],3)) )
axes.loglog( N[:,1], L_2[:,1], '-m', label=r'$||\epsilon-\epsilon^{*}||_{2}$')
axes.loglog(N[:,1], e_fit(N[:,1]),'--m', label='slope: ' + str(np.round(slope_e[0],3)) )
axes.loglog( N[:,2], L_2[:,2], '-g', label=r'$||p-p^{*}||_{2}$')
axes.loglog(N[:,2], p_fit(N[:,2]),'--g', label='slope: ' + str(np.round(slope_p[0],3)) )
axes.grid(True, which="both", ls="-")
axes.set_title(mesh_title)

axes.legend()
axes.set_xlabel(r'Field Size $N$')
axes.set_ylabel(r'Error $\epsilon$')
fig.tight_layout()

fig.savefig('Mandel_Problem_Mesh_Conv.png',dpi = 300)
plt.show()


