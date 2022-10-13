#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize as opt
plt.rcParams.update({'font.size': 19})

# ==============================================================================
# [Displacement, Trace Strain, Pressure]

mesh_title = "Mesh Convergence"
N_StT1 = np.array( [[    450,      81,      49],
			  		[   1922,     289,     225],
	   	  			[   7938,    1089,     961],
			   		[  32258,    4225,    3969] ] )

L_2_StT1 = np.array( [[0.0156383, 0.188726, 0.170848],
					  [0.00411089, 0.0464749, 0.0471813],
   					  [0.00105209, 0.011514, 0.0122751],
					  [0.000265451, 0.00286479, 0.00311006]] )
					  
slope_u_StT1 = np.polyfit(np.log10(N_StT1[:,0]),np.log10(L_2_StT1[:,0]),1)
slope_e_StT1 = np.polyfit(np.log10(N_StT1[:,1]),np.log10(L_2_StT1[:,1]),1)
slope_p_StT1 = np.polyfit(np.log10(N_StT1[:,2]),np.log10(L_2_StT1[:,2]),1)

poly1d_fn_u_StT1 = np.poly1d(slope_u_StT1) 
poly1d_fn_e_StT1 = np.poly1d(slope_e_StT1) 
poly1d_fn_p_StT1 = np.poly1d(slope_p_StT1) 

u_fit_StT1 = lambda x: np.power(10,poly1d_fn_u_StT1(np.log10(x)))
e_fit_StT1 = lambda x: np.power(10,poly1d_fn_e_StT1(np.log10(x)))
p_fit_StT1 = lambda x: np.power(10,poly1d_fn_p_StT1(np.log10(x)))
					  					 

temp_title = "Temporal Convergence"


N_S2Tt = np.array( [ [     10,      10,      10],
					 [     20,      20,      20],
					 [     40,      40,      40],
        			 [     80,      80,      80] ] )

L_2_S2Tt = np.array( [ [0.000100628, 0.000695106, 0.00187617],
					   [5.03739e-05, 0.000348023, 0.000938252],
					   [2.51628e-05, 0.000173855, 0.000468513],
					   [1.25706e-05, 8.68549e-05, 0.000234024] ] )

slope_u_S2Tt = np.polyfit(np.log10(N_S2Tt[:,0]),np.log10(L_2_S2Tt[:,0]),1)
slope_e_S2Tt = np.polyfit(np.log10(N_S2Tt[:,1]),np.log10(L_2_S2Tt[:,1]),1)
slope_p_S2Tt = np.polyfit(np.log10(N_S2Tt[:,2]),np.log10(L_2_S2Tt[:,2]),1)

poly1d_fn_u_S2Tt = np.poly1d(slope_u_S2Tt) 
poly1d_fn_e_S2Tt = np.poly1d(slope_e_S2Tt) 
poly1d_fn_p_S2Tt = np.poly1d(slope_p_S2Tt) 

u_fit_S2Tt = lambda x: np.power(10,poly1d_fn_u_S2Tt(np.log10(x)))
e_fit_S2Tt = lambda x: np.power(10,poly1d_fn_e_S2Tt(np.log10(x)))
p_fit_S2Tt = lambda x: np.power(10,poly1d_fn_p_S2Tt(np.log10(x)))
					   
# ==============================================================================
# Mesh Plot
fig, axes = plt.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)

axes.loglog( N_StT1[:,0], L_2_StT1[:,0], '-b', label=r'$||u-u^{*}||_{2}$')
axes.loglog(N_StT1[:,0], u_fit_StT1(N_StT1[:,0]),'--b', label='slope: ' + str(np.round(slope_u_StT1[0],3)) )

axes.loglog( N_StT1[:,1], L_2_StT1[:,1], '-m', label=r'$||\epsilon-\epsilon^{*}||_{2}$')
axes.loglog(N_StT1[:,1], e_fit_StT1(N_StT1[:,1]),'--m', label='slope: ' + str(np.round(slope_u_StT1[0],3)) )

axes.loglog( N_StT1[:,2], L_2_StT1[:,2], '-g', label=r'$||p-p^{*}||_{2}$')
axes.loglog(N_StT1[:,2], p_fit_StT1(N_StT1[:,2]),'--g', label='slope: ' + str(np.round(slope_u_StT1[0],3)) )

axes.grid(True, which="both", ls="-")
axes.set_title(mesh_title)

axes.legend()
axes.set_xlabel("# unknowns")
axes.set_ylabel("Error")
fig.tight_layout()

fig.savefig('Mesh_Conv_2DTL.png',dpi = 300)
plt.show()

# ==============================================================================
# Temporal Plot
fig, axes = plt.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)

axes.loglog( N_S2Tt[:,0], L_2_S2Tt[:,0], '-b', label=r'$||u-u^{*}||_{2}$')
axes.loglog(N_S2Tt[:,0], u_fit_S2Tt(N_S2Tt[:,0]),'--b', label='slope: ' + str(np.round(slope_u_S2Tt[0],3)) )

axes.loglog( N_S2Tt[:,1], L_2_S2Tt[:,1], '-m', label=r'$||\epsilon-\epsilon^{*}||_{2}$')
axes.loglog(N_S2Tt[:,1], e_fit_S2Tt(N_S2Tt[:,1]),'--m', label='slope: ' + str(np.round(slope_u_S2Tt[0],3)) )

axes.loglog( N_S2Tt[:,2], L_2_S2Tt[:,2], '-g', label=r'$||p-p^{*}||_{2}$')
axes.loglog(N_S2Tt[:,2], p_fit_S2Tt(N_S2Tt[:,2]),'--g', label='slope: ' + str(np.round(slope_u_S2Tt[0],3)) )

axes.grid(True, which="both", ls="-")
axes.set_title(temp_title)

axes.legend()
axes.set_xlabel(r'$\frac{1}{dt}$')
axes.set_ylabel("Error")
fig.tight_layout()

fig.savefig('Temp_Conv_2DQT.png',dpi = 300)
plt.show()

