#!/usr/bin/env python

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize as opt
import h5py

# ==============================================================================
def generate_z_value(x,y,z):
	triang = tri.Triangulation(x,y)
	interpolator = tri.LinearTriInterpolator(triang, z)
	Xi, Yi = np.meshgrid(xi, yi)
	zi = interpolator(Xi, Yi)
	return zi

# ==================================================================
# Input Parameters
plt.rcParams.update({'font.size': 12})
ts = 2
ngridx = 21
ngridy = 21
levs = 14
linewidths = 0.25
# ==================================================================

# Load file
f = h5py.File('./output/sol_2dtl.h5','r')

t = f['time'][:]
t = t.ravel()

geo = f['geometry/vertices']

VF = f['vertex_fields']

U = VF['solution_displacement']
P = VF['solution_pressure']
E = VF['solution_tracestrain']

U_E = VF['Exact Solution_displacement']
P_E = VF['Exact Solution_pressure']
E_E = VF['Exact Solution_tracestrain']

# ========================================================================== 80
# Generate grid values
xi = np.linspace(0,1,ngridx)
yi = np.linspace(0,1,ngridy)

# Linearly interpolate data (x, y) on a grid defined by (xi, yi)
x = geo[:][:,0]
y = geo[:][:,1]

# ========================================================================== 80
# Determine z value

# Displacement
x_disp = generate_z_value(x,y,U[ts,:,0])
y_disp = generate_z_value(x,y,U[ts,:,1])

x_disp_exact = generate_z_value(x,y,U_E[ts,:,0])
y_disp_exact = generate_z_value(x,y,U_E[ts,:,1])

# Pressure
pres = generate_z_value(x,y,P[ts,:])
pres_exact = generate_z_value(x,y,P_E[ts,:])

# ========================================================================== 80
# Generate plot, X Displacement
fig, axes = plt.subplots(ncols=1, nrows=3)
fig.set_size_inches(12,10)

axes[0].contour(xi, yi, x_disp, levels=levs, linewidths=linewidths)
cntr0 = axes[0].contourf(xi, yi, x_disp, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr0, ax=axes[0])
axes[0].set(xlim=(0, 1), ylim=(0, 1))
axes[0].set_title(r'$u_{x}$')

axes[1].contour(xi, yi, x_disp_exact, levels=levs, linewidths=linewidths)
cntr1 = axes[1].contourf(xi, yi, x_disp_exact, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr1, ax=axes[1])
axes[1].set(xlim=(0, 1), ylim=(0, 1))
axes[1].set_title(r'$u_{x}^{MMS}$')

axes[2].contour(xi, yi, np.abs(x_disp - x_disp_exact), levels=levs, linewidths=linewidths)
cntr2 = axes[2].contourf(xi, yi, np.abs(x_disp - x_disp_exact), levels=levs, cmap="RdBu_r")
fig.colorbar(cntr2, ax=axes[2])
axes[2].set(xlim=(0, 1), ylim=(0, 1))
axes[2].set_title(r'$|u_{x} - u_{x}^{MMS}|$')

fig.tight_layout()
plt.subplots_adjust(hspace=0.5)
fig.savefig('mms_2dtl_x_disp.png',dpi = 300)
plt.show()

# ========================================================================== 80
# Generate plot, Y Displacement
fig, axes = plt.subplots(ncols=1, nrows=3)
fig.set_size_inches(12,10)

axes[0].contour(xi, yi, y_disp, levels=levs, linewidths=linewidths)
cntr0 = axes[0].contourf(xi, yi, y_disp, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr0, ax=axes[0])
axes[0].set(xlim=(0, 1), ylim=(0, 1))
axes[0].set_title(r'$u_{y}$')

axes[1].contour(xi, yi, y_disp_exact, levels=levs, linewidths=linewidths)
cntr1 = axes[1].contourf(xi, yi, y_disp_exact, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr1, ax=axes[1])
axes[1].set(xlim=(0, 1), ylim=(0, 1))
axes[1].set_title(r'$u_{y}^{MMS}$')

axes[2].contour(xi, yi, np.abs(y_disp - y_disp_exact), levels=levs, linewidths=linewidths)
cntr2 = axes[2].contourf(xi, yi, np.abs(y_disp - y_disp_exact), levels=levs, cmap="RdBu_r")
fig.colorbar(cntr2, ax=axes[2])
axes[2].set(xlim=(0, 1), ylim=(0, 1))
axes[2].set_title(r'$|u_{y} - u_{y}^{MMS}|$')

fig.tight_layout()
plt.subplots_adjust(hspace=0.5)
fig.savefig('mms_2dtl_y_disp.png',dpi = 300)
plt.show()

# ========================================================================== 80
# Generate plot, Y Displacement
fig, axes = plt.subplots(ncols=1, nrows=3)
fig.set_size_inches(12,10)

axes[0].contour(xi, yi, pres, levels=levs, linewidths=linewidths)
cntr0 = axes[0].contourf(xi, yi, pres, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr0, ax=axes[0])
axes[0].set(xlim=(0, 1), ylim=(0, 1))
axes[0].set_title(r'$p$')

axes[1].contour(xi, yi, pres_exact, levels=levs, linewidths=linewidths)
cntr1 = axes[1].contourf(xi, yi, pres_exact, levels=levs, cmap="RdBu_r")
fig.colorbar(cntr1, ax=axes[1])
axes[1].set(xlim=(0, 1), ylim=(0, 1))
axes[1].set_title(r'$p^{MMS}$')

axes[2].contour(xi, yi, np.abs(pres - pres_exact), levels=levs, linewidths=linewidths)
cntr2 = axes[2].contourf(xi, yi, np.abs(pres - pres_exact), levels=levs, cmap="RdBu_r")
fig.colorbar(cntr2, ax=axes[2])
axes[2].set(xlim=(0, 1), ylim=(0, 1))
axes[2].set_title(r'$|p - p^{MMS}|$')

fig.tight_layout()
plt.subplots_adjust(hspace=0.5)
fig.savefig('mms_2dtl_pres.png',dpi = 300)
plt.show()








