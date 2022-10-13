#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import numpy 
from netCDF4 import Dataset
import scipy.ndimage
import matplotlib.pyplot as plt

# mesh = Dataset('mesh_hex.exo','r')

# SPE 10 Parameters
nx = 60
ny = 220
nz = 85

dx = 20.
dy = 10. 
dz = 2. 

# Select which layer in xy plane
zval = 0

# phi = numpy.loadtxt('spe_phi.dat').ravel()
# phi = phi.reshape([nz,ny,nx])

# perm = numpy.loadtxt('spe_perm.dat').ravel() * 9.869233e-16 

#k_x = perm[0: nz*ny*nx].reshape([nz,ny,nx])
#k_y = perm[nz*ny*nx:2*nz*ny*nx].reshape([nz,ny,nx])

# k_x = k_x[zval,:,:]
# k_y = k_y[zval,:,:]

# k_tensor = numpy.zeros([nx*ny,4])
# k_tensor[:,0] = k_x.flatten()
# k_tensor[:,1] = k_y.flatten()
# k_tensor[:,2] = numpy.zeros(nx*ny)
# k_tensor[:,3] = numpy.zeros(nx*ny)

# phi = phi[zval,:,:]

# Checkerboard
factor = 20
nx_check = numpy.int32(nx/factor)
ny_check = numpy.int32(ny/factor)

check = numpy.indices([ny_check,nx_check]).sum(axis=0) % 2
check = scipy.ndimage.zoom(check, factor, order=0)
plt.imsave('check.png',check)

# Dimensions
DOMAIN_X = 150.0e+3
DOMAIN_Y = 100.0e+3

# Two Dimensional Values


fluid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
porosity = 0.5 * numpy.ones(nx*ny) # kg / m**3
biot_coefficient = 1.0 * numpy.ones(nx*ny)
fluid_viscosity = 1.0 * numpy.ones(nx*ny) # Pa*s
shear_modulus = 0.5 * numpy.ones(nx*ny) # Pa
fluid_bulk_modulus = 2e9 * numpy.ones(nx*ny) # Pa
drained_bulk_modulus = 10e9 * numpy.ones(nx*ny) # Pa
solid_bulk_modulus = 11039657020.4 * numpy.ones(nx*ny) # Pa




# End of file

