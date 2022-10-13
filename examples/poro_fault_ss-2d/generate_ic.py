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

# Dimensions
DOMAIN_X = 1.0
DOMAIN_Y = 1.0

# Discretization
dx = 0.005
dy = 0.005 

nx = numpy.int64(DOMAIN_X / dx)
ny = numpy.int64(DOMAIN_Y / dy)

# Grid
x = numpy.linspace(-DOMAIN_X/2,DOMAIN_X/2,nx+1)
y = numpy.linspace(-DOMAIN_Y/2,DOMAIN_Y/2,ny+1)

xx, yy = numpy.meshgrid(x,y)

# Initial time
t0 = 0.0


# Two Dimensional Values
u_x = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy - 1)
u_y = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy)
p = t0*(yy*yy - 1)
trace_strain = numpy.nan_to_num(yy / numpy.abs(yy), nan=0.0) * -2 * yy
L_x = numpy.ones(x.size) * 0
L_y = numpy.ones(y.size) * t0
p_f = numpy.ones(x.size) * 0.5 * t0

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

class GenerateDB(object):


    def run(self):
        """Generate the database.
        """
        
        # Domain, centroid
        x = numpy.linspace(-DOMAIN_X/2,DOMAIN_X/2,nx+1)
        y = numpy.linspace(-DOMAIN_Y/2,DOMAIN_Y/2,ny+1)
        xx, yy = numpy.meshgrid(x,y)
        xy = numpy.column_stack((xx.flatten(), yy.flatten()))
        
        # xx = x.reshape([1, x.shape[0]]) * numpy.ones([ny ,1])
        # yy = y.reshape([y.shape[0], 1]) * numpy.ones([1, nx])
        # xy[:, 0] = numpy.ravel(xx)
        # xy[:, 1] = numpy.ravel(numpy.transpose(yy))
        # xy[:, 0] = x2.ravel()
        # xy[:, 1] = y2.ravel()

        # x2, y2 = numpy.meshgrid(x,y)                
        # xy = numpy.column_stack((
        #       x2.flatten(),
        #       y2.flatten()
        # ))
        
        
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        displacement_xGrid = {"name": "displacement_x",
                        "units": "m",
                        "data": numpy.ravel(u_x)}

        displacement_yGrid = {"name": "displacement_y",
                               "units": "m",
                               "data": numpy.ravel(u_y)}        

        pressureGrid = {"name": "pressure",
                               "units": "Pa",
                               "data": numpy.ravel(p)}

        trace_strainGrid = {"name": "trace_strain",
                            "units": "none",
                            "data": numpy.ravel(trace_strain)}

        data = {"num-x": xx.shape[0],
                "num-y": yy.shape[0],
                "points": xy,                
                "x": xx.flatten(),
                "y": yy.flatten(),
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_xGrid, displacement_yGrid, pressureGrid,
                trace_strainGrid]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("MMS_ic.spatialdb")
        io.write(data)        
        return

# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file

