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

from time import time
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
xy = numpy.column_stack((xx.flatten(), yy.flatten()))


# Initial time
t0 = 0.0


# Two Dimensional Values
u_x = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy - 1)
u_y = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy)
# p = t0*(yy*yy - 1)
p = (yy*yy - 1)
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

# Boundary Values

# xneg
xneg = numpy.column_stack((xx[:,0].flatten(), yy[:,0].flatten()))
u_x_xneg = u_x[:,0]
u_y_xneg = u_y[:,0]
p_xneg = p[:,0]
trace_strain_xneg = trace_strain[:,0]

# xpos
xpos = numpy.column_stack((xx[:,-1].flatten(), yy[:,-1].flatten()))

u_x_xpos = u_x[:,-1]
u_y_xpos = u_y[:,-1]
p_xpos = p[:,-1]
trace_strain_xpos = trace_strain[:,-1]

# yneg
yneg = numpy.column_stack((xx[0,:].flatten(), yy[0,:].flatten()))
u_x_yneg = u_x[0,:]
u_y_yneg = u_y[0,:]
p_yneg = p[0,:]
trace_strain_yneg = trace_strain[0,:]

# ypos
ypos = numpy.column_stack((xx[-1,:].flatten(), yy[-1,:].flatten()))
u_x_ypos = u_x[-1,:]
u_y_ypos = u_y[-1,:]
p_ypos = p[-1,:]
trace_strain_ypos = trace_strain[-1,:]

class GenerateDB(object):
    def run(self):
        """Generate the databases.
        """
        # displacement
        self.xneg_disp()
        self.xpos_disp()
        self.yneg_disp()
        self.ypos_disp()

        # pressure
        self.xneg_pres()
        self.xpos_pres()
        self.yneg_pres()
        self.ypos_pres()

        # time history
        self.time_history()

        return

    def xneg_disp(self):
        """Generate the database for xneg displacement.
        """
        # X Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        displacement_x = {"name": "initial_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(u_x_xneg)}

        displacement_y = {"name": "initial_amplitude_y",
                               "units": "m",
                               "data": numpy.ravel(u_y_xneg)}

        data = {"num-x": xneg.shape[0],
                "num-y": 1,
                "points": xneg,                
                "x": x[0],
                "y": y.flatten(),
                "coordsys": cs,
                "data_dim": 1,
                "values": [displacement_x, displacement_y]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xneg_disp.spatialdb")
        io.write(data)        
        return

    def xpos_disp(self):
        """Generate the database for xpos displacement.
        """
        # X Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        displacement_x = {"name": "initial_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(u_x_xpos)}

        displacement_y = {"name": "initial_amplitude_y",
                               "units": "m",
                               "data": numpy.ravel(u_y_xpos)}

        data = {"num-x": xpos.shape[0],
                "num-y": 1,
                "points": xpos,                
                "x": x[-1],
                "y": y.flatten(),
                "coordsys": cs,
                "data_dim": 1,
                "values": [displacement_x, displacement_y]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xpos_disp.spatialdb")
        io.write(data)        
        return

    def yneg_disp(self):
        """Generate the database for yneg displacement.
        """
        # Y Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        displacement_x = {"name": "initial_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(u_x_yneg)}

        displacement_y = {"name": "initial_amplitude_y",
                               "units": "m",
                               "data": numpy.ravel(u_y_yneg)}

        data = {"num-x": 1,
                "num-y": yneg.shape[0],
                "points": yneg,                
                "x": x.flatten(),
                "y": y[0],
                "coordsys": cs,
                "data_dim": 1,
                "values": [displacement_x, displacement_y]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_yneg_disp.spatialdb")
        io.write(data)        
        return

    def ypos_disp(self):
        """Generate the database for ypos displacement.
        """
        # Y Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        displacement_x = {"name": "initial_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(u_x_ypos)}

        displacement_y = {"name": "initial_amplitude_y",
                               "units": "m",
                               "data": numpy.ravel(u_y_ypos)}

        data = {"num-x": 1,
                "num-y": ypos.shape[0],
                "points": ypos,                
                "x": x.flatten(),
                "y": y[-1],
                "coordsys": cs,
                "data_dim": 1,
                "values": [displacement_x, displacement_y]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_ypos_disp.spatialdb")
        io.write(data)        
        return

    def xneg_pres(self):
        """Generate the database for xneg pressure.
        """
        # X Neg Pressure
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'm'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        pressure = {"name": "initial_amplitude",
                        "units": "Pa",
                        "data": numpy.zeros(p_xneg.shape[0]).ravel()}

        time_history_start = {"name": "time_history_start_time",
                        "units": "second",
                        "data": numpy.zeros(p_xneg.shape[0]).ravel()}

        time_history_amplitude = {"name": "time_history_amplitude",
                        "units": "Pa",
                        "data": numpy.ravel(p_xneg)}

        data = {"num-x": xneg.shape[0],
                "num-y": 1,
                "points": xneg,                
                "x": x[0],
                "y": y.flatten(),
                "coordsys": cs,
                "data_dim": 1,
                "values": [pressure, time_history_start, time_history_amplitude]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xneg_pres.spatialdb")
        io.write(data)        
        return

    def xpos_pres(self):
        """Generate the database for xpos pressure.
        """
        # X Pos Pressure
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'm'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        pressure = {"name": "initial_amplitude",
                        "units": "Pa",
                        "data": numpy.zeros(p_xpos.shape[0]).ravel()}

        time_history_start = {"name": "time_history_start_time",
                        "units": "second",
                        "data": numpy.zeros(p_xpos.shape[0]).ravel()}

        time_history_amplitude = {"name": "time_history_amplitude",
                        "units": "Pa",
                        "data": numpy.ravel(p_xpos)}

        data = {"num-x": xpos.shape[0],
                "num-y": 1,
                "points": xpos,                
                "x": x[-1],
                "y": y.flatten(),
                "coordsys": cs,
                "data_dim": 1,
                "values": [pressure, time_history_start, time_history_amplitude]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xpos_pres.spatialdb")
        io.write(data)        
        return

    def yneg_pres(self):
        """Generate the database for yneg pressure.
        """
        # Y Neg Pressure
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'm'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        pressure = {"name": "initial_amplitude",
                        "units": "Pa",
                        "data": numpy.zeros(p_yneg.shape[0]).ravel()}

        time_history_start = {"name": "time_history_start_time",
                        "units": "second",
                        "data": numpy.zeros(p_yneg.shape[0]).ravel()}

        time_history_amplitude = {"name": "time_history_amplitude",
                        "units": "Pa",
                        "data": numpy.ravel(p_yneg)}

        data = {"num-x": 1,
                "num-y": yneg.shape[0],
                "points": yneg,                
                "x": x.flatten(),
                "y": y[0],
                "coordsys": cs,
                "data_dim": 1,
                "values": [pressure, time_history_start, time_history_amplitude]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_yneg_pres.spatialdb")
        io.write(data)        
        return

    def ypos_pres(self):
        """Generate the database for ypos pressure.
        """
        # Y Pos Pressure
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'm'
        cs.inventory.spaceDim = 2
        cs._configure()
        
        pressure = {"name": "initial_amplitude",
                        "units": "Pa",
                        "data": numpy.zeros(p_ypos.shape[0]).ravel()}

        time_history_start = {"name": "time_history_start_time",
                        "units": "second",
                        "data": numpy.zeros(p_ypos.shape[0]).ravel()}

        time_history_amplitude = {"name": "time_history_amplitude",
                        "units": "Pa",
                        "data": numpy.ravel(p_ypos)}

        data = {"num-x": 1,
                "num-y": ypos.shape[0],
                "points": ypos,                
                "x": x.flatten(),
                "y": y[-1],
                "coordsys": cs,
                "data_dim": 1,
                "values": [pressure, time_history_start, time_history_amplitude]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_ypos_pres.spatialdb")
        io.write(data)        
        return

    def time_history(self):
        """Generate the time history database.
        """
        tmax = 1
        nsteps = 2001
        time = numpy.zeros(nsteps+1)
        time[:-1] = numpy.linspace(0,tmax,nsteps)
        time[-1] = 999.000
        amplitude = time.copy()
        from spatialdata.spatialdb.TimeHistoryIO import write
        write(time=time, amplitude=amplitude, units="second", filename="time_history.timedb")        
        return

# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file

