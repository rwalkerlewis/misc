
import numpy
import numpy as np

x1 = numpy.arange(-0.1, 1.01, 0.01)
y1 = numpy.arange(-0.1, 1.01, 0.01)
z1 = numpy.arange(-0.1, 1.01, 0.01)
x, y, z = numpy.meshgrid(x1, y1, z1)

xyz = numpy.zeros((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
xyz_data = numpy.ones((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
xyz[:, 0] = x.ravel()
xyz[:, 1] = y.ravel()
xyz[:, 2] = z.ravel()

R_func = lambda theta,phi: numpy.array( [ [numpy.sin(theta)*numpy.cos(phi),   numpy.sin(phi)*numpy.sin(phi),         numpy.cos(theta)],
					                      [numpy.cos(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.sin(phi),        -numpy.sin(theta)],
                                          [                -numpy.sin(phi),                  numpy.cos(phi),                        0] ] )				   

npts = 15
ntpts = 2
dim = 3

theta = numpy.deg2rad(numpy.linspace(0,30,15))
phi = numpy.deg2rad(numpy.linspace(0,14,15))



sigma_rr = np.ones([2,15])
sigma_pp = np.ones([2,15])
sigma_tt = np.ones([2,15])




zeros = numpy.zeros(15)

#sigma_r = np.array( [[sigma_rr, zeros, zeros],
#                     [zeros, sigma_phiphi, zeros],
#                     [zeros, zeros, sigma_thetatheta]])


R = numpy.zeros([npts, dim, dim])
R[:,0,0] =  numpy.sin(theta)*numpy.cos(phi)
R[:,0,1] =  numpy.sin(phi)*numpy.sin(phi)
R[:,0,2] =  numpy.cos(theta)
R[:,1,0] =  numpy.cos(theta)*numpy.cos(phi)
R[:,1,1] =  numpy.cos(theta)*numpy.sin(phi)
R[:,1,2] = -numpy.sin(theta)
R[:,2,0] = -numpy.sin(phi)
R[:,2,1] =  numpy.cos(phi)
R[:,2,2] =  numpy.zeros(theta.size)

sigma_sph = numpy.zeros([ntpts, npts, dim, dim])
sigma_sph[:,:,0,0] = sigma_rr[:,:]
sigma_sph[:,:,1,1] = sigma_pp[:,:]
sigma_sph[:,:,2,2] = sigma_tt[:,:]
#R = numpy.array( [ [numpy.sin(theta)*numpy.cos(phi),   numpy.sin(phi)*numpy.sin(phi),         numpy.cos(theta)],
#				   [numpy.cos(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.sin(phi),        -numpy.sin(theta)],
#				   [                -numpy.sin(phi),                  numpy.cos(phi),  numpy.zeros(theta.size)] ] )
				   




