#!/usr/bin/env/ python

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import h5py
import numpy

np.set_printoptions(suppress=True)
# ==============================================================================
# Computational Values
ITERATIONS = 88
EPS = 1e-20

# ==============================================================================
# Physical properties
p_shear_modulus = 3.0
p_solid_density = 2500
p_fluid_density = 1000
p_fluid_bulk_modulus = 8.0
p_solid_bulk_modulus = 10.0
p_drained_bulk_modulus = 4.0
p_biot_coefficient = 0.6
p_porosity = 0.1
p_isotropic_permeability = 1.5
p_fluid_viscosity = 1.0

G = p_shear_modulus
rho_s = p_solid_density
rho_f = p_fluid_density
K_fl = p_fluid_bulk_modulus
K_sg = p_solid_bulk_modulus
K_d = p_drained_bulk_modulus
alpha = p_biot_coefficient
phi = p_porosity
k = p_isotropic_permeability
mu_f = p_fluid_viscosity
P_0 = 1.0
R_0 = 1.0
ndim = 3


M    = 1.0 / ( phi / K_fl + (alpha - phi) /K_sg)
kappa = k/mu_f
K_u = K_d + alpha*alpha*M
S = (3*K_u + 4*G) / (M*(3*K_d + 4*G)) #(1/M) + ( (3*alpha*alpha) / (3*K_d + 4*G) )#
c = kappa / S
nu = (3*K_d - 2*G) / (2*(3*K_d + G))
nu_u = (3*K_u - 2*G) / (2*(3*K_u + G))
U_R_inf = -1.*(P_0*R_0*(1.-2.*nu))/(2.*G*(1.+nu))
U_R_zero = -1.*(P_0*R_0*(1.-2.*nu_u))/(2.*G*(1.+nu_u))
eta = alpha*(1-2*nu) / (2*(1-nu))

# For Center
grav = 9.80665
S_c = phi*(1/K_fl) + (alpha - phi)*(1/K_sg)
eta_c = ( (K_d + 4/3*G)/(2*G) ) * (1 + (K_d*S)/(alpha*alpha))
gamma_f_c = S_c*phi*rho_f*grav + (1-phi)*rho_s*grav
c_v_c = (k / gamma_f_c) * ( (K_d + 4/3*G)/(alpha*alpha + S*(K_d + 4/3*G) ) )

# ==============================================================================
# Input parameter grid
x1 = numpy.arange(-0.0, 1.01, 0.01)
y1 = numpy.arange(-0.0, 1.01, 0.01)
z1 = numpy.arange(-0.0, 1.01, 0.01)
x, y, z = numpy.meshgrid(x1, y1, z1)

xyz = numpy.zeros((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
xyz_data = numpy.ones((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
xyz[:, 0] = x.ravel()
xyz[:, 1] = y.ravel()
xyz[:, 2] = z.ravel()

# ==============================================================================
# Generate positive solutions to characteristic equation

def cryer_zeros_python(nu, nu_u, n_series=50):

    f      = lambda x: np.tan(np.sqrt(x)) - (6*(nu_u - nu)*np.sqrt(x))/(6*(nu_u - nu) - (1 - nu)*(1 + nu_u)*x) # Compressible Constituents 

    a_n = np.zeros(n_series) # initializing roots array
    xtol = 1e-30
    rtol = 1e-15
    for i in range(1,n_series+1):
        a = np.square(i*np.pi) - (i+1)*np.pi
        b = np.square(i*np.pi) + (i+1)*np.pi
        # print('a: ',a)
        # print('b: ',b) 
        f_c = 10
        f_c_old = 0
        rtol_flag = False
        c = 0
        it = 0
        while np.abs(f_c) > xtol and rtol_flag == False:
            c = (a + b) / 2
            f_c = f(c)

            # print('c: ',c)
            # print('f(c):',f_c)
                        
            if f(a)*f_c < 0:
                a = a
                b = c
            elif f(b)*f_c < 0:
                a = c
                b = b     
            else:
                print('Bisection method failed')
                # print('a: ',a)
                # print('b: ',b)
                break
            if np.abs(np.abs(f_c_old) - np.abs(f_c)) < rtol:
                rtol_flag = True                    
            it += 1
            # print('f(c): ',f_c)
            # print('rtol: ',np.abs(np.abs(f_c_old) - np.abs(f_c)))
            # print('rtol flag: ',rtol_flag)
            f_c_old = f_c
        # print('n: ',i)                
        # print('c: ',c)
        # print('f(c):',f_c)
        # print('iter: ',it)
        a_n[i-1] = c               
        
    return(a_n)    
                
def cryer_zeros(nu, nu_u,n_series=50):
    """
    This is somehow tricky, we have to solve the equation numerically in order to
    find all the positive solutions to the equation. Later we will use them to 
    compute the infinite sums. Experience has shown that 200 roots are more than enough to
    achieve accurate results. Note that we find the roots using the bisection method.
    """
    f      = lambda x: np.tan(np.sqrt(x)) - (6*(nu_u - nu)*np.sqrt(x))/(6*(nu_u - nu) - (1 - nu)*(1 + nu_u)*x) # Compressible Constituents 
#    f      = lambda x: np.tan(np.sqrt(x)) - (2*(1-2*nu)*np.sqrt(x))/(2*(1-2*nu) - (1-nu)*x) # Incompressible Constituents

    a_n = np.zeros(n_series) # initializing roots array
    for i in range(1,len(a_n)+1):
        a1 = np.square(i*np.pi) - (i+1)*np.pi
        a2 = np.square(i*np.pi) + (i+1)*np.pi        
        a_n[i-1] = opt.bisect( f,                           # function
                               a1,                          # left point 
                               a2,                          # right point (a tiny bit less than pi/2)
                               xtol=1e-30,                  # absolute tolerance
                               rtol=1e-15                   # relative tolerance
                           )  
    
    return a_n

def pressure(locs, tsteps, x_n):
    """
    Compute pressure field at locations.
    """
    (npts, dim) = locs.shape
    ntpts = tsteps.shape[0]
    pressure = np.zeros((ntpts, npts), dtype=np.float64)
    
    center = np.where(~locs.any(axis=1))[0]
    R = np.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
    R_star = R.reshape([R.size,1]) / R_0
    x_n.reshape([1,x_n.size])

    E = np.square(1-nu)*np.square(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    
    t_track = 0

    for t in tsteps:
        t_star = (c*t)/(R_0**2)
       
        pressure[t_track,:] = np.sum( ( (18*np.square(nu_u-nu) ) / (eta*E) ) * \
                                    ( (np.sin(R_star*np.sqrt(x_n))) / (R_star*np.sin(np.sqrt(x_n)) ) - 1 ) * \
                                    np.exp(-x_n*t_star) , axis=1) 

        # Account for center value
        #pressure[t_track,center] = np.sum( (8*eta*(np.sqrt(x_n) - np.sin(np.sqrt(x_n)))) / ( (x_n - 12*eta + 16*eta*eta)*np.sin(np.sqrt(x_n)) ) * np.exp(-x_n * t_star) )
        pressure[t_track,center] = np.sum( ( (18*np.square(nu_u-nu) ) / (eta*E) ) * \
                                    ( (np.sqrt(x_n)) / (np.sin(np.sqrt(x_n)) ) - 1 ) * \
                                    np.exp(-x_n*t_star)) 
        t_track += 1

    return pressure
    
def displacement_sph(locs, tsteps, x_n):
    """
    Compute displacement field at locations, spherical coordinates.
    """
    (npts, dim) = locs.shape
    ntpts = tsteps.shape[0]
    disp = np.zeros((ntpts, npts), dtype=np.float64)    
    
    center = np.where(~locs.any(axis=1))[0]
    R = np.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
    R_star = R.reshape([R.size,1]) / R_0
    x_n.reshape([1,x_n.size])
    
    E = np.square(1-nu)*np.square(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    
    t_track = 0
    
    for t in tsteps:
        t_star = (c*t)/(R_0**2)
        disp[t_track, :] = R_star.ravel() - np.nan_to_num(np.sum(((12*(1 + nu)*(nu_u - nu)) / \
                                           ((1 - 2*nu)*E*R_star*R_star*x_n*np.sin(np.sqrt(x_n))) ) * \
                                    (3*(nu_u - nu) * (np.sin(R_star*np.sqrt(x_n)) - R_star*np.sqrt(x_n)*np.cos(R_star*np.sqrt(x_n))) + \
                                    (1 - nu)*(1 - 2*nu)*R_star*R_star*R_star*x_n*np.sin(np.sqrt(x_n))) * \
                                    np.exp(-x_n*t_star),axis=1))
                                        
        t_track += 1
    
    return disp

def displacement_crt(locs, tsteps, x_n):
    """
    Compute displacement field at locations, cartesian coordinates.
    """
    (npts, dim) = locs.shape
    ntpts = tsteps.shape[0]
    disp = np.zeros((ntpts, npts, dim), dtype=np.float64)    
    
    center = np.where(~locs.any(axis=1))[0]
    R = np.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
    theta = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( numpy.sqrt(locs[:,0]**2 + locs[:,1]**2) / locs[:,2] ) ) )
    phi = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( locs[:,1] / locs[:,0] ) ) )    
    R_star = R.reshape([R.size,1]) / R_0
    x_n.reshape([1,x_n.size])
    
    E = np.square(1-nu)*np.square(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    
    t_track = 0
    
    for t in tsteps:
        t_star = (c*t)/(R_0**2)
        r_exact_N = R_star.ravel() - np.nan_to_num(np.sum(((12*(1 + nu)*(nu_u - nu)) / \
                                           ((1 - 2*nu)*E*R_star*R_star*x_n*np.sin(np.sqrt(x_n))) ) * \
                                    (3*(nu_u - nu) * (np.sin(R_star*np.sqrt(x_n)) - R_star*np.sqrt(x_n)*np.cos(R_star*np.sqrt(x_n))) + \
                                    (1 - nu)*(1 - 2*nu)*R_star*R_star*R_star*x_n*np.sin(np.sqrt(x_n))) * \
                                    np.exp(-x_n*t_star),axis=1))
                                    
        disp[t_track, :, 0] = (r_exact_N*U_R_inf)*numpy.cos(phi)*numpy.sin(theta)
        disp[t_track, :, 1] = (r_exact_N*U_R_inf)*numpy.sin(phi)*numpy.sin(theta)
        disp[t_track, :, 2] = (r_exact_N*U_R_inf)*numpy.cos(theta)
        t_track += 1    
    return disp

def displacement_small_time(locs, tsteps, x_n):
    """
    Compute displacement field at locations.
    """
    (npts, dim) = locs.shape
    ntpts = tsteps.shape[0]
    disp = np.zeros((ntpts, npts), dtype=np.float64)    
    
    center = np.where(~locs.any(axis=1))[0]
    R = np.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
    R_star = R.reshape([R.size,1]) / R_0
    x_n.reshape([1,x_n.size])
    
    E = np.square(1-nu)*np.square(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    
    t_track = 0
    
    for t in tsteps:
        t_star = (c*t)/(R_0**2)
        disp[t_track, :] = R_star.ravel() + np.nan_to_num( ( 12*(nu_u-nu) ) / ( np.sqrt(np.pi)*(1-nu)*(1+nu_u) )  * \
                                            R_star.ravel()*np.sqrt(t_star) )
                                       
        t_track += 1
    
    return disp

def E_func(x_n, nu, nu_u):
    E = (1-nu)*(1-nu)*(1+nu_u)*(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    return E
    
def stressCalc(locs, tsteps, x_n):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    ntpts = tsteps.shape[0]
    stress = np.zeros((ntpts, npts, 6), dtype=np.float64)    
    
    center = np.where(~locs.any(axis=1))[0]
    R = np.sqrt(locs[:,0]*locs[:,0] + locs[:,1]*locs[:,1] + locs[:,2]*locs[:,2])
    R_star = R.reshape([R.size,1]) / R_0
    x_n.reshape([1,x_n.size])
    theta = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( numpy.sqrt(locs[:,0]**2 + locs[:,1]**2) / locs[:,2] ) ) )
    phi = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( locs[:,1] / locs[:,0] ) ) )
        
    E = np.square(1-nu)*np.square(1+nu_u)*x_n - 18*(1+nu)*(nu_u-nu)*(1-nu_u)
    
    t_track = 0
    
    for t in tsteps:
        t_star = (c*t)/(R_0**2)

        # Normalized Radial Stress
        sigma_RR =  -1.0 - numpy.nan_to_num( numpy.sum( ( (12*(nu_u - nu) ) / \
                                                          (E*R_star*R_star*R_star*x_n*numpy.sin(numpy.sqrt(x_n))) ) * \
                                                        ( 6*(nu_u - nu)*(numpy.sin(R_star*numpy.sqrt(x_n)) - R_star*numpy.sqrt(x_n)*numpy.cos(R_star*numpy.sqrt(x_n))) + \
                                                          -(1 - nu)*(1 + nu_u)*R_star*R_star*R_star*numpy.sin(numpy.sqrt(x_n)) ) * \
                                                        numpy.exp(-x_n*t_star),axis=1))

        # Normalized Circumferential Stress
        sigma_pp =  -1.0 - numpy.nan_to_num( numpy.sum( ( (12*(nu_u - nu)) / \
                                                          (E*R_star*R_star*R_star*x_n*numpy.sin(numpy.sqrt(x_n))) ) * \
                                                        ( 3*(nu_u - nu) * ((R_star*R_star*x_n - 1.0) * (numpy.sin(R_star*numpy.sqrt(x_n))) + R_star*numpy.sqrt(x_n)*numpy.cos(R_star*numpy.sqrt(x_n))) + \
                                                          -(1 - nu)*(1 + nu_u)*R_star*R_star*R_star*x_n*numpy.sin(numpy.sqrt(x_n)) ) * \
                                                        numpy.exp(-x_n*t_star), axis=1))    
    	
        #sigma_thetatheta = sigma_phiphi
        sigma_tt = sigma_pp

        stress[t_track, :, 0] = sigma_RR*P_0        # sxx
        stress[t_track, :, 1] = sigma_tt*P_0        # syy
        stress[t_track, :, 2] = sigma_pp*P_0        # szz
        stress[t_track, :, 3] = numpy.zeros(npts)   # sxy
        stress[t_track, :, 4] = numpy.zeros(npts)   # syz
        stress[t_track, :, 5] = numpy.zeros(npts)   # sxz


        # sigma_sph = numpy.array([ [         sigma_RR*P_0, numpy.zeros(npts), numpy.zeros(npts)],
        #                        	  [numpy.zeros(npts),  sigma_thetatheta*P_0, numpy.zeros(npts)],
        #                           [numpy.zeros(npts), numpy.zeros(npts),      sigma_phiphi*P_0] ])


        # A = numpy.array( [ [numpy.sin(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.cos(phi),   -numpy.sin(phi)],
        #                    [numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)*numpy.sin(phi),    numpy.cos(phi)],
        #                    [               numpy.cos(theta),               -numpy.sin(theta), numpy.zeros(npts)] ] ) 

        # B = numpy.array( [ [numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi),   numpy.cos(theta)],
        #                    [numpy.cos(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.sin(phi),  -numpy.sin(theta)],
        #                    [                -numpy.sin(phi),                  numpy.cos(phi),  numpy.zeros(npts)] ] )    
    
#        sigma_car = numpy.einsum('ij...,jk...->ik...',np.einsum('ij...,jk...->ik...',A,sigma_sph), B)

#        stress[t_track, :, 0] = sigma_car[0,0,:] # sxx
#        stress[t_track, :, 1] = sigma_car[1,1,:] # syy
#        stress[t_track, :, 2] = sigma_car[2,2,:] # szz
#        stress[t_track, :, 3] = sigma_car[0,1,:] # sxy
#        stress[t_track, :, 4] = sigma_car[1,2,:] # syz
#        stress[t_track, :, 5] = sigma_car[0,2,:] # sxz
        t_track += 1

    return stress
    
    
# ==============================================================================

# ==============================================================================
f = h5py.File('./output/step00_hex-poroelastic.h5','r')

t = f['time'][:]
t = t.ravel()

U = f['vertex_fields/displacement'][:]
P = f['vertex_fields/pressure'][:]
E = f['vertex_fields/trace_strain'][:]
C_S = f['vertex_fields/cauchy_stress'][:]

pos = f['geometry/vertices'][:]
R = np.sqrt(pos[:,0]*pos[:,0] + pos[:,1]*pos[:,1] + pos[:,2]*pos[:,2])
# Time steps
# ts = 0.0028666667  # sec
# ts = 0.028666667  # sec
# nts = 2
# tsteps = np.arange(0.0, ts * nts, ts) + ts # sec


# t = f['time'][:]
# t = t.ravel()

# U = f['vertex_fields/displacement'][:]
# P = f['vertex_fields/pressure'][:]
# E = f['vertex_fields/trace_strain'][:]

# pos = f['geometry/vertices'][:]
# R = np.sqrt(pos[:,0]*pos[:,0] + pos[:,1]*pos[:,1] + pos[:,2]*pos[:,2])
# theta_sph = np.nan_to_num( np.arctan( np.nan_to_num( np.sqrt(pos[:,0]**2 + pos[:,1]**2) / pos[:,2] ) ) )
# phi_sph = np.nan_to_num( np.arctan( np.nan_to_num( pos[:,1] / pos[:,0] ) ) )

# # Transform position to spherical coordinates
# pos_sph = np.zeros(pos.shape)
# #U_sph = np.zeros(U.shape)
# # (r, theta, phi)
# pos_sph[:,0] = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
# pos_sph[:,1] = np.nan_to_num( np.arctan( np.nan_to_num( np.sqrt(pos[:,0]**2 + pos[:,1]**2) / pos[:,2] ) ) )
# pos_sph[:,2] = np.nan_to_num( np.arctan( np.nan_to_num( pos[:,1] / pos[:,0] ) ) )

# #pos_sph = np.nan_to_num(pos_sph, nan=0.0)
# #U_sph = np.nan_to_num(U_sph, nan=0.0)

# t_N = (c*t) / R_0**2
# P_N = P / P_0

# U_R = -np.sqrt(U[:,:,0]*U[:,:,0] + U[:,:,1]*U[:,:,1] + U[:,:,2]*U[:,:,2])

# U_R_inf = -( (P_0*R_0*(1-2*nu))/(2*G*(1+nu)))
# U_R_N = U_R/U_R_inf

zeroArray = cryer_zeros_python(nu,nu_u,ITERATIONS)
# P_exact_N = np.reshape(pressure(pos, t, zeroArray),[t.shape[0],pos.shape[0],1])
# U_exact_R_N = np.reshape(displacement(pos, t, zeroArray),[t.shape[0],pos.shape[0],1])
# U_exact_R = U_exact_R_N * U_R_inf
# U_exact_R_N_ST = np.reshdisplacement_sph = displacement_sph.reshape(ntpts,101,101,101)ape(displacement_small_time(pos, t, zeroArray),[t.shape[0],pos.shape[0],1])

# U_exact_R_ST = U_exact_R_N * U_R_zero
# U_exact_N_ST = (U_exact_R_N_ST * U_R_zero) / U_R_inf

# # Cartesian analytical data
# U_exact = np.zeros(U.shape)
# U_exact[:,:,0] = U_exact_R[:,:,0] * np.cos(pos_sph[:,2]) * np.sin(pos_sph[:,1])
# U_exact[:,:,1] = U_exact_R[:,:,0] * np.sin(pos_sph[:,2]) * np.sin(pos_sph[:,1])
# U_exact[:,:,2] = U_exact_R[:,:,0] * np.cos(pos_sph[:,1])


# # Select linear samples

# #x_slice = np.where(~pos[:,1:].any(axis=1))[0]

# x_slice = np.nonzero(np.logical_and(pos[:,1] == 0.0, pos[:,2] == 0.0))[0]
# x_arr1 = R[x_slice]
# x_arr1inds = x_arr1.argsort()
# x_slice = x_slice[x_arr1inds[::1]]

# y_slice = np.nonzero(np.logical_and(pos[:,0] == 0.0, pos[:,2] == 0.0))[0]
# y_arr1 = R[y_slice]
# y_arr1inds = y_arr1.argsort()
# y_slice = y_slice[y_arr1inds[::1]]

# #z_slice = np.where(~pos[:,:2].any(axis=1))[0]

# z_slice = np.nonzero(np.logical_and(pos[:,0] == 0.0, pos[:,1] == 0.0))[0]
# z_arr1 = R[z_slice]
# z_arr1inds = z_arr1.argsort()
# z_slice = z_slice[z_arr1inds[::1]]

# center = np.where(~pos.any(axis=1))[0]


# # Graph time snapshots
# t_steps = t.ravel()
# n_graph_steps = 10
# t_step_array = np.linspace(0,t_steps.size,n_graph_steps).astype(np.int)
# t_step_array[0] += 2
# t_step_array[-1] -= 1
# n_steps = t_N.size


# cm_numeric = ['red','orange','green','blue','indigo', 'violet']
# cm_analytic = ['red','orange','green','blue','indigo', 'violet']

# ==============================================================================
# Stress Calc

# stress_sph = stressCalc(xyz, tsteps, zeroArray).reshape(2,101,101,101,6)
# stress_sph_tensor = np.zeros([2,101,101,101,3,3])
# stress_sph_tensor[:,:,:,:,0,0] = stress_sph[:, :,:,:, 0]    # sxx
# stress_sph_tensor[:,:,:,:,1,1] = stress_sph[:, :,:,:, 1]    # syy
# stress_sph_tensor[:,:,:,:,2,2] = stress_sph[:, :,:,:, 2]    # szz

(npts, dim) = pos.shape
ntpts = tsteps.size
# R_val = np.sqrt(np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2))

C_S_tensor = np.zeros([ntpts,3,3,npts])
C_S_tensor[:,0,0,:] = C_S[:,:,0]
C_S_tensor[:,0,1,:] = C_S[:,:,3]
C_S_tensor[:,0,2,:] = C_S[:,:,5]
C_S_tensor[:,1,0,:] = C_S[:,:,3]
C_S_tensor[:,1,1,:] = C_S[:,:,1]
C_S_tensor[:,1,2,:] = C_S[:,:,4]
C_S_tensor[:,2,0,:] = C_S[:,:,5]
C_S_tensor[:,2,1,:] = C_S[:,:,4]
C_S_tensor[:,2,2,:] = C_S[:,:,2]


stress_sph = stressCalc(pos, tsteps,  zeroArray)
pressure_sph = pressure(pos, tsteps, zeroArray)
disp_sph = displacement_sph(pos, tsteps, zeroArray)
disp_crt = displacement_crt(pos, tsteps, zeroArray)

# for i in np.arange(ntpts):
    # stress_sph[i,R_val > 1.0001,:] = 0

stress_sph_tensor = np.zeros([ntpts,3,3,npts])
stress_sph_tensor[:,0,0,:] = stress_sph[:, :, 0]    # sxx
stress_sph_tensor[:,1,1,:] = stress_sph[:, :, 1]    # syy
stress_sph_tensor[:,2,2,:] = stress_sph[:, :, 2]    # szz
# Shear components are zero 

# Test using einsum
sigma_RR = stress_sph[0,:,0]
sigma_TT = stress_sph[0,:,1]
sigma_PP = stress_sph[0,:,2]

sigma_sph = numpy.array([ [     sigma_RR*P_0, numpy.zeros(npts), numpy.zeros(npts)],
                       	  [numpy.zeros(npts),      sigma_TT*P_0, numpy.zeros(npts)],
                          [numpy.zeros(npts), numpy.zeros(npts),      sigma_PP*P_0] ])

theta = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( numpy.sqrt(pos[:,0]**2 + pos[:,1]**2) / pos[:,2] ) ) )
phi = numpy.nan_to_num( numpy.arctan( numpy.nan_to_num( pos[:,1] / pos[:,0] ) ) )

A = numpy.array( [ [numpy.sin(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.cos(phi),   -numpy.sin(phi)],
                   [numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)*numpy.sin(phi),    numpy.cos(phi)],
                   [               numpy.cos(theta),               -numpy.sin(theta), numpy.zeros(npts)] ] ) 

B = numpy.array( [ [numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi),   numpy.cos(theta)],
                   [numpy.cos(theta)*numpy.cos(phi), numpy.cos(theta)*numpy.sin(phi),  -numpy.sin(theta)],
                   [                -numpy.sin(phi),                  numpy.cos(phi),  numpy.zeros(npts)] ] )    

sigma_car = numpy.einsum('ij...,jk...->ik...',np.einsum('ij...,jk...->ik...',A,sigma_sph), B)


# Test using ensum including time axis
A_tensor = B.T
B_tensor = A.T

A_time = np.zeros([ntpts,3,3,npts])
B_time = np.zeros([ntpts,3,3,npts])

for i in np.arange(ntpts):
    A_time[i,:,:,:] = A[:,:,:]
    B_time[i,:,:,:] = B[:,:,:]

stress_crt_tensor = np.zeros([ntpts,3,3,npts])

stress_crt_tensor[:,:,:,:] = numpy.einsum('hijl,hjkl->hikl',np.einsum('hijl, hjkl->hikl',A_time,stress_sph_tensor[:,:,:,:]), B_time)

# disp_sph[:,R_val > 1.0001] = 0


# disp_sph = disp_sph.reshape(ntpts,101,101,101)

# pressure_sph[:,R_val > 1.0001] = 0
# pressure_sph = pressure_sph.reshape(ntpts,101,101,101)

# stress_sph_R = stress_sph[:,:,0].reshape(ntpts,101,101,101)


# Plot Cartesian Slice
#plt.imshow(stress_crt_tensor[1,0,0,:].reshape([101,101,101])[0,:,:])
#plt.show()

# ====================================================================
# Generate graph of stress

# Select linear samples
t_N = (c*t) / R_0**2

#x_slice = np.where(~pos[:,1:].any(axis=1))[0]

x_slice = np.nonzero(np.logical_and(pos[:,1] == 0.0, pos[:,2] == 0.0))[0]
x_arr1 = R[x_slice]
x_arr1inds = x_arr1.argsort()
x_slice = x_slice[x_arr1inds[::1]]

y_slice = np.nonzero(np.logical_and(pos[:,0] == 0.0, pos[:,2] == 0.0))[0]
y_arr1 = R[y_slice]
y_arr1inds = y_arr1.argsort()
y_slice = y_slice[y_arr1inds[::1]]

#z_slice = np.where(~pos[:,:2].any(axis=1))[0]

z_slice = np.nonzero(np.logical_and(pos[:,0] == 0.0, pos[:,1] == 0.0))[0]
z_arr1 = R[z_slice]
z_arr1inds = z_arr1.argsort()
z_slice = z_slice[z_arr1inds[::1]]

center = np.where(~pos.any(axis=1))[0]


# Define snapshots as per Cheng and Detournay paper
t_N_step_array = np.array([0.001, 0.01, 0.1, 0.5, 1.0, 2.0])
t_step_array_exact = (t_N_step_array*R_0*R_0) / c

t_step_array = np.zeros(t_step_array_exact.size)

for item in np.arange(0,t_step_array_exact.size):
    t_step_array[item] = np.abs(t - t_step_array_exact[item]).argmin()

t_step_array[0] = 0
t_step_array = t_step_array.astype(np.int)    


cm_numeric = ['red','black','green','blue','indigo', 'violet']
cm_analytic = ['red','black','green','blue','indigo', 'violet']

# Pore Pressure at Center
fig, ax = plt.subplots()
fig.set_size_inches(15,10)
ax.semilogx(t_N[:], P_N[:,center,0], color=cm_numeric[0], label='Numerical')
ax.semilogx(t_N[::50], P_exact_N[::50,center,0], color=cm_analytic[0],marker='^', linestyle=' ', label='Analytical')

ax.grid()
ax.legend(loc='lower left')
ax.set(xlabel='Normalized Time, t*', ylabel='Normalized Pressure, P*', title="Cryer's Problem: Normalized Pressure at Center")
fig.tight_layout()
fig.savefig('output/cryer_pressure_at_center_hex.png',dpi = 300)
fig.show()

# ==============================================================================
# Normalized Pore Pressure along x axis
fig, ax = plt.subplots()
fig.set_size_inches(15,10)


tstep = 71
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[-1], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[-1], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
#t0 t_N = 0.01
tstep = 179
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[0], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[0], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )

#t2 t_N = 0.10
tstep = 716 # np.int(n_steps/5 * 1)
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[1], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[1], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )

#t5 t_N = 0.25
tstep = 1791 # np.int(n_steps/5 * 2)
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[2], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[2], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )

#t50 t_N = 0.5
tstep = 3583 # np.int(n_steps/5 * 4)
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[4], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[4], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )

#t70 t_N = 1.0
tstep = 7166
ax.plot(R[x_slice], P_N[tstep, x_slice,0], color=cm_numeric[5], label='Numerical, t* = ' + np.str(np.around(t_N[tstep],3) ) )
ax.plot(R[x_slice], P_exact_N[tstep, x_slice,0], color=cm_analytic[5], marker='^', linestyle=' ',  label='Analytical, t* = ' + np.str(np.around(t_N[tstep],3) ) )


ax.grid()
ax.legend(loc='lower left')
ax.set(xlabel='Normalized Radial Distance, R*', ylabel='Normalized Pressure, P*', title="Cryer's Problem: Normalized Pressure Along Radial Axis")
fig.tight_layout()
fig.savefig('output/cryer_pressure_along_x_axis_norm_hex.png',dpi = 300)
fig.show()
