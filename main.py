# -*- coding: utf-8 -*-
"""
AERO4100 Design Assignment 

Zecchaeus Noller 
"""

# Loading in modules 
import numpy as np
import matplotlib.pyplot as plt

try:
    import IPython
    shell = IPython.get_ipython()
    shell.enable_matplotlib(gui='qt5') #qt5, inline
except:
    print('Unable to open plotting window')

try:
    from greens_theorem import greens_theorem_area, greens_theorem_centroid
except:
    raise RuntimeError("Local Green's theorem script could not be imported")

span = 9.108 #m
plan_area = 15.221 #m2
chord = plan_area/span
semispan = span/2

# Skin

# Loading in aerofoil data 
naca = np.loadtxt("naca0313.txt")
naca_len = len(naca)
half_len = int(naca_len/2) + 1
naca *= chord # Scaling loop from unit length to average chord length
for column in 0, 1: # Reversing direction of points in 
    naca[half_len:,column] = naca[half_len:,column][::-1] 
    
# Defining inner and outer skin surfaces
x_outer = naca[:,0]
z_outer = naca[:,1]
def split(vector_loop):
    half_len = int(len(vector_loop))
    upper = vector_loop[:half_len]
    lower = vector_loop[half_len:]
    return upper, lower
x_outer_upper = split(x_outer)[0]
x_outer_lower = split(x_outer)[1]
z_outer_upper = split(z_outer)[0]
z_outer_lower = split(z_outer)[1]
xsection_area = np.abs(greens_theorem_area(x_outer,z_outer))

def plot_aerofoil():
    _, ax = plt.subplots()
    ax.plot(naca[:,0], naca[:,1])
    plt.title(f"Area = {xsection_area:.3f} m2")
    plt.ylim(-1, 1)
    plt.show()

def inset_offset(xvec,yvec,thickness):
    lenx = len(xvec)
    leny = len(yvec)
    if lenx != leny:
        raise RuntimeError('Vectors must have same length')
    if lenx < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    xvec_set = np.zeros(lenx+1)
    yvec_set = np.zeros(leny+1)
    for i in range(1,lenx-1):
        try:
            m = (yvec[i+1]-yvec[i-1])/(xvec[i+1]-xvec[i-1])
            mdash = -1/m
        except:
            mdash = 0
        xvec_set[i] = np.sqrt((thickness**2)/((mdash**2)+1)) + xvec[i]
        yvec_set[i] = ((abs(yvec[i]))/(yvec[i])) * (xvec_set[i] - xvec[i])*abs(mdash) + yvec[i]
    set_curve = np.vstack((xvec_set.T,yvec_set.T))
    return set_curve

wing_thickness = 0.005 #m
inner_surface = inset_offset(x_outer,z_outer,wing_thickness)

def plot_innerouter_surface():
    _, ax = plt.subplots()
    ax.plot(inner_surface[0,:], inner_surface[1,:],color='orange')
    ax.plot(naca[:,0], naca[:,1],color='blue')
    plt.title('Inner surface and outer surface')
    plt.ylim(-0.2, 0.2)
    plt.show()

# Method 1: Surface
def plot_aerofoil_surface():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    
    x = x_outer_upper 
    
    span_num = 10
    y = np.linspace(0,semispan,span_num)
    
    xx, yy = np.meshgrid(x,y)
    
    z = z_outer_upper
    zz,_ = np.meshgrid(z,y)
    
    ax.plot_surface(xx, yy, zz)
    
    ax.set_title('NACA0313 Airfoil')
    ax.set_zlim(-semispan,semispan)
    ax.set_xlim(-semispan,semispan)
    plt.show()
    
    