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
    shell.enable_matplotlib(gui='inline')
except:
    print('Unable to open plotting window')

try:
    from greens_theorem import greens_theorem_area
except:
    raise RuntimeError("Local Green's theorem script could not be imported")
    
span = 9.108 #m
plan_area = 15.221 #m2
chord = plan_area/span

# Loading in aerofoil data 
naca = np.loadtxt("naca0313.txt")
naca_len = len(naca)
half_len = int(naca_len/2) + 1
naca *= chord # Scaling loop from unit length to average chord length
for column in 0, 1: # Reversing direction of points in 
    naca[half_len:,column] = naca[half_len:,column][::-1] 
x_points = naca[:,0]
y_points = naca[:,1]
xsection_area = greens_theorem_area(x_points,y_points)

_, ax = plt.subplots()
ax.plot(naca[:,0], naca[:,1])
plt.title(f"Area = {xsection_area:.3f} m2")
plt.ylim(-1, 1)
plt.show()


