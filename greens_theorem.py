# -*- coding: utf-8 -*-
"""
Code for calculating area of closed shape with discrete points 

Zecchaeus Noller 

September 2022
"""

import numpy as np

def greens_theorem_area(xvec,yvec):
    if len(xvec) != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if len(xvec) < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    mysum = 0
    for i in range(len(xvec)): 
        if i == (len(xvec) - 1):
            # Ensuring loop is closed
            dy = yvec[0] - yvec[i] 
            dx = xvec[0] - xvec[i] 
        else:
            dy = yvec[i+1] - yvec[i]
            dx = xvec[i+1] - xvec[i]
            mysum += dy*(xvec[i] + 0.5*dx) # Derived sum equation for Green's theorem with uniform force field
    return mysum # Area is expected to be physical i.e. strictly positive

def greens_theorem_sma(xvec,yvec):
    if len(xvec) != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if len(xvec) < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    mysum = 0
    
def greens_theorem_centroid(xvec,yvec):
    if len(xvec) != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if len(xvec) < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    area = greens_theorem_area(xvec,yvec)
    my_xsum = my_ysum = 0
    for i in range(1,len(xvec)-1):
        my_xsum += (xvec[i+1]+xvec[i])*(xvec[i]*yvec[i+1] - xvec[i+1]*yvec[i])
        my_ysum += (yvec[i+1]+yvec[i])*(xvec[i]*yvec[i+1] - xvec[i+1]*yvec[i])
    my_xsum /= 6*area
    my_ysum /= 6*area
    return my_xsum, my_ysum
  
# # Testing code
# if __name__ == '__main__':
#     input('Press `Enter` to plot funny shape')
#     x_points = np.array((0,1,2,2,1,0))
#     y_points = np.array((-1,-1,0,1,2,2))
    
#     import matplotlib.pyplot as plt
#     try:
#         import IPython
#         shell = IPython.get_ipython()
#         shell.enable_matplotlib(gui='inline') #to plot in different window change 'inline' to 'qt5'
#     except:
#         print('Unable to open plotting window')  
        
#     _, ax = plt.subplots()
#     x_points_new = np.append(x_points,x_points[0])
#     y_points_new = np.append(y_points,y_points[0])
    
#     ax.plot(x_points_new, y_points_new, label="Closed loop")
#     plt.title(f"Area = {greens_theorem_area(x_points,y_points)}")
#     plt.show()
    
#     input('Press `Enter` to plot circle')
#     radius = 5
#     x_points = np.linspace(-5,5,50)
#     y_points = np.sqrt(radius**2 - x_points**2)
#     x_points = np.append(x_points,x_points[::-1])
#     y_points = np.append(y_points,-y_points[::-1])
    
#     _, ax = plt.subplots()    
#     ax.plot(x_points, y_points, label="Circle")
#     plt.title(f"Area = {greens_theorem_area(x_points,y_points)}")
#     plt.show()
