# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 07:28:16 2022

@author: znoll
"""

import numpy as np

def circumference(xvec,yvec):
    # NOTE: this module and greens_theorem.py do not check for consistent curl of the loop
    if len(xvec) != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if len(xvec) < 3:
        raise RuntimeError('Vectors must form closed shape')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    mysum = 0
    dx = 0
    dy = 0
    for i in range(len(xvec)): 
        # Ensuring loop is closed
        if i == (len(xvec) - 1):
            dx = xvec[i] - xvec[0]
            dy = xvec[i]- yvec[0]
        else:
            dx = xvec[i+1] - xvec[i]
            dy = yvec[i+1] - yvec[i]
        mysum += np.sqrt(dx**2 + dy**2)
    return abs(mysum) # Area is expected to be physical i.e. strictly positive

# x_points = np.array((0,1,2,2,1,0))
# y_points = np.array((-1,-1,0,1,2,2))

# print(f"Circumference is {circumference(x_points,y_points)}")

# Testing code
if __name__ == '__main__':
    input('Press `Enter` to plot funny shape')
    x_points = np.array((0,1,2,2,1,0))
    y_points = np.array((-1,-1,0,1,2,2))
    
    # Zeke Noller idiom
    import matplotlib.pyplot as plt
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='inline') #gui='inline' plots in the IDE and gui='qt5' plots in a separate window
    except:
        print('Unable to plot in desired window')  
        
    _, ax = plt.subplots()
    x_points_new = np.append(x_points,x_points[0])
    y_points_new = np.append(y_points,y_points[0])
    
    ax.plot(x_points_new, y_points_new, label="Closed loop")
    plt.title(f"Circumference is {circumference(x_points,y_points):.2f}")
    plt.show()
    
    input('Press `Enter` to plot circle')
    radius = 5
    x_points = np.linspace(-radius,radius,num=101)
    y_points = np.sqrt(radius**2 - (x_points**2))
    # x_points = np.append(x_points,x_points[::-1])
    # y_points = np.append(y_points,-y_points[::-1])
    print(x_points)
    
    _, ax = plt.subplots()    
    ax.plot(x_points, y_points, label="Circle")
    plt.title(f"Circumference is {circumference(x_points,y_points):.2f}")
    plt.show()
