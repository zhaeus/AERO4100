# -*- coding: utf-8 -*-
"""
Code for calculating area of closed shape with discrete points 

Zecchaeus Noller 

September 2022

All of these formulae are adapted from:
https://leancrew.com/all-this/2018/01/greens-theorem-and-section-properties/ 

"""

import numpy as np

def green_area(xvec,yvec):
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
    return mysum # Note that area doesn't have to be positive for use in subsequent formulae 
    
def green_centroid(xvec,yvec):
    if len(xvec) != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if len(xvec) < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    area = green_area(xvec,yvec)
    my_xsum = my_ysum = 0
    for i in range(0,len(xvec)-1):
        my_xsum += (xvec[i+1]+xvec[i])*(xvec[i]*yvec[i+1] - xvec[i+1]*yvec[i])
        my_ysum += (yvec[i+1]+yvec[i])*(xvec[i]*yvec[i+1] - xvec[i+1]*yvec[i])
    my_xsum /= 6*area
    my_ysum /= 6*area
    return my_xsum, my_ysum

def green_sma(xvec,yvec):
    # Gives Ixx, Iyy, Ixy relative to coordinate system from initial point in canonical directions
    lenx = len(xvec)
    if lenx != len(yvec):
        raise RuntimeError('Vectors must have same length')
    if lenx < 3:
        raise RuntimeError('Vectors must form closed polygon')
    if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
        np.append(xvec,xvec[0])
        np.append(yvec,yvec[0])
    lenx = len(xvec)
    sumxx = sumyy = sumxy = 0
    A = green_area(xvec,yvec)
    cx, cy = green_centroid(xvec,yvec)
    for i in range(lenx-1):
        sumxx += (yvec[i]**2 + yvec[i]*yvec[i+1] + yvec[i+1]**2)*(xvec[i]*yvec[i+1] 
                                                                  - xvec[i+1]*yvec[i])
        sumyy += (xvec[i]**2 + xvec[i]*xvec[i+1] + xvec[i+1]**2)*(xvec[i]*yvec[i+1]  
                                                                  - xvec[i+1]*yvec[i])
        sumxy += (xvec[i]*yvec[i+1] + 2*xvec[i]*yvec[i] 
                  + 2*xvec[i+1]*yvec[i+1] + xvec[i+1]*yvec[i])*(xvec[i]*yvec[i+1] - xvec[i+1]*yvec[i])
    Ixx = abs(sumxx/12 - A*cy**2)
    Iyy = abs(sumyy/12 - A*cx**2)
    Ixy = abs(sumxy/24 - A*cx*cy)
    return Ixx, Iyy, Ixy

# Testing code
if __name__ == '__main__':
    input('Press `Enter` to plot funny shape')
    x_points = np.array((0,1,2,2,1,0))
    y_points = np.array((-1,-1,0,1,2,2))
    
    import matplotlib.pyplot as plt
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='inline') #to plot in different window change 'inline' to 'qt5'
    except:
        print('Unable to open plotting window')  
        
    _, ax = plt.subplots()
    x_points_new = np.append(x_points,x_points[0])
    y_points_new = np.append(y_points,y_points[0])
    
    ax.plot(x_points_new, y_points_new, label="Closed loop")
    plt.title(f"Area = {green_area(x_points,y_points)}, \
              Centroid = {green_centroid(x_points,y_points)}, \
              Ixx = {green_sma(x_points_new,y_points_new)[0]},\
              Iyy = {green_sma(x_points_new,y_points_new)[1]}")
    plt.show()
    
    input('Press `Enter` to plot circle')
    radius = 5
    x_points = np.linspace(-5,5,50)
    y_points = np.sqrt(radius**2 - x_points**2)
    x_points = np.append(x_points,x_points[::-1])
    y_points = np.append(y_points,-y_points[::-1])
    
    _, ax = plt.subplots()    
    ax.plot(x_points, y_points, label="Circle")
    plt.title(f"Area = {green_area(x_points,y_points)}")
    
    input('Press `Enter` to plot rectangle')
    x_points = np.array((-1,-1,1,1,-1))
    y_points = np.array((0,4.5,4.5,0,0))
    
    import matplotlib.pyplot as plt
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='inline') #to plot in different window change 'inline' to 'qt5'
    except:
        print('Unable to open plotting window')  
        
    _, ax = plt.subplots()
    # x_points_new = np.append(x_points,x_points[0])
    # y_points_new = np.append(y_points,y_points[0])
    
    ax.plot(x_points, y_points, label="Closed loop")
    plt.title(f"Area = {abs(green_area(x_points,y_points))}, \
              Centroid = {green_centroid(x_points,y_points)}, \
              Ixx = {green_sma(x_points,y_points)[0]},\
              Iyy = {green_sma(x_points,y_points)[1]}")
    plt.show()
    
