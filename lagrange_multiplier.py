# -*- coding: utf-8 -*-
"""

Script for optimising f(x1,x2,...xn) while subject to constraints g_j(x1,x2,...,xn) >= 0

"""

import numpy as np
import matplotlib.pyplot as plt
try:
    import IPython
    shell = IPython.get_ipython()
    # To plot in Spyder window, uncomment 'inline' line:
    # shell.enable_matplotlib(gui='inline')
    # To plot out of Spyder window, unccoment 'qt5' line:
    shell.enable_matplotlib(gui='qt5')

except:
    print('Unable to open plotting window')  
   
## First and second vectorised partial derivative function ##
def partial(order,function_list,function_index,x_vec,x_index):
    
    # order determines whether this is the first partial derivative order the second partial derivative
    # function_vec is a column vector of functions that take clumn vector as input
    # function_index is the index for which function the partial derivative is being taken on
    # x_vec is 2-column vector containing 
    # Current x (and lambda) values in the first column, current left and right distances in second
    
    foo = function_list[function_index]
    h_now = x_vec[x_index,1]
    
    len_x = np.shape(x_vec)[0]
    x_static = x_vec[:,0].T
    num_central_diff = 5 # The minimum number of points for the 3rd order \
                                      #central difference is 5
    floor_half = int(np.floor(num_central_diff / 2))
    
    x_vals = np.zeros((len_x,num_central_diff))
    for row in range(len_x): #x_i th position
        for column in range(num_central_diff):
            left_shift = column - floor_half
            x_vals[row,column] = x_static[row]
            if row == x_index:
                x_vals[row,column] = x_static[row] + h_now*left_shift
    
    if order == 1:
        partial_diff = (foo(x_vals[:,floor_half+1]) - foo(x_vals[:,floor_half-1]))/ (2*h_now)
        
    elif order == 2: 
        partial_diff = (foo(x_vals[:,floor_half-1]) - 2*foo(x_vals[:,floor_half]) \
                        + foo(x_vals[:,floor_half+1])) / (h_now**2)
    elif order == 3:
        partial_diff = (foo(x_vals[:,floor_half+2]) - 2*foo(x_vals[:,floor_half+1]) \
                        + 2*foo(x_vals[:,floor_half-1]) -foo(x_vals[:,floor_half-2]))\
                        /((2*h_now**3))
    
    return partial_diff

def foo(x_vec):
    x = x_vec[0]
    y = x_vec[1]
    return x**4 + y**4

def goo(x_vec):
    x = x_vec[0]
    y = x_vec[1]
    try:
        function = -2*x + 1
    except:
        x_vec += 1e-3
        x = x_vec[0]
        y = x_vec[1]
        function = -2*x + 1
    return function

### Troubleshooting partial diff function ###
function_list = [foo,goo]

## Plotting ##
leftright = 2
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
x = y = np.linspace(-leftright,leftright,20)
x, y = np.meshgrid(x,y)
z = x**4 + y**4
ax.plot_surface(x, y, z,color='red')
w = -2*x + 1
ax.plot_surface(x,y,w,color='blue')
# ax.axes.set_ylim(bottom=-leftright, top=leftright) 
# ax.axes.set_zlim3d(bottom=-20, top=20)
# ax.axes.set_xlim(bottom=-leftright, top=leftright)

def KKT(function_list,x_init,tol):

    x_vec = x_init
    minimise = function_list[0]
    loop = True
    
    while loop==True:
        # 1. function_list is list contain function names
        # 2. First entry of function_list is the function to be minimised
        # 3. All subsequent entries in function_vec are constraints in form g_j(x1,x2,...)>=0 
        # 4. x_vec is column vector 
        # 5. First column of init_vec contains initial points for numerical analysis
        # 6. Second column of init_vec contains grid distance **MAYBE CHANGE TO BE DYNAMIC**
    
        # equal to number of lambda
        num_constraints = len(function_list) - 1 
        total_num_var= np.shape(x_vec)[0] # rows of x_vec
        num_decision_vars = total_num_var - num_constraints 
        if num_constraints + 1 > num_decision_vars:
            raise RuntimeError('Cannot solve - too many constraints!')  
        
        ## Evaluating function vector ##
        function_vec = np.zeros(total_num_var)
        # Partial derivative wrt decision variables:
        for i in range(num_decision_vars):
            function_vec[i] += partial(1,function_list,0,x_vec,i) # First entry is always partial of \
                                                                  # cost function with respect to x_i
            for j in range(num_constraints):
                function_vec[i] += x_vec[num_decision_vars + j,0] * partial(1,function_list,j+1,x_vec,i)
            #                                                     Subsequent entries are lambda_j * 
            #                                                     partial g_j wrt x_i
                constraint_foo = function_list[j+1]  
                function_vec[j + num_decision_vars] = constraint_foo(x_vec[0:num_decision_vars,0])                                   
            
        print(f"Function vector is {function_vec}")
        
        ## Evaluating J matrix ##
        # Saving memory for J matrix to come
        J= np.zeros((total_num_var,total_num_var))
        # Top left square matrix
        top_left = np.identity((num_decision_vars))
        for i in range(num_decision_vars):
            top_left[i,i] *= partial(2,function_list,0,x_vec,i)
            for j in range(num_constraints):
                top_left[i,i] += -x_vec[num_decision_vars+j,0]\
                                *partial(2,function_list,j+1,x_vec,i)
        # Top right rectangular matrix 
        top_right = np.zeros((num_decision_vars,num_constraints))
        for j in range(num_constraints):
            for i in range(num_decision_vars):
                top_right[i,j] = -partial(1,function_list,j+1,x_vec,i)
        # bottom_left rectangular matrix
        bottom_left = -top_right.T
        bottom_right = np.zeros((num_constraints,num_constraints))
        # bottom_right = np.identity((num_constraints))
        top_matrices = np.hstack((top_left,top_right))
        bottom_matrices = np.hstack((bottom_left,bottom_right))
        J = np.vstack((top_matrices,bottom_matrices))
        print('Jacobian is:')
        print(J)
        
        change = np.linalg.solve(J,function_vec)        
        x_new = x_vec[:,0] - change
        # Find some way to have an integer number of spars or whatever
        
        print(f"The minimal value is {minimise(x_new)}")
        
        print(f"The variable vector is {x_new}")
        
        goo = function_list[1]
        if goo(x_new) > 0:
            print('Safe')
        else: 
            print('Not safe')
        print(f"The constraint function is {goo(x_new)}")
        
        for i in range(num_decision_vars,total_num_var):
            if x_new[i] < 0:
                print('Redundant lambda - optimisation halted')
                loop = False # Find some way to make this nicer
        if np.abs(minimise(x_new)-minimise(x_vec[:,0]))/(minimise(x_vec[:,0])) < tol:
            loop = False
            print('Stable solution found')
            
        print(f"*"*50)
        
        x_vec[:,0] = np.copy(x_new)
        hist_x = np.array(())
        hist_y = np.copy(hist_x)
        ### Jank, please remove once done ###
        if num_decision_vars == 2:
            hist_x = np.append(hist_x,x_vec[0,0])
            hist_y = np.append(hist_y, x_vec[1,0])
            hist_vec = np.array((hist_x[-1],hist_y[-1]))
            if loop == True:
                ax.scatter(hist_x[-1],hist_y[-1],minimise(hist_vec),color='black')
            else: 
                ax.scatter(hist_x[-1],hist_y[-1],minimise(hist_vec),color='green')
 
    # ax.axes.set_ylim3d(bottom=min(hist_y), top=max(hist_y)) 
    # ax.axes.set_zlim3d(bottom=-5, top=5)
    # ax.axes.set_xlim3d(bottom=min(hist_x), top=max(hist_x))
    plt.show()           
    pass
    
# Analogue to defining appropriate mass and stress functions

f_list = [foo,goo]

# x = 1*np.ones((3,2))
x = leftright*np.ones((3,2))
for value in range(len(x)):
    x[value,0] = np.random.random_sample()
    
# x[0,0] = -2
# x[1,0] = np.sqrt(2)
# x[2,0] = 1
x[:,1] = 0.01

# for num in [1,2,3]:
#     print(f"The {num}th partial derivative is {partial(num,f_list,0,x,0)}")
 

KKT(f_list,x,tol=0.005)  
