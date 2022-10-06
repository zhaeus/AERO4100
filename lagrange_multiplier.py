# -*- coding: utf-8 -*-
"""

Script for optimising f(x1,x2,...xn) while subject to constraint g(x1,x2,...,xn) = 0

"""

import numpy as np
import matplotlib.pyplot as plt


### Troubleshooting and plotting ### 
# def plotting():
#     # Idiom for plotting outside of Spyder window 
#     import matplotlib.pyplot as plt
#     try:
#         import IPython
#         shell = IPython.get_ipython()
#         shell.enable_matplotlib(gui='qt5') #qt5, inline
#     except:
#         print('Unable to open plotting window')
        
#     num_x = num_y = 20
#     x = np.linspace(-2,2,num=num_x)
#     y = np.linspace(-2,2,num=num_y)
#     x,y = np.meshgrid(x,y)
    
#     minimise = x**3 - x**2 - 2*x - y**3 + 1      
    
    
def partial(order,function_list,function_index,x_vec,x_index):
    
    # order determines whether this is the first partial derivative order the second partial derivative
    # function_vec is a column vector of functions that take clumn vector as input
    # function_index is the index for which function the partial derivative is being taken on
    # x_vec is 2-column vector containing 
    # Current x (and lambda) values in the first column, current left and right distances in second
    
    foo = function_list[function_index]
    x_now = x_vec[x_index,0] 
    h_now = x_vec[x_index,1]
    
    x_left = np.copy(x_vec)
    x_left[x_index] -= h_now
    x_mid = np.copy(x_vec)
    x_mid[x_index] += 0
    x_right = np.copy(x_vec)
    x_right[x_index] += h_now
    
    if order == 1:
        partial_diff = (foo(x_right) - foo(x_left))/ (2*h_now)
        
    elif order == 2: 
        partial_diff = (foo(x_right) - 2*foo(x_mid) + foo(x_left)) / (h_now**2)
    
    return partial_diff[0]
    

def KKT(function_list,x_init):
    num_steps = 10
    iteration = 0 
    x_vec = x_init
    while iteration<num_steps:
        iteration += 1
        # 1. function_list is list contain function names
        # 2. First entry of function_list is the function to be minimised
        # 3. All subsequent entries in function_vec are constraints in form g_j(x1,x2,...)>=0 
        # 4. x_vec is column vector 
        # 5. First column of init_vec contains initial points for numerical analysis
        # 6. Second column of init_vec contains grid distance **MAYBE CHANGE TO BE DYNAMIC**
    
        num_constraints = len(function_list)-1  # equal to number of lambda
        total_num_var= np.shape(x_vec)[0] # rows of x_vec
        num_decision_vars = total_num_var - num_constraints 
        if num_constraints > num_decision_vars:
            raise RuntimeError('Cannot solve; too many constraints!')  
        
        ## Evaluating function vector ##
        function_vec = np.zeros(total_num_var)
        # Partial derivative wrt decision variables:
        for i in range(num_decision_vars):
            function_vec[i] += partial(1,function_list,0,x_vec,i) # First entry is always partial of \
                                                                  # cost function with respect to x_i
            for j in range(1,len(function_list)):
                function_vec[i] += -x_vec[num_decision_vars + j - 1][0] * partial(1,function_list,j,x_vec,i)
            #                                                     Subsequent entries are lambda_j * 
            #                                                     partial g_j wrt x_i
        for i in range(num_constraints):
            constraint_foo = function_list[i+1]
            function_vec[i + num_decision_vars] += constraint_foo(x_vec[0])
        
        ## Evaluating J matrix ##
        # Saving memory for J matrix to come
        J= np.zeros((total_num_var,total_num_var))
        # Top left square matrix
        top_left = np.identity((num_decision_vars))
        for i in range(num_decision_vars):
            top_left[i,i] *= partial(2,function_list,0,x_vec,i)
        # Top right rectangular matrix 
        top_right = np.zeros((num_decision_vars,num_constraints))
        for j in range(num_constraints):
            for i in range(num_decision_vars):
                top_right[i,j] = partial(1,function_list,j+1,x_vec,i)
        # bottom_left rectangular matrix
        bottom_left = -top_right.T 
        bottom_right = np.zeros((num_constraints,num_constraints))
        top_matrices = np.hstack((top_left,top_right))
        bottom_matrices = np.hstack((bottom_left,bottom_right))
        J = np.vstack((top_matrices,bottom_matrices))
        # print(np.size(J))
        # print(np.size(x_vec[:,0]))
        
        change = np.linalg.solve(J,function_vec)
        # print(np.shape(change))
        x_vec[:,0] -= change 
        print(x_vec[:,0])
        
        
        ### very jank, please delete ###
        if num_decision_vars == 2 and iteration == 1:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        
            x = y = np.linspace(-10,10,20)
            xx, yy = np.meshgrid(x,y)
            
            zz = xx**2+yy**2
            ax.plot_surface(xx, yy, zz,color='red')
            
            ww = -xx-yy + 3
            ww = ax.plot_surface(xx,yy,ww,color='blue')
            
            ax.set_xlim(-10,10)
            ax.set_ylim(-10,10)
            plt.show()
pass
    
    
def foo(x_vec):
    x = x_vec[0]
    y = x_vec[1]
    return x**2 + y**2

def goo(x_vec):
    x = x_vec[0]
    y = x_vec[1]
    return -x -y + 3

f_list = []
f_list.append(foo)
f_list.append(goo)

x = 1*np.ones((3,2))
x[:,1] = 0.1
# print(x)
# print(partial(2,f_list,0,x,1))
    
KKT(f_list,x)  
    
    
    
    
"""

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(x,y,foo,color='blue')

ax.set_title('Function surface')
# ax.set_zlim(-lim,lim)
# ax.set_xlim(-lim,lim)
plt.show() """