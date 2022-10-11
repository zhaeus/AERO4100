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
    
try: 
    from circumference import arclength
except:
    raise RuntimeError("Local arc length script could not be imported")





#### Wing skin ####

span = 9.108 #m
semi_span = span/2
num_span = 50
span_space = np.linspace(0,semi_span,num=num_span)

plan_area = 15.221 #m2
chord = plan_area/span
semispan = span/2

### Materials and load conditions ###
airplane_mass = 3620 #kg
gravity = 9.81 #m/s2
force = 2.5*gravity*airplane_mass
force_linear = force/semi_span

def plot_shear_diagram():
    span_vec = np.array((0,semi_span))
    shear_vec = np.array((-force,0))
    _, ax = plt.subplots()
    ax.plot(span_vec,shear_vec,color='blue',label='Shear force')
    ax.set_xlabel('Distance along span [m]')
    ax.set_ylabel('Shear force [N]')
    plt.title('Shear force diagram for pull-up maneouvre')
    # plt.ylim(-0.2,0.2)
    # plt.legend()
    plt.show()

def plot_bendingmoment_diagram():
    num_span = 10
    span_vec = np.linspace(0,semi_span,num=num_span)
    moment_vec = 0.5*force*span_vec**2 - force*semi_span/2
    _, ax = plt.subplots()
    ax.plot(span_vec,moment_vec,color='blue',label='Bending moment')
    ax.set_xlabel('Distance along span [m]')
    ax.set_ylabel('Bending moment [N.m]')
    plt.title('Bending moment diagram for pull-up maneouvre')
    # plt.ylim(-0.2,0.2)
    # plt.legend()
    plt.show()

# Loading in aerofoil data 
naca = np.loadtxt("naca0313.txt")
naca_len = len(naca)
half_len = int(naca_len/2) + 1
naca *= chord # Scaling loop from unit length to average chord length
for column in 0, 1: # Reversing direction of points in 
    naca[half_len:,column] = naca[half_len:,column][::-1] 
    
# Defining inner and outer skin surfaces
x_outer = naca[:,0]
x_outer = np.append(x_outer,x_outer[0])
z_outer = naca[:,1]
z_outer = np.append(z_outer,z_outer[0])
def split(vector_loop): # Fix this shit up
    half_len = int(len(vector_loop))
    upper = vector_loop[:half_len]
    lower = vector_loop[half_len:]
    return upper, lower
x_outer_upper = split(x_outer)[0]
x_outer_lower = split(x_outer)[1]
z_outer_upper = split(z_outer)[0]
z_outer_lower = split(z_outer)[1]
xsection_area = np.abs(greens_theorem_area(x_outer,z_outer))

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
            # m = 0.5*((yvec[i]-yvec[i-1])*(xvec[i+1]-xvec[i]) + (yvec[i+1]-yvec[i])*(xvec[i]-xvec[i-1]))
            mdash = -1/m
        except:
            mdash = 0
        sign = (2*bool(yvec[i]>yvec[i-1])-1)
        if mdash == 0:
            xvec_set[i] = sign*thickness + xvec[i]
            yvec_set[i] = yvec[i]
        xvec_set[i] = np.sqrt((thickness**2)/((mdash**2)+1))*sign + xvec[i]
        yvec_set[i] = (xvec_set[i] - xvec[i])*mdash + yvec[i] 
    set_curve = np.vstack((xvec_set.T,yvec_set.T))
    return set_curve

def plot_aerofoil():
    _, ax = plt.subplots()
    ax.plot(naca[:,0], naca[:,1])
    plt.title(f"Area = {xsection_area:.3f} m2")
    plt.ylim(-1, 1)
    plt.show()
    
skin_thickness = 0.0015 #m
inner_surface = inset_offset(x_outer,z_outer,skin_thickness)
x_inner = inner_surface[0,:]
z_inner = inner_surface[1,:]

def plot_innerouter_surface():
    _, ax = plt.subplots()
    ax.plot(x_inner,z_inner,color='orange',label='Inner')
    ax.plot(x_outer,z_outer,color='blue',label='Outer - NACA0313')
    plt.title('Inner surface and outer surface')
    # plt.ylim(-chord/2, chord/2)
    plt.ylim(-0.2,0.2)
    plt.legend()
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
    
    ax.plot_surface(xx, yy, zz,color='red')
    
    ax.set_title('NACA0313 Airfoil')
    ax.set_zlim(-semispan/3,semispan/3)
    ax.set_xlim(-semispan,semispan)
    plt.show()


# Parametrising skin cross-sectional area

skin_outer_enclosedarea = abs(greens_theorem_area(x_outer,z_outer))
skin_inner_enclosedarea = abs(greens_theorem_area(x_inner,z_inner))
skin_area = abs(skin_outer_enclosedarea - skin_inner_enclosedarea)

def skin_area(skin_thickness_var):
    inner_surface = inset_offset(x_outer,z_outer,skin_thickness_var)
    x_inner_var = inner_surface[0,:]
    z_inner_var = inner_surface[1,:]
    
    skin_outer_enclosedarea_constant = abs(greens_theorem_area(x_outer,z_outer))
    skin_inner_enclosedarea_var = abs(greens_theorem_area(x_inner_var,z_inner_var))
    skin_area = abs(skin_outer_enclosedarea_constant - skin_inner_enclosedarea_var)
    
    return skin_area

skin_thickness_testvar = skin_thickness
# print(f"Skin area with thickness of {skin_thickness_testvar} m is {skin_area(skin_thickness_testvar)}")  


### Box section ###
x_box_left = 0.175
x_box_right = 1.1
box_width = abs(x_box_right - x_box_left)
z_box_top = 0.1
z_box_bottom = -0.05
box_height = abs(z_box_top - z_box_bottom)
x_wingbox_inner = np.array((x_box_left,x_box_right,x_box_right,x_box_left,x_box_left))
z_wingbox_inner = np.array((z_box_top,z_box_top,z_box_bottom,z_box_bottom,z_box_top))

wing_thickness_var = 1.5e-3 #m
web_thickness_var = 1.5e-3 #m
def wingskin_inner(wing_thickness,web_thickness_var):
    x_wingbox_outer = np.array((x_box_left-web_thickness_var,
                                x_box_right+web_thickness_var,
                                x_box_right+web_thickness_var,
                                x_box_left-web_thickness_var,
                                x_box_left-web_thickness_var))
    
    z_wingbox_outer = np.array((z_box_top+wing_thickness_var,
                                z_box_top+wing_thickness_var,
                                z_box_bottom-wing_thickness_var,
                                z_box_bottom-wing_thickness_var,
                                z_box_top+wing_thickness_var))
    return x_wingbox_outer, z_wingbox_outer


    
box_skin_innerenclosed_area = abs(greens_theorem_area(x_wingbox_inner,
                                                      z_wingbox_inner))
box_skin_outerenclosed_area = abs(greens_theorem_area(wingskin_inner(wing_thickness_var,
                                                                     web_thickness_var)[0],
                                                      wingskin_inner(wing_thickness_var,
                                                                     web_thickness_var)[1]))
box_skin_area = abs(box_skin_outerenclosed_area - box_skin_innerenclosed_area)

box_skin_perimeter = 2*(box_width + box_height + 2*wing_thickness_var + 2*web_thickness_var)

### Area, Ixx ... formulas ###
def web_area(web_thickness):
    # Note box_height is constant, web_thickness is variable
    return box_height * web_thickness
def spar_cap_area(cap_width,cap_thickness):
    # Assuming square cap i.e. each side is equal
    return 2*cap_width*cap_thickness - cap_thickness**2

cap_side_var = 30e-3 #m
cap_thickness_var = 2e-3 #m

def rotate(xpoints,ypoints,angle):
    # Expect angle in degrees
    angle = np.radians(angle)
    mat = np.zeros((2,2))
    mat[0,0] = np.cos(angle)
    mat[1,0] = -np.sin(angle)
    mat[0,1] = np.sin(angle)
    mat[1,1] = np.cos(angle)
    points = np.vstack((xpoints,ypoints))
    return mat @ points

def spar_cap_coords(cap_side,cap_thickness,angle):
    # Note rotation is clockwise and about (0,0)
    
    # x_corner = corner_coords[0]
    # y_corner = corner_coords[1]
    
    spar_xpoints = np.array((0,0,cap_thickness,cap_thickness,cap_side,cap_side,0))
    spar_zpoints = np.array((0,cap_side,cap_side,cap_thickness,cap_thickness,0,0))
    
    spar_xpoints, spar_zpoints = rotate(spar_xpoints,spar_zpoints,angle)
    
    # _, ax = plt.subplots()
    # ax.plot(spar_xpoints,spar_zpoints,color='black',label='Spar cap')
    
    # plt.title('Spar cap')
    # # plt.ylim(-chord/2, chord/2)
    # # plt.ylim(-0.2,0.2)
    # plt.legend()
    # plt.show()  
    
    return spar_xpoints, spar_zpoints

def web_spar_I_A(t_web,a_spar,t_spar):
    web_I = 2*(1/12)*t_web*(box_height/2)**3
    sparx = spar_cap_coords(a_spar,t_spar,180)[0]
    sparz = spar_cap_coords(a_spar,t_spar,180)[1]
    spar_area_single = abs(greens_theorem_area(sparx,sparz))
    spar_centroid_z = greens_theorem_centroid(sparx,sparz)[1] + box_height/2
    
    A = spar_cap_area(a_spar,t_spar) + t_web*(box_height**3)/(6*(spar_centroid_z**2))
    I = 2*(web_I + 2*spar_area_single*spar_centroid_z**2)

    return I, A


# Stringers and skin 
stringer_height_var = 30e-3 #m
stringer_width_var = 30e-3 #m
stringer_thickness_var = 2e-3 #m

def plot_box2d():
    _, ax = plt.subplots()
    # Aerofoil
    ax.plot(x_inner,z_inner,color='orange',label='Inner')
    ax.plot(x_outer,z_outer,color='blue',label='Outer mold line - NACA0313')
    
    # Skin and web
    # ax.plot(x_wingbox_inner,
    #         z_wingbox_inner,
    #         color='black',label='Wing box inner')
    ax.plot(wingskin_inner(wing_thickness_var,web_thickness_var)[0],
            wingskin_inner(wing_thickness_var,web_thickness_var)[1],
            color='black',label='Wing box inner surface')
    
    # Spar caps (corners)
    for index in range(4):
        angle = (index+1)*90 
        ax.plot(wingskin_inner(wing_thickness_var,web_thickness_var)[0][index] 
                + spar_cap_coords(cap_side_var,cap_thickness_var,angle)[0],
                wingskin_inner(wing_thickness_var,web_thickness_var)[1][index] 
                + spar_cap_coords(cap_side_var,cap_thickness_var,angle)[1],
                color='black',linewidth=5)
    
    
    plt.title('Wing cross section')
    # plt.ylim(-chord/2, chord/2)
    plt.ylim(-0.2,0.2)
    plt.legend()
    plt.show()

def plot_box3d():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    
    span_num = 20
    y = np.linspace(0,semispan,span_num)
    
    # Skin
    skinx = x_outer_upper 
    skinxx, skinyy = np.meshgrid(skinx,y)
    skinz = z_outer
    skinzz,_ = np.meshgrid(skinz,y)
    ax.plot_surface(skinxx, skinyy, skinzz,color='red')
    
    # Box
    box_x = x_wingbox
    box_xx, box_yy = np.meshgrid(box_x,y)
    box_z = z_wingbox
    box_zz,_ = np.meshgrid(box_z,y)
    ax.plot_surface(box_xx, box_yy, box_zz,color='black')
    
    ax.set_title('NACA0313 Airfoil')
    ax.set_zlim(-semispan/5,semispan/5)
    ax.set_xlim(-chord/2,1.5*chord)
    plt.show()

    