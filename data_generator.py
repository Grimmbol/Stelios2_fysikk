from iptrack import *
from scipy.optimize import curve_fit

import numpy as np

#Numerical constants
euler_step_len = 0.0005
fric_const = 0.00145

glob_C = 1

#Constants relating to files
tracker_path = "Tracker"
data_loc = "energidata"
num_data_files = 35
max_duration = 11 #in seconds

#Physical constant values
mass_ball = 0.0029 #+- 0.05 grams 
grav_constant = 9.8665
radius_ball = 0.0375 #-+ 0.05 mm 
#The ball is modeled with an infinitly thin shell
inertia_ball = (3/2)*(mass_ball)*radius_ball**2
spacing_track_supports = 0.20 #+- 0.3

def set_global(C):
    global glob_C
    glob_C = C

#The exponential function we are trying to fit.
def exp_func(t,C,lam):
    return C*np.exp(lam*t)

#Try to use this for fit, unneccesary to use least squares to find C
def exp_func_const(t, lam):
    return glob_C*np.exp(lam*t)

#Returns a list of tuples containing the error of each parameter
def find_error(avg, folder, number):
    
    
    for i in range(1, number+1):
        data = np.loadtxt(folder+"/v"+str(i)+"_track", skiprows=2)
        
    
    
#Returns a list with the normal force for each data point in the given data
def generate_normal_force(data):
    polynomial = np.polyfit(data[1], data[2], 15)

    normal_series = []
    for lineno in range(len(data[0])):
        y, dydx, d2ydx2, alpha, R = trvalues(polynomial, data[1][lineno])
        normal_force = (((mass_ball*data[3][lineno]**2) / R)
                       +(mass_ball*grav_constant*np.cos(alpha)))
        normal_series.append(normal_force)
    return normal_series

#Returns a list of the friction force giving the ball rotation per point in the
#input data
def generate_friction_force(data):
    polynomial = np.polyfit(data[1], data[2], 15)

    friction_series = []
    for lineno in range(len(data[0])):
        y, dydx, d2ydx2, alpha, R = trvalues(polynomial, data[1][lineno])
        friction = ((mass_ball*grav_constant*np.sin(alpha))
                    -(mass_ball
                      *((mass_ball*grav_constant*np.sin(alpha)-fric_const*data[3][lineno])
                        /(mass_ball+(inertia_ball/radius_ball**2)))))
        friction_series.append(friction)
    return friction_series

#Iterates over the chosen number of files in increasing order, taking the
#average as it works its way through
def generate_avg_data(folder, number):
   
    #A list of datasets
    temp = []
    count_list = []
    for i in range(1, number+1):
        current_set = ([],[],[],[])
        data = np.loadtxt(folder+"/v"+str(i)+"_track", skiprows=2)
        line_counter = 0
        for line in data:
            if line_counter >= len(count_list):
                count_list.append(1)
            else:
                count_list[line_counter] += 1
            
            try:
                vel = (np.sqrt(
                              ((line[1]-current_set[1][line_counter-1])**2)+
                              ((line[2]-current_set[2][line_counter-1])**2))
                       /(line[0]-current_set[0][line_counter-1]))
            except:
                vel = 0
            try:
                current_set[0][line_counter] += line[0]
                current_set[1][line_counter] += line[1]
                current_set[2][line_counter] += line[2]
                current_set[3][line_counter] += vel

            except:
                current_set[0].append(line[0])
                current_set[1].append(line[1])
                current_set[2].append(line[2])
                current_set[3].append(vel)
        
            line_counter += 1
        temp.append(current_set)
    print(count_list)
    
    #Stores the final sequence, properly averaged
    result = ([0 for i in range(len(count_list))],
              [0 for i in range(len(count_list))],
              [0 for i in range(len(count_list))],
              [0 for i in range(len(count_list))])

    #Take average
    #Usual format: t, x, y, v
    for set in temp:
       for i in range(len(set[0])):
           result[0][i] += set[0][i]
           result[1][i] += set[1][i]
           result[2][i] += set[2][i]
           result[3][i] += set[3][i]
       

    for i in range(len(count_list)):
        result[0][i] /= count_list[i]
        result[1][i] /= count_list[i]
        result[2][i] /= count_list[i]
        result[3][i] /= count_list[i]
    
    return(result)

#Generates numeric data using eulers explicit method. Averages all results
def generate_numeric_data(folder, number):

    #A tuple to store the numeric results. Format (t, x, y, v)
    num_res = ([],[],[],[])
   
    #A loop doing the euler estimation, and averaging the results
    for i in range(1, number+1):
        new_num_res = ([],[],[],[])
        tdata = np.linspace(0,max_duration,max_duration/euler_step_len)
        tracker_data = np.loadtxt(folder+"/v"+str(i)+"_track", skiprows=2)

        polynomial = iptrack(folder+"/v"+str(i)+"_track", 15)

        #Set up our start positions
        cur_x = tracker_data[0][1]
        cur_y = tracker_data[0][2]
        cur_v = 0
        
        for point in tdata:
            y, dydx, d2ydx2, alpha, R = trvalues(polynomial, cur_x)
            
            new_num_res[0].append(point)
            new_num_res[1].append(cur_x)
            new_num_res[2].append(cur_y)
            new_num_res[3].append(cur_v)
            
            #Update position values
            cur_x = cur_x + np.cos(alpha)*cur_v*euler_step_len
            cur_y = cur_y - np.sin(alpha)*cur_v*euler_step_len
            
            #Estimate next velocity (always paralell with the track)
            cur_v = (cur_v 
                  + (mass_ball*grav_constant*np.sin(alpha)-(fric_const*cur_v))
                  / (mass_ball + (inertia_ball/(radius_ball**2)))
                  *  euler_step_len)
        
        #Average out the new result with the old ones
        if(len(num_res[0]) == 0):
            num_res = new_num_res
        else:
            diff = len(new_num_res[0])-len(num_res[0])
            for i in range(min(len(new_num_res[0]), len(num_res[0]))):
                num_res[0][i] += new_num_res[0][i]#;num_res[0][i] /= 2
                num_res[1][i] += new_num_res[1][i]#;num_res[1][i] /= 2
                num_res[2][i] += new_num_res[2][i]#;num_res[2][i] /= 2
                num_res[3][i] += new_num_res[3][i]#;num_res[3][i] /= 2
            
            
    return num_res

def exp_fit_avg(folder, number):
    lam, err = (0, 0)
    global glob_C
    for i in range(1, number+1):
        data = np.loadtxt("./" + str(folder)+"/data_"+str(i),skiprows=1)
        ydata = []
        tdata = [] 

        for line in data:
            ydata.append(line[1])
            tdata.append(line[0])
        
        result = curve_fit(exp_func_const, xdata=tdata, ydata=ydata)
        #Update global constant
        lam += result[0][0]
        err += np.sqrt(np.diag(result[1]))

    #Taking average of params
    lam = lam/number
    err = err/number
    
    #Neatly package our results
    result_vec = (lam, err)

    return result_vec
