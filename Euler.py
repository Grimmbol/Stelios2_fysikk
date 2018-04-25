# -*- coding: utf-8 -*-

from data_generator import *
import math
from iptrack import *
import matplotlib.pyplot as plt
import numpy as np

def polynomial(coefs, x):
    result = 0
    
    exp_counter = len(coefs)-1
    print(exp_counter, len(coefs))
    for i in range(len(coefs)):
        result += coefs[i] * (x**exp_counter)
        exp_counter -= 1

    return result

def plot_energy_fit_with_scatter(fitted):
    #Draw the regression function
    t = np.linspace(0, max_duration, 500)
    y = potential_energy(exp_func_const(t, fitted[0]))
    
    plt.plot(t, y, label="regression")
    
    #Generate scatterplot of all data,
    for i in range(1, num_data_files+1):
        data = np.loadtxt("./" + data_loc +"/data_"+str(i),skiprows=1)
        ydata = []
        tdata = []

        for line in data:
            ydata.append(potential_energy(line[1]))
            tdata.append(line[0])

        plt.scatter(tdata, ydata)

    #Print the standard error
    print("Standard error for lambda is " + str(np.sqrt(np.diag(fitted[1]))))
    print("Lambda is " + str(fitted[0]))

def plot_height_with_scatter(fitted):
     #Draw the regression function
    t = np.linspace(0, max_duration, 500)
    y = exp_func_const(t, fitted[0])

    plt.plot(t, y, label="regression")
    
    #Generate scatterplot of all data,
    for i in range(1, num_data_files+1):
        data = np.loadtxt("./" + data_loc +"/data_"+str(i),skiprows=1)
        ydata = []
        tdata = []

        for line in data:
            ydata.append(line[1])
            tdata.append(line[0])

        plt.scatter(tdata, ydata)

def plot_y_fit(fitted):
     t = np.linspace(0,10,500)
     y = exp_func(t, fitted[0], fitted[1])
     
     fig = plt.figure()
     fig.plot(t, y, label="regression")

def plot_normal_force(data):
    normal_series = generate_normal_force(data)
    
    #plt.plot(data[0], data[2])    
    plt.plot(data[0], normal_series)

def plot_friction_force(data):
    friction_series = generate_friction_force(data)

    #plt.plot(data[0], data[2])
    plt.plot(data[0], friction_series)

def plot_error_y(data, errors):
    over = ([],[],[],[])
    under = ([],[],[],[])
    for i in range(len(data[0])):
        over[0].append(data[0][i])
        over[1].append(data[1][i] + errors[0][i])
        over[2].append(data[2][i] + errors[1][i])
        over[3].append(data[3][i] + errors[2][i])

        under[0].append(data[0][i])
        under[1].append(data[1][i] - errors[0][i])
        under[2].append(data[2][i] - errors[1][i])
        under[3].append(data[3][i] - errors[2][i])

    plt.plot(data[0], data[2])
    plt.plot(data[0], over[2])
    plt.plot(data[0], under[2])

#TODO plot:
# - The total energy and the regression
# - The actual path of the ball and the numeric estimation
# - The y over time plot for both the numeric and the observed path
# - The normal force over time, plotted for both numeric and observed
# - The friction force over time, for both numeric and observed
# - The air drag over time, for both numeric and observed

def main():
    #First, generate the data needed. These methods are defined in data_generator.py
    avg_data, errors = generate_avg_data("Tracker", 36)
    #num_res = generate_avg_num_data("Tracker", 36)
    
    #print(errors)
    
    print(unc_W(500))
   

    set_global(avg_data[2][0])
    print("C is " + str(avg_data[2][0]) + " with standard error " + str(errors[1][0]))
    fitted = exp_fit_avg("energidata", num_data_files)
    
    #avg_data_normal = generate_normal_force(avg_data)
    #avg_num_data_normal = generate_normal_force(num_res)

    #avg_data_friction = generate_friction_force(avg_data)
    #avg_num_data_friction  = generate_friction_force(num_res)

    #plt.plot(num_res[0], num_res[2])
    #plt.plot(num_res[0], avg_num_data_normal)
    #plt.plot(num_res[0], avg_num_data_friction)

    #plot_error_y(avg_data, errors)
    #print(errors[3])
    
    plot_energy_fit_with_scatter(fitted)
    #plot_height_with_scatter(fitted)
    #Finally, show the plots
    plt.show()
    
main()
