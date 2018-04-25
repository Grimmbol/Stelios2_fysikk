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

#TODO plot:
# - The total energy and the regression
# - The actual path of the ball and the numeric estimation
# - The y over time plot for both the numeric and the observed path
# - The normal force over time, plotted for both numeric and observed
# - The friction force over time, for both numeric and observed
# - The air drag over time, for both numeric and observed

def main():
    #First, generate the data needed. These methods are defined in data_generator.py
    avg_data = generate_avg_data("Tracker", 36)
    #num_res = generate_avg_num_data("Tracker", 3)
    set_global(avg_data[2][0])
    
    fitted = exp_fit_avg("energidata", num_data_files)
    
    #plt.plot(avg_data[0], avg_data[2])
    
    plot_energy_fit_with_scatter(fitted)
    #plot_height_with_scatter(fitted)
    #Finally, show the plots
    plt.savefig("Lin_Reg_W(t).svg")
    plt.show()
    
main()
