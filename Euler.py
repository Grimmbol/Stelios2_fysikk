# -*- coding: utf-8 -*-

from data_generator import *
import math
from iptrack import *
import matplotlib.pyplot as plt
import numpy as np

glob_C = 3.94

def polynomial(coefs, x):
    result = 0
    
    exp_counter = len(coefs)-1
    print(exp_counter, len(coefs))
    for i in range(len(coefs)):
        result += coefs[i] * (x**exp_counter)
        exp_counter -= 1

    return result

def potential_energy(h):
    return mass_ball * grav_constant * h

def plot_energy_fit_with_scatter(fitted):
    #Draw the regression function
    t = np.linspace(0, 10, 500)
    y = potential_energy(exp_func(t, fitted[0], fitted[1]) )

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

def plot_y_fit(fitted):
     t = np.linspace(0,10,500)
     y = exp_func(t, fitted[0], fitted[1])

     plt.plot(t, y, label="regression")

def plot_normal_force(data):
    normal_series = generate_normal_force(data)

    #plt.plot(data[0], data[2])    
    plt.plot(data[1], normal_series)

def plot_friction_force(data):
    friction_series = generate_friction_force(data)

    #plt.plot(data[0], data[2])
    plt.plot(data[1], friction_series)


def plot_friction(track_data, euler_data):
    pass

#TODO plot:
# - The total energy and the regression
# - The actual path of the ball and the numeric estimation
# - The y over time plot for both the numeric and the observed path
# - The normal force over time, plotted for both numeric and observed
# - The friction force over time, for both numeric and observed
# - The air drag over time, for both numeric and observed
# Make the regression only depend on lambda, with the constant set by the tracker
# data

def main():
    #First, generate the data needed. These methods are defined in data_generator.py
    avg_data = generate_avg_data("Tracker", 1)
    fitted = exp_fit_avg("energidata", num_data_files)
    num_res = generate_numeric_data("Tracker", 1)
    
    #Then call plotting subroutines as needed
    plot_normal_force(num_res)
    plot_friction_force(num_res)
    
    #plt.plot(avg_data[0], avg_data[2])
    #plt.plot(num_res[0], num_res[2])

    #Finally, show the plots
    plt.show()

main()
