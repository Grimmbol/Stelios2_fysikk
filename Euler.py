import math
from iptrack import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#The ball mass is a guess, I don't have my notes here at the time of writing
mass_ball = 0.010
grav_constant = 9.8665

datapath = "full_data"

def polynomial(coefs, x):
    result = 0
    
    exp_counter = len(coefs)-1
    print(exp_counter, len(coefs))
    for i in range(len(coefs)):
        result += coefs[i] * (x**exp_counter)
        exp_counter -= 1

    return result

#Take the vector returned by the trvalues function as first parameter, x as the second
def euler(vals, x):
    pass

#The exponential function we are trying to fit.
#NOTE: When using the tracker data, we get a fit of the maximum ball height,
#not the potential energy
def exp_func(t,C,lam):
    return C*np.exp(lam*t)

def potential_energy(h):
    return mass_ball * grav_constant * h

#Take a folder, and a number of files to read. Expected naming: data_n
def exp_fit_avg(folder, number):
    C = 0
    lam = 0
    cov_c = 0
    cov_lam = 0
    for i in range(1, number+1):
        data = np.loadtxt("./" + str(folder)+"/data_"+str(i),skiprows=1)
        ydata = []
        tdata = [] 

        for line in data:
            ydata.append(line[1])
            tdata.append(line[0])
            
        result = curve_fit(exp_func, xdata=tdata, ydata=ydata)
        print(result)
        C += result[0][0]
        lam += result[0][1]

    #Taking average of params
    C = C/number
    lam = lam/number

    print("Found exponential based on "+ str(number) + " datasets: ")
    print(str(C) + "e^" + str(lam) + "x")
    
    #Neatly package our results
    #First the parameters themselves
    result_vec = []
    result_vec.append(C)
    result_vec.append(lam)
    
    #Then the covariances
    result_vec.append

    return result_vec
       

def main():

    degree = 15

    """
    coefs = iptrack.iptrack("full_data", degree)
    print(coefs)
    values = iptrack.trvalues(coefs, -0.5);
    print(values)
    """
    fitted = exp_fit_avg("energidata", 35)

    exp_count = degree
    """
    for i in range(degree+1):
        print(str(coefs[i])+"*x^("+str(exp_count)+")+")
        exp_count -= 1
    """
    #Plotting
    """
    x = np.linspace(-0.5,0.5, 500)
    y = polynomial(coefs, x)

    plt.plot(x, y, label = "plotkin")
    plt.show()
    """

    x1 = np.linspace(0, 10, 500)
    y1 = exp_func(x1, fitted[0], fitted[1])

    points = np.loadtxt("energidata/data_1", skiprows=1)
    pointsx = []
    pointsy = []

    for point in points:
        pointsx.append(point[0])
        pointsy.append(point[1])

    plt.plot(x1, y1, label="regression")
    plt.plot(pointsx, pointsy)

    plt.show()

    

main()
