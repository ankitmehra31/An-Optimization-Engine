import math
import random
import matplotlib as mp1
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as smp
import time
from scipy.misc import derivative
from IPython.display import clear_output
#User input
print("Format of Equation \n For Polynomial - a*x^n+b*x^(n-1)+c   \n For trignometric - sin(a*x^2)+cos(b*x^4) \n For logarithmic - ln(ax)-ln(b*x^2) \n For Exponential - e^(a*x^2)-e^(b*x) \n where 'a','b','c' are constants and n, n-1 are powers \n NOTE-Use proper brackets for the priority")

#---------------------------------------------------------------------------------------------------------------------------
try:
    equation = input("Enter the Equation for evaluation:")
    limit_ask=int(input("Do you want to go with the\n 1.Default limit(0,1) \n 2.You want to input limit"))
    if limit_ask==1:  #default limit
        xlimit1=0
        xlimit2=1
    elif limit_ask==2: #user generated limit
        xlimit1 = float(input("Enter 1st limit")) #Minimum value of x
        xlimit2 = float(input("Enter 2nd limit")) #Maximum value of x
    else:
        print("Invalid input")

    # Parsing
    e = 2.71828
    equation = equation.replace('^', '**')
    equation = equation.replace('sin(', 'np.sin(')
    equation = equation.replace('cos(', 'np.cos(')
    equation = equation.replace('tan(', 'np.tan(')
    equation = equation.replace('cot(', 'np.cot(')
    equation = equation.replace('sec(', 'np.sec(')
    equation = equation.replace('cosec(', 'np.cosec(')
    equation = equation.replace('ln(', 'np.log(')
    equation = equation.replace('exp(', 'np.exp(')
    equation = equation.replace('e^(', 'np.exp(')

    #solving the equation
    def func(x):
        x = x
        return eval(equation)

    #plotting tne graph for the equation
    x_axis = np.linspace(xlimit1, xlimit2, 1000)
    y_axis = []
    for i in x_axis:
        y_axis.append(func(i))
    plt.plot(x_axis, y_axis, [xlimit1, xlimit2], [0, 0], 'k')
    # plt.plot([-10 , 10] , [0 , 0] , 'k')
    plt.title("Graph of f(x) wrt to x")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid()
    plt.show()
except SyntaxError:    #handling Syntax error
    print("Please refer to the format in which the equation needs to be input")
    quit()

#------------------------------------------------------------------------------------------------------------------
def bisect(a, b):     #Bisection method algorithm
    if (func(a) * func(b) >= 0):
        print("You have not assumed right a and b\n")
        return

    c = a
    n=0
    while ((b - a) >= 0.01):
        n+=1
        # Find middle point
        c = (a + b) / 2

        # Check if middle point is root
        if (func(c) == 0.0):
            break

        # Decide the side to repeat the steps
        if (func(c) * func(a) < 0):
            b = c
        else:
            a = c

    print("The value of root is : ", "%.4f" % c)


def bisection():  #taking input and calling the bisect method
    for i in range(1000):
        a = float(input("First approximation root for bisection method:"))
        b = float(input("Second approximation root for bisection method:"))
        n = int(input("Number of iteration"))
        if func(a) * func(b) < 0:
            bisect(a, b)
            break
        else:
            print("Given approximate value do not bracket the root")
            print("Try again with different values")
    plt.plot([a, a], [0, func(a)])
    plt.annotate('a', xy=(a - 0.01, -0.2))

    plt.plot([b, b], [0, func(b)])
    plt.annotate('b', xy=(b - 0.01, -0.2))

    x_itterated = []
    for i in range(1, n + 1):
        plt.plot(x_axis, y_axis, [xlimit1, xlimit2], [0, 0], 'k')
        plt.plot([a, a], [0, func(a)])
        plt.annotate('a', xy=(a - 0.01, -0.2))

        plt.plot([b, b], [0, func(b)])
        plt.annotate('b', xy=(b - 0.01, -0.2))
        x = (a + b) / 2
        plt.plot(x, func(x), 'ro', label='x')
        plt.plot([x, x], [0, func(x)], 'k')
        print("Iteration number:", i, "|" "x=", x, "|", "f(x)=", func(x))
        time.sleep(1)
        x_itterated.append(x)
        clear_output(wait=True)    #used to plot multiple iteration plots in an animation format
        plt.plot([x, x], [0, func(x)])
        plt.annotate('x', xy=(x - 0.01, -0.2))
        plt.show()

        if func(x) == 0:
            print("The approximate root=", x)
            break
        elif func(a) * func(x) > 0:
            a = x
        else:
            b = x
        time.sleep(1)

#--------------------------------------------------------------------------------------------------------------
#Derivative done using The first principle and central difference is used as the error in that has degree 2
def derivFunc(x):
    d=0.0001
    num=func(x+d)-func(x-d)
    den=2*d
    der=num/den
    return der


def newtonRaphson(x):
    h = func(x) / derivFunc(x)
    while abs(h) >= 0.0001:
        h = func(x) / derivFunc(x)

        # x(i+1) = x(i) - f(x) / f'(x)
        x = x - h

    print("The value of the root is : ",
          "%.4f" % x)

#------------------------------------------------------------------------------------------------------------------
#Functions to be used in the Golden section Method
def check_pos(x1,x2): #checking the position
    if x2<x1:
        label='right'
    else:
        label=''
    return label
def update_interior(xl,xu):  #updating the value of upper and lower limit after each iteration
    d=((np.sqrt(5)-1)/2)*(xu-xl)
    x1=xl+d
    x2=xu-d
    return x1,x2

#FINDING MAXIMUM FUNCTION
def find_max(xl,xu,x1,x2,label):
    fx1=func(x1)
    fx2=func(x2)
    if fx2>fx1 and label=='right':
        xl=xl
        xu=x1
        new_x=update_interior(xl,xu)
        x1=new_x[0]
        x2=new_x[1]
        xopt=x2
    else:
        xl=x2
        xu=xu
        new_x=update_interior(xl,xu)
        x1=new_x[0]
        x2=new_x[1]
        xopt=x1
    return xl,xu,xopt

#FINDING MINIMUM FUNCTION
def find_min(xl,xu,x1,x2,label):
    fx1=func(x1)
    fx2=func(x2)
    if fx2>fx1 and label=='right':
        xl=x2
        xu=xu
        new_x=update_interior(xl,xu)
        x1=new_x[0]
        x2=new_x[1]
        xopt=x1
    else:
        xl=xl
        xu=x1
        new_x=update_interior(xl,xu)
        x1=new_x[0]
        x2=new_x[1]
        xopt=x2
    return xl,xu,xopt


# PLOTTING FUNCTION for Golden section method
def plot_graph(xl, xu, x1, x2):
    clear_output(wait=True)

    # plot sinus graph
    plt.plot(x, y)
    plt.plot([xlimit1, xlimit2], [0, 0], 'k')

    # plot x1 point
    plt.plot(x1, func(x1), 'ro', label='x1')
    plt.plot([x1, x1], [0, func(x1)], 'k')

    # plot x2 point
    plt.plot(x2, func(x2), 'bo', label='x2')
    plt.plot([x2, x2], [0, func(x2)], 'k')

    # plot xl line
    plt.plot([xl, xl], [0, func(xl)])
    plt.annotate('xl', xy=(xl - 0.01, -0.2))

    # plot xu line
    plt.plot([xu, xu], [0, func(xu)])
    plt.annotate('xu', xy=(xu - 0.01, -0.2))

    # plot x1 line
    plt.plot([x1, x1], [0, func(x1)], 'k')
    plt.annotate('x1', xy=(x1 - 0.01, -0.2))

    # plot x2 line
    plt.plot([x2, x2], [0, func(x2)], 'k')
    plt.annotate('x2', xy=(x2 - 0.01, -0.2))

    # y-axis limit
    plt.ylim([-10, 10])
    plt.show()

#main function for the Golden section method
def golden_search(xl, xu,mode, et):
    it = 0
    e = 1
    while e >= et:
        new_x = update_interior(xl, xu)
        x1 = new_x[0]
        x2 = new_x[1]
        fx1 = func(x1)
        fx2 = func(x2)
        label = check_pos(x1, x2)
        clear_output(wait=True)
        plot_graph(xl, xu, x1, x2)  # PLOTTING
        plt.show()

        #SELECTING AND UPDATING BOUNDARY-INTERIOR POINTS
        if mode=='max':
            new_boundary = find_max(xl, xu, x1, x2, label)
        elif mode=='min':
            new_boundary=find_min(xl,xu,x1,x2,label)
        else:
            print('Please define min/max mode')
            break #exit if mode not min or max
        xl = new_boundary[0]
        xu = new_boundary[1]
        xopt = new_boundary[2]

        it += 1
        print('Iteration: ', it)
        r = (np.sqrt(5) - 1) / 2  # GOLDEN RATIO
        e = ((1 - r) * (abs((xu - xl) / xopt))) * 100  # Error
        print('Error:', e)
        time.sleep(1)


print("What operation do you want to perform \n 1.Find maximum or minimum value \n 2.Find root")
n=int(input())
if n==1:
    x = np.linspace(xlimit1, xlimit2, 100)
    y = func(x)
    # EXECUTING GOLDEN SEARCH FUNCTION
    m=int(input("Enter what you want to find 1.Max or 2.Min"))
    if m==1:
        golden_search(-10, 10,"max", 0.05)
    elif m==2:
        golden_search(-10, 10,"min", 0.05)
    else:
        print("Invalid Input")
elif n==2:
    k=int(input("Enter method \n 1.Newton Raphson \n 2.Bisection method"))
    if k==1:
        x0 = int(input("Enter Initial value for Newton Raphson"))  # Initial guess asked from user
        newtonRaphson(x0)
    elif k==2:
        bisection()
    else:
        print("Invalid Input")
else:
    print("Invalid Input")
