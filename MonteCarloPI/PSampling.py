"""
Estimating Pi With Monte Carlo Sampling:

Using Monte-Carlo random sampling to estimate pi by taking the ratio of
points generated inside and outside of a circle nested in a sphere

Written by Andrew Projansky: 8/22/2022

***
Need to figure out best way to slide between multiple graphs of different depths
***
"""
import matplotlib.pyplot as plt
import random
import numpy as np

"""
monte_run: for number of sampling points, randomly generates position and
checks if inside or outside a circle
----------
Inputs:
    depth: int
        integer count of the number of random points to be generated

Outputs:
    xpos_I: list
        list of x values for points inside circle
    ypos_I: list
        list of y values for points inside circle, corresponding to xpos_I
    xpos_O: list
        list of x values for points outside circle
    ypos_O: list
        list of y values for points outside circle, corresponding to xpos_O
    pi: float
        approximation of pi based on 4 times the ration of points inside circle
        versus total points 


"""
def monte_run(depth):

    xpos_I = []
    ypos_I = []
    xpos_O = []
    ypos_O = []
    for i in range(depth):
        point_gen = [random.random(), random.random()]
        dist = np.sqrt((0.5-point_gen[0])**2 + (0.5-point_gen[1])**2)
        if dist < 0.5:
            xpos_I.append(point_gen[0])
            ypos_I.append(point_gen[1])
        else:
            xpos_O.append(point_gen[0])
            ypos_O.append(point_gen[1])
    pi = 4*(len(xpos_I)/(len(xpos_I)+len(xpos_O)))
    print(pi)
    return xpos_I, ypos_I, xpos_O, ypos_O, pi

def display(xpos_I, ypos_I, xpos_O, ypos_O, pi, depth):

    circle = plt.Circle((0.5, 0.5), 0.5, alpha=0.4, color = 'green')
    square = plt.Rectangle((0,0), 1, 1, alpha=0.1, color = 'gray')
    fig, ax = plt.subplots()

    ax.add_patch(circle)
    ax.add_patch(square)
    ax.scatter(xpos_I, ypos_I, color='r')
    ax.scatter(xpos_O, ypos_O, color='b')
    ax.set_aspect('equal', adjustable='box')

    plt.title('Pi = ' + str(pi))
    plt.show()

###
depth = 1000
xpos_I, ypos_I, xpos_O, ypos_O, pi = monte_run(depth)
display(xpos_I, ypos_I, xpos_O, ypos_O, pi, depth)
