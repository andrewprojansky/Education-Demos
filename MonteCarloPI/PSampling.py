"""
Estimating Pi With Monte Carlo Sampling:

Using Monte-Carlo random sampling to estimate pi by taking the ratio of
points generated inside and outside of a circle nested in a sphere

Written by Andrew Projansky: 8/22/2022

***
Need to update this so positions and pi are list, can slide between depths
***
"""
import matplotlib as plt
import random

def monte_run(depth):

    xpos_I = []
    ypos_I = []
    xpos_O = []
    ypos_O = []
    for i in range(depth):
        point_gen = [random.rand(), random.rand()]
        dist = (0.5-point_gen[0])**2 + (0.5-point_gen[1])**2
        if dist < 0.5:
            xpos_I.append(point_gen[0])
            ypos_I.append(poiont_gen[1])
        else:
            xpos_O.append(point_gen[0])
            ypos_O.append(poiont_gen[1])
    pi = 4*(len(xpos_I)/(len(xpos_I)+len(xpos_O)))
    return xpos_I, ypos_I, xpos_O, ypos_O, pi

def display(xpos_I, ypos_I, xpos_O, ypos_O, pi, depth):

    circle = plt.Circle((0.5, 0.5), 0.5)
    square = plt.Rectangle((0,0), 1)
    fig, ax = plt.subplots()

    ax.add_patch(circle)
    ax.add_patch(square)
    ax.scatter(xpos_I, ypos_I, color='b')
    ax.scatter(xpos_O, ypos_O, color='b')

    plt.show()
%%%
depth = 100
xpos_I, ypos_I, xpos_O, ypos_O, pi = monte_run(depth)

display(xpos_I, ypos_I, xpos_O, ypos_O, pi, depth)
