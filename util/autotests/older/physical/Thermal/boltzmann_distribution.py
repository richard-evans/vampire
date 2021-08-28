#!/usr/bin/env python3

import numpy
import matplotlib.pyplot as plt

k_u = 1.0e-23
T = 1

def Prob(angle):
    global k_u, T

    k_B = 1.38064852e-13
    sintheta = numpy.sin(angle)

    return sintheta * numpy.exp(-k_u*sintheta**2/(k_B*T))

def z_angle(vector):
    unit_vector = vector / numpy.linalg.norm(vector[0:3])
    return numpy.arccos(numpy.dot(unit_vector[0:3], (0, 0, 1)))

def normalise_height(array):
    inv_max_height = 1.0 / numpy.amax(array)
    array[:] = [val*inv_max_height for val in array]

def get_bins():
    voutput = open("output")
    angle_bins = numpy.zeros(181)

    line = voutput.readline()
    while line:
        if line[:1] != "#" and len(line) > 1:
            spin = [float(i) for i in line.split()[0:3]]
            angle = int(round(z_angle(spin)*180/numpy.pi))
            angle_bins[angle] += 1
        line = voutput.readline()
    voutput.close()

    return angle_bins



distribution = [Prob(i*numpy.pi/180) for i in range(0, 181)]
normalise_height(distribution)

angle_bins = get_bins()
rev_bins = angle_bins[::-1]
angle_bins = [x+y for x,y in zip(angle_bins, rev_bins)]

normalise_height(angle_bins)


max_error = 0.0
for i in range(0, 181):
    error = abs(angle_bins[i] - distribution[i])
    if error > max_error:
        max_error = error

print(max_error)

plt.plot(angle_bins)
plt.plot(distribution)

plt.show()
