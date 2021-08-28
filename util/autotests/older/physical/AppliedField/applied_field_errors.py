#!/usr/bin/env python3

import numpy as numpy

lam   = 1.0            # damping
gamma = 1.76e11        # gyromagnetic ratio
H     = 10.0           # applied field

def sech(x):
    return 1.0 / numpy.cosh(x)

def S(t):
    """
    returns the analytic result for a single spin initially along the x-axis
    with an applied field along the z-axis at time t
    """
    global lam, gamma, H

    arg = t * gamma * H / (1 + lam**2)
    larg = lam * arg

    Sx = sech(larg) * numpy.cos(arg)
    Sy = sech(larg) * numpy.sin(arg)
    Sz = numpy.tanh(larg)
    return numpy.array([Sx, Sy, Sz])

def parse_input():
    global H
    vinput = open("input")
    lines = vinput.readlines()
    for line in lines:
        if line[:1] != "#" and "=" in line:
            key, value = line.split("=")
            key = key.strip()
            if key == "sim:applied-field-strength":
                H = float(value.split()[0])
                return

parse_input()

voutput = open("output")
lines = voutput.readlines()
max_error = 0.0

for line in lines:
    if line[:1] != "#" and len(line) > 1:
        data = [float(i) for i in line.split()]
        time = data[0]
        spin = data[1:4]
        analytic_spin = S(time)
        errors = [abs(a-b) for a,b in zip(spin, analytic_spin)]
        print(time, errors[0], errors[1], errors[2])

        max_error = max(max_error, max(errors))

print("# maximum error = ", max_error)
