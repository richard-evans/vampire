#!/usr/bin/env python3

import numpy

k_u = 1.0e-23
T = 300.0

def Prob(angle):
    global k_u, T

    k_B = 1.38064852e-13
    sintheta = numpy.sin(angle)

    return sintheta * numpy.exp(-k_u*sintheta**2/(k_B*T))

