#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

# the program takes conformal field         theta*R
# and the conformal time derivative if it   d(theta*R)/d(eta)
# for the moment we just R=eta for Radiation and R=1 in flat space

def ICgen(R,ic='flat',sigma=0.1,bst=1):
    delta = np.linspace(0,2,1000)
    if ic == 'flat':
        theta = bst*delta/delta
        vheta = theta
        theta = theta*R
    elif ic == 'axitm':
        theta = bst*np.exp(-delta/sigma)
        vheta = theta
        theta = theta*R
    elif ic == 'axitv':
        theta = delta*0
        vheta = bst*np.exp(-delta/sigma)*R
    x = np.column_stack((delta, theta, vheta))
    np.savetxt('ics.txt', x, delimiter=' ', fmt='%.8e %.8e %.8e')   # X is an array

IcGen(1,'axitv',0.3,30)
