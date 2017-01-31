#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

from math import *
import numpy as np

from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import fmin_slsqp
import scipy.interpolate as inter

from scipy.integrate import simps

from scipy.interpolate import *
import matplotlib.pyplot as plt
from _cubic import *

x = np.arange(10)
y = np.sin(x)

cs1 = CubicSpline(x, y, bc_type="natural")
cs2 = CubicSpline(x, y, bc_type="clamped")
cs3 = CubicSpline(x, y, bc_type="not-a-knot")

#cs4 = CubicSpline(x, y, bc_type="periodic")
xs = np.arange(-0.5, 9.6, 0.1)

result = inter.spline(x, y, xs, order=3, kind='smoothest')



plt.figure(figsize=(6.5, 4))
plt.plot(x, y, 'o', label='data')


#plt.plot(xs, np.sin(xs), label='true')
plt.plot(xs, cs1(xs), label="natural")
plt.plot(xs, cs2(xs), label="clamped")
plt.plot(xs, cs3(xs), label="not-a-knot")
plt.plot(xs, result, label="spline1D")
#plt.plot(xs, cs4(xs), label="S")
#plt.plot(xs, cs(xs, 1), label="S'")
#plt.plot(xs, cs(xs, 2), label="S''")
#plt.plot(xs, cs(xs, 3), label="S'''")
plt.xlim(-0.5, 9.5)
plt.legend(loc='lower left', ncol=2)
plt.show()
