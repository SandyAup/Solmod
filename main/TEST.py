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

x = np.arange(10)
y = np.sin(x)

#cs1 = inter.CubicSpline(x, y, bc_type="natural")
cs1 = inter.CubicSpline(x, y, bc_type=((2, 0.0), (2, 0.0)))
cs2 = inter.CubicSpline(x, y, bc_type="clamped")
cs3 = inter.CubicSpline(x, y, bc_type="not-a-knot")


#cs4 = CubicSpline(x, y, bc_type="periodic")
xs = np.arange(-2.5, 12., 0.1)

value = cs1(xs[0])
#value = CubicSpline.__call__(cs1, xs[0])
print value

#-----------------------------------------------

result_test = inter.spline(x, y, xs, order=3, kind='smoothest')

result = []

for i in range(0, len(xs)):

	if xs[i] > x[0] and xs[i] < x[9] :
		result.append(inter.spline(x, y, xs[i], order=3, kind='smoothest'))

	elif xs[i] < x[0] :
		x_tmp = x[0] + 0.005
		y_tmp = inter.spline(x, y, x_tmp, order=3, kind='smoothest')
		B = (y_tmp - y[0]) / (x_tmp - x[0])
		logA = y[0] - (B * x[0])
		result.append(logA + B * xs[i])

	elif xs[i] > x[9] : 
		x_tmp = x[9] - 0.005
		y_tmp = inter.spline(x, y, x_tmp, order=3, kind='smoothest')
		B = (y[9] - y_tmp) / (x[9] - x_tmp)
		logA = y[9] - (B * x[9])
		result.append(logA + B * xs[i])


'''
	# Extrapolation loi de puissance : y = A x^B
	elif ekn_log < xk_log_min :
		ekn_tmp = np.log10(R_to_Ekn(1.005, a_cr, mGeV, z_cr))
		y_tmp = inter.spline(xk_log, parIS, ekn_tmp, order=3, kind='smoothest')
		B = (y_tmp - parIS[0]) / (ekn_tmp - xk_log[0])
		logA = parIS[0] - (B * xk_log[0])
		result = logA + B * ekn_log

	elif ekn_log > xk_log_max :
		ekn_tmp = np.log10(R_to_Ekn(799.995, a_cr, mGeV, z_cr))
		y_tmp = inter.spline(xk_log, parIS, ekn_tmp, order=3, kind='smoothest')
		B = (parIS[5] - y_tmp) / (xk_log[5] - ekn_tmp)
		logA = parIS[5] - (B * xk_log[5])
		result = logA + B * ekn_log
'''
#-----------------------------------------------



plt.figure(figsize=(6.5, 4))
plt.plot(x, y, 'o', label='data')


#plt.plot(xs, np.sin(xs), label='true')
plt.plot(xs, cs1(xs), label="natural")
#plt.plot(xs, cs2(xs), label="clamped")
plt.plot(xs, cs3(xs), label="not-a-knot")
plt.plot(xs, result, label="spline1D")
#plt.plot(xs, cs4(xs), label="S")
#plt.plot(xs, cs(xs, 1), label="S'")
#plt.plot(xs, cs(xs, 2), label="S''")
#plt.plot(xs, cs(xs, 3), label="S'''")
plt.xlim(-2.5, 12.)
plt.legend(loc='lower left', ncol=2)
plt.show()
