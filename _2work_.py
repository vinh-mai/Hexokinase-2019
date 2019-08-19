#======================#
#		MAIN FILE	   #
#======================#

# Import some needed packages

from scipy.integrate import odeint
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from _model_ import *
#from conservation import *
#
#
# Define numerical solutions
# Step size
h = 0.01
#
# Solution
def solution(model, initial_cond, t0, t1, p):
    # p: parameters
    t = np.arange(t0, t1, h)
    sl = integrate.odeint(model, initial_cond, t, args=p)
    return(sl)

#
# Time grid
t0 = 0.0
t1 = 45.0
t = np.arange(t0, t1, h)

# Solutons to the full model
ys = []
yss = []

for i in range(len(y0)):
    ys.append(solution(MODEL, y0[i], t0, t1, p0))
    yss.append(solution(MODEL, y0[i], t0, t1, p0s))
    
ysa = []
for i in range(len(xl)):
    ysa.append(solution(MODEL, y0[1], t0, t1, p0l[i]))
#
#=====Conservation Test=====#
#plt.figure(figsize=(7, 5))
#plt.plot(t, glc, label='G')
#plt.plot(t, atp, label='ATP')
#plt.plot(t, Pi, label='Pi')
#plt.xlabel('Time ($s$)')
#plt.ylabel('$G6P$ concentration $mM$')


#=========================#
# Numerical solutions #
# and simplified model #

list_color = ['red', 'green', 'blue',\
			  'indigo', 'crimson', 'black', 'maroon']

#=========================#
# 
plt.figure(figsize=(7,5))
plt.plot(t, ys[0][:,3], color=list_color[0], linestyle='--',\
		 label='$[P_i]$ = ' + str(Pi[0]) + ' $mM$')
plt.plot(t, yss[0][:,3], color=list_color[0], linestyle='-.', \
		 label='$[P_i]$ = ' + str(Pi[0]) + ' $mM$')
for i in range(1, len(Pi)):
    plt.plot(t, ys[i][:,3], color=list_color[i], \
    		label='$[P_i]$ = ' + str(Pi[i]) + ' $mM$')
    plt.plot(t, yss[i][:,3], color=list_color[i], \
    		label='$[P_i]$ = ' + str(Pi[i]) + ' $mM$', \
    		linestyle='-.')
plt.xlabel('Time ($s$)')
plt.ylabel('$G6P$ concentration $mM$')
plt.legend()

#=========================#
# SA illustration #

sa_color = ['g', 'r', 'b']
sa_labels = ['$0.7k_{-9}$', '$1.0k_{-9}$', '$1.3k_{-9}$']

fig=plt.figure(figsize=(7,5))
for i in range(len(xl)):
    plt.plot(t, ysa[i][:,3], color=sa_color[i],\
    	     label=sa_labels[i])
plt.xlabel('Time ($s$)')
plt.ylabel('$G6P$ concentration $mM$')
plt.legend()

plt.show()
#=======================================================#
#=						END								#
#=======================================================#
